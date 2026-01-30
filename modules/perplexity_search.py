#!/usr/bin/env python3

"""
Author: Diego Fuentes
Contact email: diegofupa@gmail.com
Barcelona
Date: 2026-01-29

Module for using Perplexity AI Chat API to cross-validate gene/variant-disease associations.
Uses structured JSON schema output for guaranteed valid responses.
Includes:
- Optional smart chunking for large candidate lists
- Async concurrent queries (when chunking enabled)
- Exponential backoff retry logic
- Graceful result aggregation and deduplication
- Conservative mode (no chunking) for cost-conscious users

For large datasets with hundreds/thousands of genes and variants.
Toggle chunking for cost vs speed trade-off.
"""

import json
import os
import time
import random
import asyncio
import logging
from typing import Optional, List, Dict

# Configure logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


# -------------------------
# Token Estimation Helpers
# -------------------------
def estimate_tokens(text: str) -> int:
    """
    Rough estimate of token count for a string.
    Uses: ~4 characters = 1 token as per Perplexity guidelines (and usual English text estimates)

    Args:
        text: String to estimate tokens for

    Returns:
        Estimated token count
    """
    return max(1, len(text) // 4)


def calculate_chunk_candidates_tokens(candidates_dict: Dict) -> int:
    """
    Calculate total tokens needed to serialize candidates dictionary.
    Includes JSON formatting overhead (~10%).

    Args:
        candidates_dict: Dictionary with genes/variants lists

    Returns:
        Estimated token count
    """
    json_str = json.dumps(candidates_dict, indent=2)
    base_tokens = estimate_tokens(json_str)
    # Add ~10% overhead for JSON structure just in case
    return int(base_tokens * 1.1)


def get_prompt_template_tokens(query: str) -> int:
    """
    Calculate fixed tokens in the prompt template (not including candidates).

    Args:
        query: Disease name

    Returns:
        Estimated token count for fixed prompt structure
    """

    template = f"""You are a biomedical genetics expert curating disease-gene associations.

Here are the CANDIDATE GENES and VARIANTS extracted from literature and ClinVar through a text-mining pipeline:

{{CANDIDATES_PLACEHOLDER}}

TASK:

1. Evaluate each gene and variant for genuine association with {query}

2. Use your knowledge of published literature (PubMed, PMC, ClinVar)

3. Consider: It is key to assess the biological relevance and existing literature for each candidate as cross-validation.

4. Categorize each gene as:

   - CONFIRMED: Strong peer-reviewed evidence in disease literature or ClinVar clinical assertions
   - UNCERTAIN: Limited evidence, conflicting reports, or indirect association in available literature
   - REJECTED: No credible association with the disease whatsoever, including spurious links in fetched data
   - NOT ASSESSED: If insufficient information is available to make a determination or a limit in the token count is reached.

5. Return structured JSON with your verdict"""

    return estimate_tokens(template)


# -------------------------
# Smart Chunking Logic
# -------------------------
def chunk_candidates(
    candidates: Dict,
    query: str,
    max_tokens: int = 2048,
    safety_margin: float = 0.2,
    verbose: bool = False,
) -> List[Dict]:
    """
    Split large candidate dictionary into smaller chunks based on token limits.

    Args:
        candidates: Full candidates dictionary
        query: Disease name (for prompt template estimation)
        max_tokens: Maximum tokens available for response + prompt
        safety_margin: Fraction of max_tokens to reserve (default 0.2 = 20%)
        verbose: Print chunking details

    Returns:
        List of candidate chunks, each respecting token limits
    """

    # Calculate fixed tokens in prompt template
    fixed_tokens = get_prompt_template_tokens(query)

    # Calculate available tokens for candidates
    available_tokens = int(max_tokens * (1 - safety_margin)) - fixed_tokens

    if available_tokens <= 0:
        logger.error(
            f"Insufficient token budget. Fixed prompt uses {fixed_tokens} tokens, "
            f"only {max_tokens * (1 - safety_margin)} available with safety margin."
        )
        raise ValueError("Token budget too small for prompt template")

    if verbose:
        logger.info(
            f"Token budget: {max_tokens} total, {fixed_tokens} fixed prompt, "
            f"{available_tokens} available for candidates (safety margin: {safety_margin*100}%)"
        )

    # Flatten all candidates into single list for chunking
    all_items = []

    for key, items in candidates.items():
        if isinstance(items, list):
            for item in items:
                item_tokens = estimate_tokens(json.dumps(item))
                all_items.append({"key": key, "value": item, "tokens": item_tokens})

    if not all_items:
        logger.warning("[WARNING] No candidates found in input dictionary")
        return [candidates]

    # Group items into chunks respecting token limit
    chunks = []
    current_chunk = {key: [] for key in candidates.keys()}
    current_chunk["_chunk_tokens"] = 0

    for item in all_items:
        item_tokens = (
            estimate_tokens(json.dumps(item["value"])) + 50
        )  # +50 for JSON overhead just in case

        # If adding this item exceeds limit, start new chunk
        if (
            current_chunk["_chunk_tokens"] + item_tokens > available_tokens
            and current_chunk["_chunk_tokens"] > 0
        ):
            # Remove tracking field before appending
            current_chunk.pop("_chunk_tokens")
            chunks.append(current_chunk)
            current_chunk = {key: [] for key in candidates.keys()}
            current_chunk["_chunk_tokens"] = 0

        # Add item to current chunk
        current_chunk[item["key"]].append(item["value"])
        current_chunk["_chunk_tokens"] += item_tokens

    # Add final chunk
    if current_chunk["_chunk_tokens"] > 0:
        current_chunk.pop("_chunk_tokens")
        chunks.append(current_chunk)

    # Clean empty keys
    clean_chunks = []
    for chunk in chunks:
        clean_chunk = {k: v for k, v in chunk.items() if v}
        if clean_chunk:
            clean_chunks.append(clean_chunk)

    if verbose:
        logger.info(
            f"Split {len(all_items)} candidates into {len(clean_chunks)} chunks"
        )
        for i, chunk in enumerate(clean_chunks):
            total_items = sum(len(v) for v in chunk.values() if isinstance(v, list))
            logger.info(f"Chunk {i+1}: {total_items} items")

    return clean_chunks


# -------------------------
# Async Perplexity API Calls
# -------------------------
async def resilient_perplexity_call(
    client,
    query: str,
    candidates_chunk: Dict,
    model: str = "sonar-pro",
    temperature: float = 0.2,
    max_tokens: int = 2048,
    max_retries: int = 3,
    verbose: bool = False,
) -> Optional[Dict]:
    """
    Single async API call to Perplexity with retry logic.

    Args:
        client: Perplexity client instance
        query: Disease name
        candidates_chunk: Chunk of candidates to validate
        model: Model to use
        temperature: LLM temperature
        max_tokens: Max response tokens
        max_retries: Retry attempts
        search_context_size: Context size for the search: Lower, Medium, High
        verbose: Print debug info

    Returns:
        Verdict dictionary or None if all retries failed
    """

    # Define JSON schema
    response_schema = {
        "type": "json_schema",
        "json_schema": {
            "name": "gene_variant_verdict",
            "strict": True,
            "schema": {
                "type": "object",
                "properties": {
                    "confirmed": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "List of genes/variants with strong evidence of disease association",
                    },
                    "uncertain": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "List of genes/variants with limited or conflicting evidence",
                    },
                    "rejected": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "List of genes/variants with no credible disease association whatsoever",
                    },
                    "not_assessed": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "List of genes/variants that were not assessed due to insufficient information or token limits",
                    },
                    "notes": {
                        "type": "string",
                        "description": "Detailed explanation of categorization rationale and evidence quality",
                    },
                },
                "required": [
                    "confirmed",
                    "uncertain",
                    "rejected",
                    "not_assessed",
                    "notes",
                ],
                "additionalProperties": False,
            },
        },
    }

    # Build prompt
    prompt = f"""You are a biomedical genetics expert curating disease-gene associations.

Here are the CANDIDATE GENES and VARIANTS extracted from literature and ClinVar through a text-mining pipeline:

{json.dumps(candidates_chunk, indent=2)}

TASK:

1. Evaluate each gene and variant for genuine association with {query}

2. Use your knowledge of published literature (PubMed, PMC, ClinVar)

3. Consider: It is key to assess the biological relevance and existing literature for each candidate as cross-validation.

4. Categorize each gene as:

   - CONFIRMED: Strong peer-reviewed evidence in disease literature or ClinVar clinical assertions
   - UNCERTAIN: Limited evidence, conflicting reports, or indirect association in available literature
   - REJECTED: No credible association with the disease whatsoever, including spurious links in fetched data
   - NOT ASSESSED: If insufficient information is available to make a determination or a limit in the token count is reached.

5. Return structured JSON with your verdict"""

    # Retry loop
    for attempt in range(max_retries):
        try:
            response = await client.chat.completions.create(
                model=model,
                messages=[{"role": "user", "content": prompt}],
                temperature=temperature,
                max_tokens=max_tokens,
                response_format=response_schema,
            )

            # Parse response
            response_text = response.choices[0].message.content
            verdict = json.loads(response_text.strip())

            if verbose:
                logger.info(
                    f"✓ Chunk processed: "
                    f"{len(verdict.get('confirmed', []))} confirmed, "
                    f"{len(verdict.get('rejected', []))} rejected, "
                    f"{len(verdict.get('uncertain', []))} uncertain, "
                    f"{len(verdict.get('not_assessed', []))} not assessed.\n"
                )

            return verdict

        except Exception as e:
            error_str = str(e).lower()
            is_retriable = any(
                keyword in error_str
                for keyword in [
                    "rate limit",
                    "429",
                    "timeout",
                    "connection",
                    "503",
                    "502",
                ]
            )

            if not is_retriable:
                logger.error(f"Non-retriable error: {str(e)[:100]}")
                logger.error(
                    f"Check whether it is the token limit for the candidates tokens: {calculate_chunk_candidates_tokens(candidates_chunk)} + {get_prompt_template_tokens(query)} fixed tokens for a maximum or {max_tokens}!"
                )
                return None

            if attempt < max_retries - 1:
                delay = (1.0 * (2**attempt)) + random.uniform(0, 1)
                logger.warning(
                    f"[WARNING] Attempt {attempt+1}/{max_retries} failed (retrying in {delay:.1f}s)"
                )
                await asyncio.sleep(delay)
            else:
                logger.error(f"Max retries exceeded for chunk")
                return None

    return None


# -------------------------
# Result Aggregation from Async Calls
# -------------------------
def aggregate_verdicts(verdicts: List[Dict], verbose: bool = False) -> Dict:
    """
    Merge multiple chunk verdicts into single result.
    Deduplicates and combines notes.

    Args:
        verdicts: List of verdict dictionaries from each chunk
        verbose: Print aggregation details

    Returns:
        Aggregated verdict dictionary
    """

    aggregated = {
        "confirmed": [],
        "uncertain": [],
        "rejected": [],
        "not_assessed": [],
        "notes": "",
        "chunk_count": len(verdicts),
        "successful_chunks": sum(1 for v in verdicts if v is not None),
    }

    all_notes = []

    for verdict in verdicts:
        if verdict is None:
            continue

        # Add items and deduplicate
        for key in ["confirmed", "uncertain", "rejected", "not_assessed"]:
            if key in verdict:
                aggregated[key].extend(verdict[key])

        # Collect notes
        if verdict.get("notes"):
            all_notes.append(verdict["notes"][:-1])

    # QA STEP: Deduplicate while preserving order
    for key in ["confirmed", "uncertain", "rejected", "not_assessed"]:
        seen = set()
        deduped = []
        for item in aggregated[key]:
            if item not in seen:
                seen.add(item)
                deduped.append(item)
        aggregated[key] = sorted(deduped)

    # Combine notes
    aggregated["notes"] = (
        " | ".join(all_notes) if all_notes else "Processed all chunks successfully."
    )

    if verbose:
        logger.info(
            f"Aggregated results: "
            f"{len(aggregated['confirmed'])} confirmed, "
            f"{len(aggregated['uncertain'])} uncertain, "
            f"{len(aggregated['rejected'])} rejected, "
            f"{len(aggregated['not_assessed'])} not assessed"
        )

    return aggregated


# -------------------------
# Main Async Function
# -------------------------
async def perplexity_API_search_async(
    query: str,
    candidates: Dict,
    api_key: str,
    verbose: bool = False,
    temperature: float = 0.2,
    max_tokens: int = 2048,
    model: str = "sonar-pro",
    max_retries: int = 3,
    max_concurrent: int = 3,
    token_limit: int = 4096,
    enable_chunking: bool = True,
) -> Dict:
    """
    Async version with optional chunking: Split large candidates into chunks, query concurrently.

    Args:
        query: Disease name
        candidates: Full candidates dictionary (may be large)
        api_key: Perplexity API key
        verbose: Print debug info
        temperature: LLM temperature
        max_tokens: Response token limit per chunk
        model: Model to use
        max_retries: Retry attempts per chunk
        max_concurrent: Max concurrent API calls
        token_limit: Total token budget (default 4096)
        search_context_size: Context size for the search: Lower, Medium, High
        enable_chunking: Enable smart chunking (default True).
                        If False, sends all candidates in single API call (cheaper but slower/limited by token)

    Returns:
        Aggregated verdict dictionary
    """

    try:
        from perplexity import AsyncPerplexity
    except ImportError:
        logger.error("AsyncPerplexity not available. Install: pip install perplexityai")
        return {
            "error": "AsyncPerplexity import failed",
            "confirmed": [],
            "uncertain": [],
            "rejected": [],
            "not_assessed": [],
            "notes": "",
        }

    # Determine chunking strategy
    if not enable_chunking:
        if verbose:
            logger.info(
                "[MODE] Chunking DISABLED - Conservative mode (single API call)"
            )
        chunks = [candidates]
    else:
        if verbose:
            logger.info(
                "[MODE] Chunking ENABLED - Aggressive mode (concurrent multi-chunk)"
            )

        # Chunk candidates smartly
        if verbose:
            logger.info(
                f"Chunking {sum(len(v) if isinstance(v, list) else 0 for v in candidates.values())} "
                f"candidates based on {token_limit} token limit..."
            )

        chunks = chunk_candidates(
            candidates,
            query,
            max_tokens=token_limit,
            safety_margin=0.2,
            verbose=verbose,
        )

    if verbose:
        logger.info(f"Created {len(chunks)} chunk(s) for processing")

    # Create async client and run queries
    async with AsyncPerplexity(api_key=api_key) as client:

        if enable_chunking and len(chunks) > 1:
            # Concurrent processing with semaphore
            semaphore = asyncio.Semaphore(max_concurrent)

            async def limited_call(chunk):
                async with semaphore:
                    return await resilient_perplexity_call(
                        client,
                        query,
                        chunk,
                        model=model,
                        temperature=temperature,
                        max_tokens=max_tokens,
                        max_retries=max_retries,
                        verbose=verbose,
                    )

            # Execute all chunks concurrently
            if verbose:
                logger.info(
                    f"Launching {len(chunks)} concurrent queries "
                    f"(max {max_concurrent} simultaneous)..."
                )

            verdicts = await asyncio.gather(
                *[limited_call(chunk) for chunk in chunks], return_exceptions=False
            )
        else:
            # Single chunk or disabled chunking - single API call
            if verbose:
                logger.info(f"Sending candidates in single API call...")

            verdict = await resilient_perplexity_call(
                client,
                query,
                chunks[0],
                model=model,
                temperature=temperature,
                max_tokens=max_tokens,
                max_retries=max_retries,
                verbose=verbose,
            )
            verdicts = [verdict]

    # Aggregate results
    if verbose:
        logger.info("Aggregating results...")

    aggregated = aggregate_verdicts(verdicts, verbose=verbose)

    # Add mode information
    aggregated["chunking_enabled"] = enable_chunking
    aggregated["mode"] = (
        "concurrent" if (enable_chunking and len(chunks) > 1) else "single"
    )

    return aggregated


# -------------------------
# Wrapper for Sync Code
# -------------------------
def perplexity_API_search(
    query: str,
    candidates: Dict,
    api_key: str,
    verbose: bool = False,
    temperature: float = 0.2,
    token_chunk: int = 2048,
    model: str = "sonar-pro",
    max_retries: int = 3,
    max_concurrent: int = 3,
    token_limit: int = 4096,
    enable_chunking: bool = True,
) -> Dict:
    """
    Synchronous wrapper around async chunked search.
    Use this in your existing Genlit.py code.

    Args:
        query: Disease name
        candidates: Full candidates dictionary
        api_key: Perplexity API key
        verbose: Print debug info
        temperature: LLM temperature
        token_chunk: Response token limit per chunk
        model: Model to use
        max_retries: Retry attempts per chunk
        max_concurrent: Max concurrent API calls
        token_limit: Total token budget
        search_context_size: Context size for the search: Lower, Medium, High
        enable_chunking: Enable smart chunking (default True).
                        Set to False for conservative cost mode (single API call).

    Returns:
        Aggregated verdict dictionary
    """

    return asyncio.run(
        perplexity_API_search_async(
            query=query,
            candidates=candidates,
            api_key=api_key,
            verbose=verbose,
            temperature=temperature,
            max_tokens=token_chunk,
            model=model,
            max_retries=max_retries,
            max_concurrent=max_concurrent,
            token_limit=token_limit,
            enable_chunking=enable_chunking,
        )
    )


# -------------------------
# Example Usage
# -------------------------

if __name__ == "__main__":
    import os
    from dotenv import load_dotenv

    load_dotenv()

    query = "Pancreatic cancer"
    candidates = {
        "genes_in_both_abstract_clinvar": ["ATM", "BRCA2", "TP53"],
        "genes_only_abstract": [
            "-GOV",
            "ACLY",
            "ACSS1",
            "AKT",
            "AMIGO",
            "AMIGO1-3",
            "AMIGO2",
            "ANGPTL4",
            "APC",
            "ARID1A",
            "ARIH2",
            "AURKA",
            "BRCA1",
            "CCR2",
            "CCR5",
            "CD34",
            "CD4",
            "CD8",
            "CDK",
            "CDKN2A",
            "CHEK1",
            "CILP2",
            "CREBBP",
            "CTBP2",
            "CTLA-4",
            "DIO3",
            "E2F",
            "EGFR",
            "ERK",
            "FGF",
            "FLT1",
            "FOXP3",
            "GLUT2",
            "HER2",
            "HK2",
            "HLA-DR",
            "HM2-CAR",
            "IL-2",
            "IL2",
            "JUN",
            "KMT2D",
            "KRAS",
            "KRASG12D",
            "KRASG12R",
            "MAP2K4",
            "MAPK",
            "MARCO",
            "MCU",
            "NRF2",
            "PARP",
            "PD-1",
            "PD-L1",
            "PD-L2",
            "PDX1",
            "PI3K",
            "PIK3CA",
            "PLP",
            "PPARGC1A",
            "PTEN",
            "RAS",
            "RB1",
            "SERPINE1",
            "SHP2",
            "SLC30",
            "SLC39",
            "SOD1",
            "SOX2",
            "TERT",
            "TIM-3",
            "TIMP-1",
            "VEGF",
            "VEGFA",
            "VIM",
            "WNT",
            "YY1",
            "ZIP",
            "ZIP6",
            "ZIP7",
            "ZIP8",
        ],
        "genes_only_clinvar": [
            "CHEK2",
            "FRY",
            "N4BP2L1",
            "N4BP2L2",
            "PALB2",
            "STK11",
            "VHL",
            "ZAR1L",
        ],
        "genes_only_full_text": [
            "-I",
            "-L1",
            "-L2",
            "2-B08",
            "A2AR",
            "A2B",
            "ABCA1",
            "ABCG1",
            "ACAT",
            "ACAT1",
            "ACOD1",
            "ACSL4",
            "ACSS2",
            "ACTA2",
            "ADIRF",
            "ADSL",
            "AFP",
            "AIF1",
            "ALDOA",
            "ALDOB",
            "ALDOC",
            "ALT",
            "AMBP",
            "AMIGO-2",
            "AMIGO1",
            "AMIGO3",
            "AMOT",
            "AMPK",
            "ANKRD22",
            "APMAP",
            "AR",
            "ARG1",
            "ARG2",
            "ARIH1",
            "ARIH2-OE",
            "AST",
            "ATF4",
            "ATP7A",
            "ATR",
            "AURKB",
            "B2M",
            "B7-H1",
            "B7-H3",
            "BAD",
            "BCL2",
            "BCR",
            "BECN1",
            "BET",
            "BME",
            "BME2",
            "BMP",
            "BRAF",
            "BRCA",
            "BREA2",
            "BSA",
            "CA9",
            "CAR",
            "CBP",
            "CCL13",
            "CCL18",
            "CCL19",
            "CCL2",
            "CCL21",
            "CCL3",
            "CCL4",
            "CCL5",
            "CCL8",
            "CCNE2",
            "CCR1",
            "CCR3",
            "CCR4",
            "CCR7",
            "CD103",
            "CD11",
            "CD14",
            "CD15",
            "CD155",
            "CD163",
            "CD19",
            "CD206",
            "CD226",
            "CD27",
            "CD274",
            "CD28",
            "CD3",
            "CD31",
            "CD33",
            "CD36",
            "CD38",
            "CD39",
            "CD3D",
            "CD3E",
            "CD3G",
            "CD4-",
            "CD40",
            "CD44",
            "CD47",
            "CD56",
            "CD66",
            "CD68",
            "CD73",
            "CD79A",
            "CD79B",
            "CD80",
            "CD86",
            "CD8A",
            "CDC20",
            "CDC25C",
            "CDF",
            "CDH1",
            "CDH2",
            "CDH23",
            "CDH5",
            "CDK1",
            "CDK2",
            "CEP17",
            "CFTR",
            "CHEK",
            "CHGA",
            "CHGB",
            "CHK1",
            "CHK2",
            "CK2",
            "CLDN",
            "CLDN5",
            "COL-1A1",
            "COL1A1",
            "COX2",
            "CSF1",
            "CTLA4",
            "CTRB1",
            "CXCL1",
            "CXCL10",
            "CXCL11",
            "CXCL12",
            "CXCL13",
            "CXCL14",
            "CXCL16",
            "CXCL5",
            "CXCL8",
            "CXCL9",
            "CXCR2",
            "CXCR4",
            "CXCR6",
            "DCN",
            "DEGA",
            "DMT1",
            "DNAJB1",
            "DNAM-1",
            "DNMT",
            "DNMT1",
            "DNMT3A",
            "DRP1",
            "DSG2",
            "E2F1",
            "E3BP",
            "ECM",
            "EGF",
            "ENO1",
            "ENT1",
            "EPCAM",
            "EPO",
            "ER",
            "ERBB2",
            "ETV1",
            "EZH2",
            "F4",
            "FABP5",
            "FAF2",
            "FAK",
            "FAM111B",
            "FAS",
            "FASN",
            "FATP2",
            "FBP1",
            "FBXW11",
            "FGF-10",
            "FGF2",
            "FGFR1-4",
            "FH",
            "FKBP9",
            "FLI1",
            "FN",
            "FN1",
            "FOXA1",
            "FR4",
            "GAPDH",
            "GCDH",
            "GCG",
            "GLUT1",
            "GLUT4",
            "GM-CSF",
            "GMPS",
            "GNAI",
            "GNAS",
            "GPCR",
            "GPI",
            "GPR155",
            "GPR31",
            "GPR35",
            "GPR43",
            "GPR81",
            "GPR84",
            "GPR91",
            "GRP91",
            "GSDMB",
            "GSK-3",
            "GTPSCS",
            "GZMB",
            "H-2",
            "H19",
            "H2A",
            "H3K18",
            "H3K27",
            "H3K9",
            "HA",
            "HB-EGF",
            "HDAC",
            "HDAC-6",
            "HDAC1",
            "HDAC2",
            "HDAC3",
            "HEPFAL",
            "HER",
            "HER-2",
            "HIF",
            "HIF1A",
            "HIF2",
            "HK",
            "HK1",
            "HLA-A",
            "HLADR",
            "HM1",
            "HM2",
            "HM2-CART",
            "HM3",
            "HM4",
            "HMGB1",
            "HO-1",
            "HSP90B1",
            "HSPC111",
            "HUWE1",
            "HYOU1",
            "ICAM-1",
            "IDH",
            "IDH1",
            "IDH2",
            "IDO",
            "IDO1",
            "IFN",
            "IFNG",
            "IG-L",
            "IGF-1",
            "IGF1",
            "IGFR1",
            "IGHG1",
            "IKK",
            "IL",
            "IL-10",
            "IL-12",
            "IL-15",
            "IL-4",
            "IL-4I1",
            "IL-4R",
            "IL-6",
            "IL-7",
            "IL-8",
            "IL10",
            "IL15",
            "IL4I1",
            "ILK",
            "IMPDH",
            "INOS",
            "INS",
            "IRAK3",
            "IRF3",
            "IRG1",
            "ITGA3",
            "ITGA5",
            "ITGA7",
            "ITGAV",
            "ITGB1",
            "ITGB3",
            "ITGB5",
            "JAK",
            "JAK1",
            "JAK2",
            "JMJD1A",
            "JNK",
            "KAT2A",
            "KCTD9",
            "KDM",
            "KDM5B",
            "KDM6B",
            "KDM8",
            "KEAP1",
            "KRT19",
            "KRT7",
            "L2HGDH",
            "LC3B",
            "LDAH",
            "LDH",
            "LDH-A",
            "LDHA",
            "LDHB",
            "LIV-1",
            "LIV1",
            "LIV1A",
            "LOX",
            "LSD",
            "LSD1",
            "LUAD",
            "LUM",
            "MAP1LC3B",
            "MAT2A",
            "MAVS",
            "MCL1",
            "MCT1",
            "MDA5",
            "MDM2",
            "MDM4",
            "MET",
            "METTL3",
            "MHC-APP",
            "MHC-I",
            "MHC-II",
            "MLH1",
            "MLPH",
            "MMP-13",
            "MMP10",
            "MMP2",
            "MMP9",
            "MMR",
            "MPC1",
            "MPC2",
            "MS4A1",
            "MS4A7",
            "MSH2",
            "MTAP",
            "MTHFD2",
            "MUC5AC",
            "MYC",
            "NAD",
            "NADPH",
            "NAMPT",
            "NFAT",
            "NFAT1",
            "NFKB1",
            "NICD1",
            "NLRP3",
            "NLRX1",
            "NMDA",
            "NMDAR",
            "NNT",
            "NORAD",
            "NOTCH",
            "NOTCH1",
            "NOTCH3",
            "NR4A2",
            "OX40",
            "OXCT1",
            "P2A",
            "P300",
            "P4HA1",
            "P4HA2",
            "P53",
            "PARK2",
            "PC",
            "PCBP2",
            "PCNA",
            "PD",
            "PD1",
            "PDGFRB",
            "PDIA3",
            "PDK1",
            "PERK",
            "PFKFB3",
            "PFKFB4",
            "PGAM1",
            "PGAM4",
            "PGK1",
            "PHGDH",
            "PI3",
            "PKG",
            "PKLR",
            "PKM",
            "PKM2",
            "PLPP1",
            "PLVAP",
            "PMS2",
            "PPAR",
            "PPP",
            "PRF1",
            "PRMT",
            "PRMT5",
            "PRPF19",
            "PRPF8",
            "PRSS1",
            "PSAT1",
            "PSC",
            "PTM",
            "PTX3",
            "PVR",
            "QDPR",
            "RAB27B",
            "RACK1",
            "RANTES",
            "RB",
            "RBMX",
            "RBR",
            "RC48",
            "REDD1",
            "REG1B",
            "RELA",
            "RET",
            "RGS5",
            "RIG-1",
            "RIG-I",
            "RING1",
            "RING2",
            "RNF14",
            "RNF144B",
            "RNF216",
            "RNF31",
            "RREB-1",
            "RREB1",
            "RTK",
            "SETDB1",
            "SHH",
            "SHMT2",
            "SIRT1",
            "SIRT5",
            "SLC",
            "SLC-30A9",
            "SLC11A2",
            "SLC13A3",
            "SLC2A1",
            "SLC30A",
            "SLC30A1",
            "SLC30A10",
            "SLC30A2",
            "SLC30A3",
            "SLC30A4",
            "SLC30A5",
            "SLC30A6",
            "SLC30A7",
            "SLC30A8",
            "SLC30A9",
            "SLC39A",
            "SLC39A1",
            "SLC39A10",
            "SLC39A11",
            "SLC39A12",
            "SLC39A13",
            "SLC39A14",
            "SLC39A2",
            "SLC39A4",
            "SLC39A5",
            "SLC39A6",
            "SLC39A7",
            "SLC39A8",
            "SLC39A9",
            "SLC43A2",
            "SLC6A6",
            "SLC6A8",
            "SLC7A11",
            "SLG",
            "SLPI5",
            "SMAD4",
            "SNAI1",
            "SNARE",
            "SOCS1",
            "SOD2",
            "SP1",
            "SPARC",
            "SREBP1C",
            "STAT",
            "STAT1",
            "STAT3",
            "STAT5",
            "STAT6",
            "STING",
            "SUCNR1",
            "SUMO",
            "SUV39H1",
            "T2A",
            "TAP",
            "TAZ",
            "TCR",
            "TDO2",
            "TET2",
            "TFEB",
            "TGF",
            "TGR5",
            "TIAR1",
            "TIM3",
            "TLR2",
            "TLR4",
            "TMB",
            "TMEM163",
            "TNC",
            "TNF",
            "TORC1",
            "TPI1",
            "TPM",
            "TPO",
            "TREM1",
            "TRIM14",
            "TRIM27",
            "TRIM56",
            "TSPAN8",
            "TWIST1",
            "UBE2C",
            "UKLF",
            "USP14",
            "VCAM-1",
            "VCAN",
            "VEGF-A",
            "VEGFR1-3",
            "VEGFR2",
            "VGLUT1",
            "VH",
            "VHH",
            "VISTA",
            "VP1",
            "VP2",
            "VP3",
            "WNT1",
            "WWP2",
            "YAP",
            "YAP1",
            "YTHDF2",
            "ZAP70",
            "ZEB1",
            "ZIP1",
            "ZIP1-14",
            "ZIP1-3",
            "ZIP10",
            "ZIP11",
            "ZIP12",
            "ZIP12-14",
            "ZIP13",
            "ZIP14",
            "ZIP2",
            "ZIP3",
            "ZIP4",
            "ZIP4-8",
            "ZIP5",
            "ZIP9",
            "ZNT-1",
            "ZNT1",
            "ZNT2",
            "ZO-1",
        ],
        "variants_in_abstract_clinvar": [],
        "variants_only_text": [],
        "variants_only_clinvar": [
            "GRCh37/hg19 13q13.1(chr13:32668492-33154022)x1",
            "NC_000016.9:g.(23637719_23640524)_(23641791_23646182)del",
            "NC_000022.10:g.(29090106_29091114)_(29099555_29105993)del",
            "NM_000051.4:c.342_343ins[PX241358.1:g.1_282]",
            "NM_000059.4(BRCA2):c.3617del (p.Gly1206fs)",
            "NM_000059.4(BRCA2):c.4621A>T (p.Lys1541Ter)",
            "NM_000059.4(BRCA2):c.9477_9478insTTGAC (p.Asn3160fs)",
            "NM_000455.5(STK11):c.157_158insT (p.Asp53fs)",
            "NM_000546.6(TP53):c.1168_1171dup (p.Asp391delinsAlaTer)",
            "NM_000546.6(TP53):c.322_332delinsATTCA (p.Gly108_Leu111delinsIleGln)",
            "NM_000551.4(VHL):c.340+1G>T",
            "NM_024675.4(PALB2):c.3247_3250del (p.Glu1083fs)",
            "NM_024675.4:c.3501_3502ins[PX241361.1:g.1_290]",
        ],
        "variants_only_full_text": ["p.G12C"],
        "pmids": [
            "41601623",
            "41601038",
            "41597204",
            "41596590",
            "41595468",
            "41594707",
            "41591978",
            "41591505",
            "41589692",
            "41589468",
            "41584360",
            "41583520",
            "41583513",
            "41583444",
            "41583364",
            "41582215",
            "41581946",
            "41577813",
            "41574789",
            "41574674",
            "41572548",
            "41572361",
            "41571832",
            "41571299",
            "41570824",
            "41570797",
            "41568754",
            "41567798",
            "41567498",
            "41565617",
            "41564862",
            "41564647",
            "41564373",
            "41564125",
            "41564103",
            "41563663",
            "41562186",
            "41560687",
            "41560327",
            "41558149",
            "41556816",
            "41556452",
            "41554705",
            "41554333",
            "41553811",
            "41550919",
            "41549130",
            "41545782",
            "41545714",
            "41545706",
        ],
    }
    start_time = time.time()
    print("=" * 80)
    print("GenLit Perplexity Cross-Validation Example")
    print("=" * 80)
    print(f"\nDisease: {query}")
    print(
        f"Candidates: {sum(len(v) if isinstance(v, list) else 0 for v in candidates.values())} total items"
    )
    print()

    api_key = os.getenv("PERPLEXITY_API_KEY")
    if not api_key:
        logger.error("PERPLEXITY_API_KEY not found in environment")
        exit(1)

    # Example 1: With chunking enabled (default, faster but more API calls)
    print("\n--- MODE 1: CHUNKING ENABLED (Concurrent, Faster) ---\n")
    results_chunked = perplexity_API_search(
        query=query,
        candidates=candidates,
        api_key=api_key,
        verbose=True,
        model="sonar",
        temperature=0.2,
        max_tokens=2048,
        max_retries=3,
        max_concurrent=3,
        token_limit=4096,
        enable_chunking=True,  # ← Chunking ON
    )
    elapsed1 = time.time() - start_time
    start_time2 = time.time()
    print("\n" + "=" * 80)
    print("Results with Chunking")
    print("=" * 80)
    print(
        f"Mode: {results_chunked.get('mode')} | Chunks: {results_chunked.get('chunk_count')} "
        f"| Successful: {results_chunked.get('successful_chunks')}"
    )
    print(
        f"Confirmed: {len(results_chunked.get('confirmed', []))} | "
        f"Rejected: {len(results_chunked.get('rejected', []))}"
    )
    print()
    print(results_chunked)
    print()

    # Example 2: With chunking disabled (conservative, fewer API calls but may hit token limits)
    print("\n--- MODE 2: CHUNKING DISABLED (Conservative, Cheaper) ---\n")
    results_no_chunk = perplexity_API_search(
        query=query,
        candidates=candidates,
        api_key=api_key,
        verbose=True,
        model="sonar",
        temperature=0.2,
        max_tokens=4096,
        max_retries=3,
        max_concurrent=3,
        token_limit=4096,
        enable_chunking=False,  # Chunking OFF
    )

    print("\n" + "=" * 80)
    print("Results without Chunking")
    print("=" * 80)
    print(
        f"Mode: {results_no_chunk.get('mode')} | Chunks: {results_no_chunk.get('chunk_count')} "
        f"| Successful: {results_no_chunk.get('successful_chunks')}"
    )
    print(
        f"Confirmed: {len(results_no_chunk.get('confirmed', []))} | "
        f"Rejected: {len(results_no_chunk.get('rejected', []))}"
    )
    print()
    elapsed2 = time.time() - start_time2
    print(results_no_chunk)
    print()
    print("\n" + "=" * 80)
    print("Comparison")
    print("=" * 80)
    print("Chunking ENABLED:")
    print(f"  - API calls: {results_chunked.get('chunk_count')} (potentially parallel)")
    print(f"  - Processing time: {elapsed1}")
    print(f"  - Cost: Higher (~{results_chunked.get('chunk_count')} API calls)")
    print()
    print("Chunking DISABLED:")
    print(f"  - API calls: {results_no_chunk.get('chunk_count')} (single call)")
    print(f"  - Processing time: {elapsed2}")
    print(f"  - Cost: Lower (1 API call)")
