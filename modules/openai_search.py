#!/usr/bin/env python3

"""
Author: Diego Fuentes
Contact email: diegofupa@gmail.com
Barcelona
Date: 2026-03-31

Module for using OpenAI Chat Completions API to cross-validate gene/variant-disease associations.
Uses structured JSON output for guaranteed valid responses.
Includes:
- Optional smart chunking for large candidate lists
- Async concurrent queries (when chunking enabled)
- Exponential backoff retry logic
- Graceful result aggregation and deduplication
- Conservative mode (no chunking) for cost-conscious users

Valid models: gpt-4o, gpt-4o-mini, gpt-4-turbo
Requires: pip install openai
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
    Uses: ~4 characters = 1 token as per OpenAI guidelines (and usual English text estimates)

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

5. For EACH candidate that you mention in any category, add an entry in the 'confidence' object where:
   - The key is exactly the candidate identifier (e.g. "ATM", "BRCA2 c.3617del").
   - The value is a single string in the format:
     "<confidence_float_0.0-1.0> | <brief rationale>"
   Example:
     "ATM": "0.95 | Multiple ClinVar pathogenic assertions and strong disease literature"

    → EVERY candidate MUST have a confidence entry matching keys in arrays (confirmed, uncertain, rejected) with a confidence rationale that should be brief but informative.

6. Return structured JSON with your verdict. Strictly follow the JSON schema provided in the system instructions. NO OMISSIONS: If token-limited, use NOT_ASSESSED for remainder"""

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
# Async OpenAI API Calls
# -------------------------
async def resilient_openai_call(
    client,
    query: str,
    candidates_chunk: Dict,
    model: str = "gpt-4o",
    temperature: float = 0.2,
    max_tokens: int = 2048,
    max_retries: int = 3,
    verbose: bool = False,
) -> Optional[Dict]:
    """
    Single async API call to OpenAI with retry logic.

    Args:
        client: OpenAI async client instance
        query: Disease name
        candidates_chunk: Chunk of candidates to validate
        model: Model to use
        temperature: LLM temperature
        max_tokens: Max response tokens
        max_retries: Retry attempts
        verbose: Print debug info

    Returns:
        Verdict dictionary or None if all retries failed
    """

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

5. For EACH candidate that you mention in any category, add an entry in the 'confidence' object where:
   - The key is exactly the candidate identifier (e.g. "ATM", "BRCA2 c.3617del").
   - The value is a single string in the format:
     "<confidence_float_0.0-1.0> | <brief rationale>"
   Example:
     "ATM": "0.95 | Multiple ClinVar pathogenic assertions and strong disease literature"

    → EVERY candidate MUST have a confidence entry matching keys in arrays (confirmed, uncertain, rejected) with a confidence rationale that should be brief but informative.

6. Return structured JSON with your verdict. Strictly follow the JSON schema provided in the system instructions. NO OMISSIONS: If token-limited, use NOT_ASSESSED for remainder"""

    # Retry loop
    for attempt in range(max_retries):
        try:
            response = await client.chat.completions.create(
                model=model,
                messages=[{"role": "user", "content": prompt}],
                temperature=temperature,
                max_tokens=max_tokens,
                response_format={"type": "json_object"},
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
        "confidence": {},
    }

    all_notes = []

    for verdict in verdicts:
        if verdict is None:
            continue

        for key in ("confirmed", "uncertain", "rejected", "not_assessed"):
            if key in verdict:
                aggregated[key].extend(verdict[key])

        if verdict.get("notes"):
            all_notes.append(verdict["notes"])

        if "confidence" in verdict and isinstance(verdict["confidence"], dict):
            for cand_id, conf_str in verdict["confidence"].items():
                # last chunk's view wins; easy to reason about downstream
                aggregated["confidence"][cand_id] = conf_str

    # dedupe lists as before
    for key in ("confirmed", "uncertain", "rejected", "not_assessed"):
        seen = set()
        deduped = []
        for item in aggregated[key]:
            if item not in seen:
                seen.add(item)
                deduped.append(item)
        aggregated[key] = deduped

    aggregated["notes"] = (
        " ".join(all_notes) if all_notes else "Processed all chunks successfully."
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


def validate_confidence_coverage(
    aggregated: dict, original_candidates: dict, verbose: bool = False
) -> dict:
    """Ensure confidence covers all string candidates (genes/variants).
    If any candidate is missing from confidence, add a NOT_ASSESSED entry with rationale.
    Args:
        aggregated: Aggregated verdict dictionary with 'confidence' field
        original_candidates: Original candidates dictionary to check coverage against
        verbose: Print details about missing confidence entries
    """
    all_candidates: list[str] = []

    for key, items in original_candidates.items():
        # Only process candidate lists (genes/variants), skip metadata dicts
        if not isinstance(items, list):
            continue
        for item in items:
            # Only enforce for simple string identifiers
            if isinstance(item, str):
                all_candidates.append(item)

    conf_keys = set(aggregated.get("confidence", {}).keys())

    # Only hash strings
    missing = set(all_candidates) - conf_keys

    if missing:
        if verbose:
            logger.warning(
                f"Missing confidence for {len(missing)} candidates "
                f"(showing up to 10): {list(missing)[:10]}"
            )
        # Auto-fill NOT_ASSESSED for missing
        for cand in missing:
            aggregated["confidence"][
                cand
            ] = "0.00 | Confidence missing from model response (expected for rejected results)"

    return aggregated


# -------------------------
# Main Async Function
# -------------------------
async def openai_API_search_async(
    query: str,
    candidates: Dict,
    api_key: str,
    verbose: bool = False,
    temperature: float = 0.2,
    max_tokens: int = 2048,
    model: str = "gpt-4o",
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
        api_key: OpenAI API key
        verbose: Print debug info
        temperature: LLM temperature
        max_tokens: Response token limit per chunk
        model: Model to use
        max_retries: Retry attempts per chunk
        max_concurrent: Max concurrent API calls
        token_limit: Total token budget (default 4096)
        enable_chunking: Enable smart chunking (default True).
                        If False, sends all candidates in single API call (cheaper but slower/limited by token)

    Returns:
        Aggregated verdict dictionary
    """

    try:
        from openai import AsyncOpenAI
    except ImportError:
        logger.error("openai package not available. Install: pip install openai")
        return {
            "error": "openai import failed",
            "confirmed": [],
            "uncertain": [],
            "rejected": [],
            "not_assessed": [],
            "notes": "",
            "confidence": {},
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
    client = AsyncOpenAI(api_key=api_key)

    if enable_chunking and len(chunks) > 1:
        # Concurrent processing with semaphore
        semaphore = asyncio.Semaphore(max_concurrent)

        async def limited_call(chunk):
            async with semaphore:
                return await resilient_openai_call(
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

        verdict = await resilient_openai_call(
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
    aggregated = validate_confidence_coverage(aggregated, candidates, verbose=verbose)

    # Add mode information
    aggregated["chunking_enabled"] = enable_chunking
    aggregated["mode"] = (
        "concurrent" if (enable_chunking and len(chunks) > 1) else "single"
    )

    return aggregated


# -------------------------
# Wrapper for Sync Code
# -------------------------
def openai_API_search(
    query: str,
    candidates: Dict,
    api_key: str,
    verbose: bool = False,
    temperature: float = 0.2,
    token_chunk: int = 2048,
    model: str = "gpt-4o",
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
        api_key: OpenAI API key
        verbose: Print debug info
        temperature: LLM temperature
        token_chunk: Response token limit per chunk
        model: Model to use
        max_retries: Retry attempts per chunk
        max_concurrent: Max concurrent API calls
        token_limit: Total token budget
        enable_chunking: Enable smart chunking (default True).
                        Set to False for conservative cost mode (single API call).

    Returns:
        Aggregated verdict dictionary
    """

    return asyncio.run(
        openai_API_search_async(
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
