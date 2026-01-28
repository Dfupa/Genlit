#!/usr/bin/env python3
"""
Author: Diego Fuentes
Contact email: diegofupa@gmail.com
Barcelona
Date:2026-01-23

Module for using Perplexity AI Chat API to cross-validate gene/variant-disease associations.
Uses structured JSON schema output for guaranteed valid responses.
Can be executed as is to demonstrate functionality.
"""

import json
import os


# ----------------------------
# Perplexity search
# ----------------------------
def perplexity_API_search(
    query: str,
    candidates: dict,
    client: str = None,
    verbose: bool = False,
    temperature: float = 0.2,
    max_tokens: int = 2048,
    model: str = "sonar-pro",
) -> dict:
    """
    Use Perplexity Chat API with JSON schema to sanity-check a gene/variant dictionary.
    Uses chat.completions.create() with structured output for guaranteed valid JSON response.

    Args:
        query: Disease name
        candidates: Dictionary returned of results, ideally with the output format of run_search()
        verbose: Print reasoning and API responses if True
        temperature: Controls randomness (0.0 = deterministic, 1.0 = creative). Default: 0.2
        max_tokens: Maximum tokens in response. Default: 2048
        model: Model to use ("sonar" or "sonar-pro"). Default: "sonar-pro"

    Returns:
        Dictionary with 'confirmed', 'uncertain', 'rejected', and 'notes' keys
    """
    # Define JSON schema for structured output
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
                    "notes": {
                        "type": "string",
                        "description": "Detailed explanation of categorization rationale and evidence quality",
                    },
                },
                "required": ["confirmed", "uncertain", "rejected", "notes"],
                "additionalProperties": False,
            },
        },
    }

    # Build structured prompt
    prompt = f"""You are a biomedical genetics expert curating disease-gene associations for: {query}

Here are the CANDIDATE GENES and VARIANTS extracted from literature and ClinVar through a text-mining pipeline:
{json.dumps(candidates, indent=2)}

TASK:
1. Evaluate each gene and variant for genuine association with {query}
2. Use your knowledge of published literature (PubMed, PMC, ClinVar)
3. Consider: It is key to assess the biological relevance and existing literature for each candidate as cross-validation.
4. Categorize each gene as:
   - CONFIRMED: Strong peer-reviewed evidence in disease literature or ClinVar clinical assertions
   - UNCERTAIN: Limited evidence, conflicting reports, or indirect association in available literature
   - REJECTED: No credible association with the disease whatsoever, including spurious links in fetched data
5. Return structured JSON with your verdict"""

    if verbose:
        print("\n[INFO] Querying Perplexity for cross-validation...\n")
        print(f"[DEBUG] Disease: {query}")
        print(f"[DEBUG] Model: {model}")
        print(f"[DEBUG] Temperature: {temperature}, Max tokens: {max_tokens}\n")

    try:
        response = client.chat.completions.create(
            model=model,
            messages=[{"role": "user", "content": prompt}],
            temperature=temperature,
            max_tokens=max_tokens,
            response_format=response_schema,
        )

        if verbose:
            print("[LOG] ✓ API Response received\n")

        # Extract JSON from structured response
        verdict = {"confirmed": [], "uncertain": [], "rejected": [], "notes": ""}

        # Parse from response.choices[0].message.content (standard format)
        if hasattr(response, "choices") and response.choices:
            response_text = response.choices[0].message.content

            if verbose:
                print(f"[DEBUG] Raw response:\n{response_text}\n")

            try:
                # Response should be pure JSON from schema
                verdict = json.loads(response_text.strip())

                # Validate structure
                required_keys = {"confirmed", "uncertain", "rejected", "notes"}
                if not all(key in verdict for key in required_keys):
                    raise ValueError(
                        f"Missing required keys. Expected: {required_keys}, Got: {set(verdict.keys())}"
                    )

                if verbose:
                    print("[LOG] ✓ Successfully parsed structured JSON response")
                    print(f"[LOG] Confirmed: {len(verdict.get('confirmed', []))} genes")
                    print(f"[LOG] Uncertain: {len(verdict.get('uncertain', []))} genes")
                    print(f"[LOG] Rejected: {len(verdict.get('rejected', []))} genes")
                    print(f"[LOG] Notes: {verdict.get('notes', '')[:-1]}\n")

                return verdict

            except json.JSONDecodeError as e:
                print(f"[ERROR] JSON parsing error: {str(e)}")
                if verbose:
                    print(f"[DEBUG] Attempted to parse: {response_text[:200]}...")

                # Fallback: try to extract JSON from malformed response
                import re

                json_match = re.search(
                    r'\{[\s\S]*"confirmed"[\s\S]*"rejected"[\s\S]*\}', response_text
                )
                if json_match:
                    try:
                        verdict = json.loads(json_match.group(0))
                        if verbose:
                            print(
                                "[LOG] ✓ Successfully parsed JSON from fallback extraction"
                            )
                        return verdict
                    except json.JSONDecodeError:
                        pass

            except ValueError as e:
                print(f"[ERROR] chema validation error: {str(e)}")
                if verbose:
                    print(f"[DEBUG] Response was: {response_text}")

        # If structured parsing failed, return empty verdict
        print("[ERROR] Could not extract structured JSON from response")
        if verbose:
            print(f"[DEBUG] Response object: {response}")

        return verdict

    except Exception as e:
        print(f"[ERROR] Perplexity API error: {str(e)}")
        if verbose:
            import traceback

            traceback.print_exc()

        return {
            "error": str(e),
            "confirmed": [],
            "uncertain": [],
            "rejected": [],
            "notes": "",
        }


# -------------------------
# Example usage
# -------------------------

if __name__ == "__main__":

    try:
        from perplexity import Perplexity
        from dotenv import load_dotenv

        load_dotenv()
        client = Perplexity(api_key=os.getenv("PERPLEXITY_API_KEY"))

    except ImportError as e:
        print(f"[ERROR]  Perplexity API not available: {str(e)}")
        exit(1)

    query = "Hidradenitis Suppurativa"
    candidates = {
        "genes_in_both_abstract_clinvar": ["PSTPIP1"],
        "genes_only_abstract": ["GH", "IGF-1"],
        "genes_only_clinvar": ["IL36RN", "NCSTN", "PSENEN"],
        "genes_only_full_text": ["GH-R", "IGF", "MEN1"],
        "variants_in_abstract_clinvar": [],
        "variants_only_text": [],
        "variants_only_clinvar": [
            "NM_003978.5(PSTPIP1):c.1034A>G (p.Tyr345Cys)",
            "NM_003978.5(PSTPIP1):c.655C>T (p.Gln219Ter)",
            "NM_003978.5(PSTPIP1):c.831G>T (p.Glu277Asp)",
            "NM_012275.3(IL36RN):c.-28+1G>A",
            "NM_012275.3(IL36RN):c.115G>T (p.Gly39Cys)",
            "NM_012275.3(IL36RN):c.179G>T (p.Gly60Val)",
            "NM_012275.3(IL36RN):c.205_212del (p.Ser69fs)",
            "NM_012275.3(IL36RN):c.266A>G (p.Tyr89Cys)",
            "NM_012275.3(IL36RN):c.338C>A (p.Ser113Ter)",
            "NM_015331.3(NCSTN):c.1229C>T (p.Ala410Val)",
            "NM_015331.3(NCSTN):c.1285C>T (p.Arg429Ter)",
            "NM_015331.3(NCSTN):c.1352+1G>C",
            "NM_015331.3(NCSTN):c.344_351del (p.Thr115fs)",
            "NM_015331.3(NCSTN):c.97G>A (p.Gly33Arg)",
            "NM_172341.4(PSENEN):c.168T>G (p.Tyr56Ter)",
        ],
        "variants_only_full_text": [],
        "pmids": ["34630684"],
    }

    results = perplexity_API_search(query, candidates, client, verbose=True)
