#!/usr/bin/env python3

import os
import re
import flair
import time
import json
import logging
import argparse
import json as json_lib
import pandas as pd
from datetime import datetime
from Bio import Entrez
from xml.etree import ElementTree as ET
from modules import *

# Configure logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

start_time = time.time()

# -------------------------
# Entity extraction func.
# -------------------------


def extract_entities(text: str, min_score: float = 0.75) -> list[str]:
    """
    Return gene and protein mentions from text.
    Uses the NLP pipeline loading a BioNER model.
    Returns a list of entities
    """
    # split text into sentences
    sentences = splitter.split(text)
    # perform NER tagging
    tagger.predict(sentences)
    # perform gene linker normalization
    gene_linker.predict(sentences)
    # perform disease linker normalization
    disease_linker.predict(sentences)
    # perform species linker normalization
    species_linker.predict(sentences)

    genes = []
    for sentence in sentences:
        # iterate over entities in sentence, filter for genes (could use links, disease, species too)
        for entity in sentence.get_labels():
            if (
                entity.value == "Gene"
                and is_valid_gene(entity)
                and has_strong_gene_link(entity, min_score)
            ):
                genes.append(entity.data_point.text)

    return genes


# -------------------------
# Helper functions for NLP entity filtering.
# -------------------------


def is_valid_gene(entity) -> bool:
    """Filter gene entities to remove common false positives and non-gene terms.
    Uses regex patterns to identify valid gene symbols (e.g., all caps, 2-8 characters, no common stop words).
    Returns True if entity is a valid gene mention, False otherwise.
    """
    ACRONYM_PATTERN = re.compile(r"^[A-Z0-9\-]{2,8}$")
    text = entity.data_point.text
    # keep acronyms only
    if not ACRONYM_PATTERN.match(text):
        return False
    return True


def has_strong_gene_link(entity, min_score: float = 0.75) -> bool:
    """Check if gene entity tagged by the NLP model has a strong link to the label
    of interest with a confidence score above the threshold.
    """
    if entity.score >= min_score:
        return True
    return False


# -------------------------
# Extract variants func.
# -------------------------


def extract_variants(text: str) -> list[str]:
    """
    Extract genetic variants from text using multiple HGVS formats and dbSNP IDs.
    Returns a list of variant strings (HGVS or rsID format) found in the text.
    Returns empty list if no variants found.

    Supported formats:
    - cDNA level: c.123A>G, c.123_124del, c.123_124ins
    - Protein level: p.Ala123Val, p.M1V, p.Gln219Ter
    - dbSNP IDs: rs12345678
    - Full HGVS: NM_003978.5(PSTPIP1):c.1034A>G (p.Tyr345Cys)
    """
    variants = []

    # Pattern 1: Full HGVS format with transcript, gene, cDNA, and protein notation
    # Example: NM_003978.5(PSTPIP1):c.1034A>G (p.Tyr345Cys)
    full_hgvs = re.compile(
        r"NM_[0-9]+\.[0-9]+\([A-Z0-9\-]+\):c\.[0-9_]+[A-Za-z]*>[A-Za-z]+"
    )

    # Pattern 2: cDNA-level HGVS notation
    # Example: c.123A>G, c.123_124del, c.123_124ins
    cdna_hgvs = re.compile(
        r"c\.[0-9]+(?:_[0-9]+)?(?:[A-Za-z]+>[A-Za-z]+|(?:del|ins|dup|inv)[A-Za-z]*)"
    )

    # Pattern 3: Protein-level HGVS notation
    # Example: p.Ala123Val, p.M1V, p.Gln219Ter
    protein_hgvs = re.compile(
        r"p\.[A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2}|p\.[A-Z][0-9]+[A-Z*]"
    )

    # Pattern 4: dbSNP IDs
    # Example: rs12345678, rs999999999
    rs_id = re.compile(r"\brs[0-9]+\b")

    # Extract all variant types with priority (most specific first)

    # Full HGVS (highest priority - most informative)
    full_matches = full_hgvs.findall(text)
    if full_matches:
        variants.extend(full_matches)

    # cDNA HGVS (if not already captured in full HGVS)
    cdna_matches = cdna_hgvs.findall(text)
    for match in cdna_matches:
        if match not in variants:  # Avoid duplicates
            variants.append(match)

    # Protein HGVS (if not already captured in full HGVS)
    protein_matches = protein_hgvs.findall(text)
    for match in protein_matches:
        if match not in variants:
            variants.append(match)

    # dbSNP IDs (if not already captured)
    rs_matches = rs_id.findall(text)
    for match in rs_matches:
        if match not in variants:
            variants.append(match)

    # Remove duplicates while preserving order
    seen = set()
    unique_variants = []
    for variant in variants:
        if variant not in seen:
            seen.add(variant)
            unique_variants.append(variant)

    return unique_variants


# -------------------------
# Interlink Pubmed/PMC + Clinvar
# -------------------------


def run_search(
    disease_term: str,
    search_terms: str = " AND ((genetics OR genomics) OR (variant OR genetic variant))",
    clinical_relevance: list[str] = ["Pathogenic", "Likely pathogenic"],
    retnumber: int = 10,
    min_score: float = 0.75,
    min_review_stars: int = 0,
) -> dict:
    """
    This function performs the search based on the disease term.
    Uses a NCBI query by default, as well as filtering by clinical relevance
    And incorporating the tagger NLP min scoring for strong evidence.
    Integrates HGNC gene deduplication via verify_genes_with_hgnc()
    to verify genes against the HUGO Gene Nomenclature Committee database.

    Returns a dictionary of cross-linked terms with VERIFIED genes
    """

    logger.info(f"Searching for: {disease_term}")

    # ===== PUBMED/PMC SEARCH =====

    results = search_and_fetch_pubmed(
        str(disease_term + search_terms), retmax=retnumber, fetch_pmc=True
    )

    # Handle empty literature results
    if not results:
        logger.warning(
            f"No PubMed articles found for query: {disease_term}{search_terms}"
        )
        pmids = []
        abstract_genes = []
        full_text_genes = []
        abstract_variants = []
        full_text_variants = []
    else:
        pmids = [r["pmid"] for r in results]
        abstract_genes = []
        full_text_genes = []
        abstract_variants = []
        full_text_variants = []

        # Extract entities from abstracts and full text
        for r in results:
            text = r["abstract"]
            if not text:
                continue

            # Extract entities from abstracts
            genes = extract_entities(text, min_score)
            abstract_genes.extend(genes)

            abstract_VAR = extract_variants(text)
            abstract_variants.extend(abstract_VAR)

            # Extract from full text if available
            if r.get("pmc_xml"):
                try:
                    full_text = ET.tostring(r["pmc_xml"], encoding="unicode")

                    # Extract entities from full text
                    full_text_g = extract_entities(full_text, min_score)
                    full_text_genes.extend(full_text_g)

                    full_VAR = extract_variants(full_text)
                    full_text_variants.extend(full_VAR)

                except Exception as e:
                    logger.warning(
                        f"Error extracting full text for PMID {r['pmid']}: {e}"
                    )
                    continue

    # ===== CLINVAR SEARCH =====

    ClinVar = search_and_fetch_clinvar(disease_term, retmax=retnumber)

    # Handle empty ClinVar results
    if not ClinVar:
        logger.warning(f"No ClinVar results found for '{disease_term}'")
        clinvar_genes = []
        clinvar_variants = []
        df = pd.DataFrame()
        fdf = pd.DataFrame()
    else:
        df = pd.concat(ClinVar, ignore_index=True)

        # Filter by ClinVar review strength (0–4 stars)
        if min_review_stars is not None and min_review_stars > 0:
            before = len(df)
            df = df[df.get("review_status_stars", 0) >= min_review_stars]
            after = len(df)
            logger.info(
                f"ClinVar review-status filter: >= {min_review_stars} stars "
                f"(Filtered from {before} → {after} final variants)"
            )

        # Filter by clinical relevance
        fdf = df[df["clinical_significance"].isin(clinical_relevance)]

        # Extract genes and variants from filtered ClinVar results
        if fdf.empty:
            logger.warning(
                f"No pathogenic/likely pathogenic variants in ClinVar for '{disease_term}'. "
                f"(Found {len(df)} total variants, but none match significance filter)"
            )
            clinvar_genes = []
            clinvar_variants = []
        else:
            try:
                # Extract gene symbols
                clinvar_genes = (
                    fdf["genes"]
                    .dropna()
                    .apply(lambda x: [g.strip() for g in x.split(",")])
                    .explode()
                    .unique()
                    .tolist()
                )

            except Exception as e:
                logger.warning(f"Error extracting ClinVar genes: {e}")
                clinvar_genes = []

            try:
                # Extract variation names
                clinvar_variants = fdf["variation_name"].dropna().unique().tolist()
            except Exception as e:
                logger.warning(f"Error extracting ClinVar variants: {e}")
                clinvar_variants = []

    # ===== HGNC GENE DEDUPLICATION & VERIFICATION =====
    logger.info("=" * 80)
    logger.info("Running HGNC gene verification and deduplication...")
    logger.info("=" * 80)

    # Collect all genes from all sources (before verification)
    all_genes_raw = list(
        set(abstract_genes) | set(full_text_genes) | set(clinvar_genes)
    )
    all_variants_raw = list(
        set(abstract_variants) | set(full_text_variants) | set(clinvar_variants)
    )

    logger.info(f"Raw genes collected: {len(all_genes_raw)}")
    logger.info(f"Raw variants collected: {len(all_variants_raw)}")

    # Call wrapper function to verify genes against HGNC database
    hgnc_verification = None
    verified_genes_canonical = []
    verified_genes_unmapped = []
    verified_genes_invalid = []

    try:
        logger.info("Calling verify_genes_with_hgnc()...")
        hgnc_verification = verify_genes_with_hgnc(all_genes_raw, all_variants_raw)

        # Extract verification results
        verified_genes_canonical = hgnc_verification.get("genes_canonical", [])
        verified_genes_unmapped = hgnc_verification.get("genes_unmapped", [])
        verified_genes_invalid = hgnc_verification.get("genes_invalid", [])

        # Log detailed verification results
        logger.info("=" * 80)
        logger.info("HGNC VERIFICATION RESULTS:")
        logger.info("=" * 80)
        logger.info(f"✓ Canonical (verified): {len(verified_genes_canonical)} genes")
        if verified_genes_canonical:
            logger.info(
                f"  Genes: {', '.join(verified_genes_canonical[:10])}"
                + (
                    f"... (+{len(verified_genes_canonical)-10})"
                    if len(verified_genes_canonical) > 10
                    else ""
                )
            )

        logger.info(
            f"⚠ Unmapped (not found in HGNC): {len(verified_genes_unmapped)} genes"
        )
        if verified_genes_unmapped:
            unmapped_names = [
                g.get("canonical", g) if isinstance(g, dict) else g
                for g in verified_genes_unmapped[:10]
            ]
            logger.info(
                f"  Genes: {', '.join(unmapped_names)}"
                + (
                    f"... (+{len(verified_genes_unmapped)-10})"
                    if len(verified_genes_unmapped) > 10
                    else ""
                )
            )

        logger.info(f"✗ Invalid (filtered out): {len(verified_genes_invalid)} entries")
        logger.info("=" * 80)

    except Exception as e:
        logger.error(f"HGNC verification failed: {e}")
        logger.debug(f"Exception details: {e}", exc_info=True)
        logger.warning("Falling back to raw genes (without HGNC verification)")

        # Fallback -> use raw genes if verification fails
        verified_genes_canonical = all_genes_raw
        hgnc_verification = None

    # ===== CALCULATE SUMMARY STATISTICS =====

    all_genes = set(verified_genes_canonical)
    all_variants = (
        set(abstract_variants) | set(full_text_variants) | set(clinvar_variants)
    )

    # Check for complete empty results
    if not all_genes and not all_variants and not pmids:
        logger.error(f"No results found from any source for: {disease_term}")
        logger.info(
            "Consider trying different search terms or checking disease name spelling"
        )

    elif not all_genes and not all_variants:
        logger.warning(f"No genes or variants found, but found {len(pmids)} articles")

    # Remove duplicates from gene and variant lists (if any)

    abstract_variants = list(set(abstract_variants))
    full_text_variants = list(set(full_text_variants))
    clinvar_variants = list(set(clinvar_variants))

    # ===== COMPILE RESULTS DICTIONARY =====

    results_dict = {
        "genes_in_both_abstract_clinvar": sorted(
            set(abstract_genes) & set(verified_genes_canonical)
        ),
        "genes_only_abstract": sorted(
            set(abstract_genes) - set(verified_genes_canonical)
        ),
        "genes_only_clinvar": sorted(
            set(clinvar_genes) - set(verified_genes_canonical)
        ),
        "genes_only_full_text": sorted(set(full_text_genes) - set(abstract_genes)),
        "genes_unmapped": verified_genes_unmapped,
        "genes_invalid": verified_genes_invalid,
        "variants_in_abstract_clinvar": sorted(
            set(abstract_variants) & set(clinvar_variants)
        ),
        "variants_only_text": sorted(set(abstract_variants) - set(clinvar_variants)),
        "variants_only_clinvar": sorted(set(clinvar_variants) - set(abstract_variants)),
        "variants_only_full_text": sorted(
            set(full_text_variants) - set(abstract_variants)
        ),
        "pmids": pmids,
        "hgnc_verification": hgnc_verification,
        # Summary statistics
        "summary": {
            "disease_query": disease_term,
            "total_articles_found": len(pmids),
            "genes_from_literature_raw": len(
                set(abstract_genes) | set(full_text_genes)
            ),
            "genes_from_clinvar_raw": len(clinvar_genes),
            "total_raw_genes": len(all_genes_raw),
            "genes_verified_canonical": len(verified_genes_canonical),
            "genes_unmapped": len(verified_genes_unmapped),
            "genes_invalid": len(verified_genes_invalid),
            "total_unique_genes": len(all_genes),
            "variants_from_literature": len(
                set(abstract_variants) | set(full_text_variants)
            ),
            "variants_from_clinvar": len(clinvar_variants),
            "total_unique_variants": len(all_variants),
            "clinvar_pathogenic_variants": len(fdf) if not fdf.empty else 0,
        },
    }

    # Log summary statistics
    summary = results_dict["summary"]
    logger.info("=" * 80)
    logger.info("FINAL SEARCH SUMMARY:")
    logger.info("=" * 80)
    logger.info(f"Disease: {summary['disease_query']}")
    logger.info(f"Articles found: {summary['total_articles_found']}")
    logger.info(f"Raw genes collected: {summary['total_raw_genes']}")
    logger.info(f"Verified canonical genes: {summary['genes_verified_canonical']}")
    logger.info(f"Unmapped/novel genes: {summary['genes_unmapped']}")
    logger.info(f"Invalid genes (filtered): {summary['genes_invalid']}")
    logger.info(f"Total variants: {summary['total_unique_variants']}")
    logger.info(
        f"Pathogenic ClinVar variants: {summary['clinvar_pathogenic_variants']}"
    )
    logger.info("=" * 80)

    return results_dict


# -------------------------
# Output formatting func.
# -------------------------


def output_results(
    disease: str,
    pipeline_results: dict,
    ai_verdict: dict = None,
    output_file: str = None,
    verbose: bool = False,
) -> None:
    """
    Format and output results as CSV or TXT table.

    Args:
        disease: The disease term searched
        pipeline_results: Dictionary from run_search()
        ai_verdict: Dictionary from AI provider cross-validation search
        output_file: Output filepath (auto-detects .csv or .txt)
        verbose: Print to console if True
    """
    # Prepare data for tabular output
    rows = []

    # Internal function to extract quickly the gene name from variant string
    def extract_gene_from_variant(variant_str: str) -> str:
        """
        Extract gene name from HGVS variant string.

        Examples:
        "NM_015331.3(NCSTN):c.1229C>T (p.Ala410Val)" → "NCSTN"
        "NM_003978.5(PSTPIP1):c.1034A>G (p.Tyr345Cys)" → "PSTPIP1"
        "rs123456" → "rs123456" (rsID, no gene extraction)

        Returns:
        Gene name if found in variant string, otherwise returns the full variant string
        """
        import re

        # Pattern: NM_XXXXXX.X(GENENAME):
        match = re.search(r"\(([A-Z0-9\-]+)\):", variant_str)
        if match:
            return match.group(1)
        return variant_str  # rsID or other format

    # Helper to get AI confidence/rationale for an entity
    def get_ai_confidence(
        entity_name: str, is_variant: bool = False
    ) -> tuple[str, str]:
        """Get confidence score and rationale for entity from ai_verdict.
        For variants, tries to extract gene name first for lookup in confidence dict.
        Args:
            entity_name: Gene name or variant string
            is_variant: If True, extract gene name from variant string first for lookup
        Returns:
            Tuple[str, str] or ("", "") if not found.
        """
        if not ai_verdict or "confidence" not in ai_verdict:
            return "", ""

        # For variants, try gene name first (as status does)
        lookup_name = (
            extract_gene_from_variant(entity_name) if is_variant else entity_name
        )

        conf_str = ai_verdict["confidence"].get(lookup_name, "")
        if conf_str:
            # Split "0.95 | Multiple ClinVar assertions"
            parts = conf_str.split(" | ", 1)
            score = parts[0] if parts else ""
            rationale = parts[1] if len(parts) > 1 else ""
            return score, rationale
        return "", ""

    # 1) Genes section - WITH AI verdict
    if ai_verdict and "error" not in ai_verdict:
        # Extract gene lists from AI verdict for quick lookup
        ai_confirmed_genes = set(ai_verdict.get("confirmed", []))
        ai_uncertain_genes = set(ai_verdict.get("uncertain", []))
        ai_rejected_genes = set(ai_verdict.get("rejected", []))

        def get_ai_status(entity_name: str, is_variant: bool = False) -> str:
            """
            Determine AI status for an entity (gene or variant).

            Args:
                entity_name: Gene name or variant string
                is_variant: If True, extract gene name from variant string first

            Returns:
                Status string: "Confirmed", "Uncertain", "Rejected", or "Not assessed"
            """
            lookup_name = (
                extract_gene_from_variant(entity_name) if is_variant else entity_name
            )
            if lookup_name in ai_confirmed_genes:
                return "Confirmed"
            elif lookup_name in ai_uncertain_genes:
                return "Uncertain"
            elif lookup_name in ai_rejected_genes:
                return "Rejected"
            else:
                return "Not assessed"

        for gene in pipeline_results.get("genes_in_both_abstract_clinvar", []):
            status = get_ai_status(gene, is_variant=False)
            conf_score, conf_rationale = get_ai_confidence(
                gene, is_variant=False
            )
            rows.append(
                {
                    "Entity Type": "Gene",
                    "Entity Name": gene,
                    "Associated Gene": "",
                    "Source": "Both",
                    "AI Status": status,
                    "AI Confidence": conf_score,
                    "AI Rationale": conf_rationale,
                    "Notes": "Found in literature (PubMed|PMC) + ClinVar",
                }
            )

        #  Repeat for other gene sections (genes_only_abstract, etc.)
        for gene_list, source, note in [
            (
                pipeline_results.get("genes_only_abstract", []),
                "Abstract Only",
                "Found only in literature (abstracts) but not in ClinVar",
            ),
            (
                pipeline_results.get("genes_only_clinvar", []),
                "ClinVar Only",
                "Found only in ClinVar but not in fetched abstracts",
            ),
            (
                pipeline_results.get("genes_only_full_text", []),
                "Full Text PMC",
                "Found only in full-text articles but not in ClinVar",
            ),
        ]:
            for gene in gene_list:
                status = get_ai_status(gene, is_variant=False)
                conf_score, conf_rationale = get_ai_confidence(
                    gene, is_variant=False
                )
                rows.append(
                    {
                        "Entity Type": "Gene",
                        "Entity Name": gene,
                        "Associated Gene": "",
                        "Source": source,
                        "AI Status": status,
                        "AI Confidence": conf_score,
                        "AI Rationale": conf_rationale,
                        "Notes": note,
                    }
                )

        # 2) Variants section
        for variant_list, source, note in [
            (
                pipeline_results.get("variants_in_abstract_clinvar", []),
                "Both",
                "Found in literature (PubMed|PMC) + ClinVar",
            ),
            (
                pipeline_results.get("variants_only_text", []),
                "Abstract Only",
                "Found only in literature (abstracts) but not in ClinVar",
            ),
            (
                pipeline_results.get("variants_only_full_text", []),
                "Full Text Only",
                "Found only in full-text articles but not in ClinVar",
            ),
            (
                pipeline_results.get("variants_only_clinvar", []),
                "ClinVar Only",
                "Found only in ClinVar but not in abstracts or full text",
            ),
        ]:
            for variant in variant_list:
                variant_str = str(variant)
                gene_from_variant = extract_gene_from_variant(variant_str)
                status = get_ai_status(variant_str, is_variant=True)
                conf_score, conf_rationale = get_ai_confidence(
                    variant_str, is_variant=True
                )
                rows.append(
                    {
                        "Entity Type": "Variant",
                        "Entity Name": variant_str,
                        "Associated Gene": gene_from_variant,
                        "Source": source,
                        "AI Status": status,
                        "AI Confidence": conf_score,
                        "AI Rationale": conf_rationale,
                        "Notes": note,
                    }
                )

    else:
        # Fallback WITHOUT AI cross-check mode activated
        for gene_list, source, note in [
            (
                pipeline_results.get("genes_in_both_abstract_clinvar", []),
                "Both",
                "Found in literature (PubMed|PMC) + ClinVar",
            ),
            (
                pipeline_results.get("genes_only_abstract", []),
                "Abstract Only",
                "Found only in literature (abstracts) but not in ClinVar",
            ),
            (
                pipeline_results.get("genes_only_clinvar", []),
                "ClinVar Only",
                "Found only in ClinVar but not in fetched abstracts",
            ),
            (
                pipeline_results.get("genes_only_full_text", []),
                "Full Text PMC",
                "Found only in full-text articles but not in ClinVar",
            ),
        ]:
            for gene in gene_list:
                rows.append(
                    {
                        "Entity Type": "Gene",
                        "Entity Name": gene,
                        "Associated Gene": "",
                        "Source": source,
                        "Notes": note,
                    }
                )

        for variant_list, source, note in [
            (
                pipeline_results.get("variants_in_abstract_clinvar", []),
                "Both",
                "Found in literature (PubMed|PMC) + ClinVar",
            ),
            (
                pipeline_results.get("variants_only_text", []),
                "Abstract Only",
                "Found only in literature (abstracts) but not in ClinVar",
            ),
            (
                pipeline_results.get("variants_only_full_text", []),
                "Full Text Only",
                "Found only in full-text articles but not in ClinVar",
            ),
            (
                pipeline_results.get("variants_only_clinvar", []),
                "ClinVar Only",
                "Found only in ClinVar but not in abstracts or full text",
            ),
        ]:
            for variant in variant_list:
                variant_str = str(variant)
                gene_from_variant = extract_gene_from_variant(variant_str)
                rows.append(
                    {
                        "Entity Type": "Variant",
                        "Entity Name": variant_str,
                        "Associated Gene": gene_from_variant,
                        "Source": source,
                        "Notes": note,
                    }
                )

    df = pd.DataFrame(rows)

    # Prepare output
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    header_text = f"Disease: {disease}\nGenerated: {timestamp}\n{'='*80}\n"

    if output_file:
        if output_file.endswith(".csv"):
            df.to_csv(output_file, index=False)
            logger.info(f"Results saved to CSV: {output_file}")
        else:  # TXT → TAB-SEPARATED
            df.to_csv(output_file, sep="\t", index=False)

            with open(output_file, "r") as f:
                content = f.read()

                header_text = (
                    f"Disease: {disease}\n"
                    f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
                    f"{'='*80}\n\n"
                    f"{content}"
                )

            with open(output_file, "w", encoding="utf-8") as f:
                f.write(header_text)

            logger.info(
                f"Results saved to TSV: {output_file} (tab-separated, CSV compatible)"
            )

            # Append AI verdict if present
            if ai_verdict:
                with open(output_file, "a", encoding="utf-8") as f:
                    f.write(f"\n{'='*80}\nAI verdict:\n")
                    f.write(json_lib.dumps(ai_verdict, indent=2) + "\n")

                logger.info(
                    f"Results saved to TSV: {output_file} (tab-separated, CSV compatible)"
                )

    # Console output
    if verbose or not output_file:
        logger.info(header_text)
        logger.info(df.to_string(index=False))
        logger.info(f"{'='*80}")
        logger.info(f"Total entities: {len(df)}")
        if ai_verdict and "error" not in ai_verdict:
            logger.debug(
                f"AI Verdict:\n{json_lib.dumps(ai_verdict, indent=2)}"
            )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="GeneLit: Extract genes and variants for a disease from literature (Pubmed/PMC) and ClinVar, using NLP from Flair with a BioNER model to extract the entities, with optional AI cross-check.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Author:
Diego Fuentes
diegofupa@gmail.com
Barcelona, 2026-01-29

BASIC USAGE:
python Genlit.py -q "Disease Name" -e "youremail@gmail.com" [options]

EXAMPLES:

# Basic run (PubMed/PMC + ClinVar, no AI cross-check)
python Genlit.py -q "Hidradenitis Suppurativa" -e "youremail@gmail.com"

# With Perplexity cross-check and verbose output to CSV
python Genlit.py -q "Type 2 Diabetes" -e "youremail@gmail.com" --ai-provider perplexity --verbose -o diabetes_results.csv

# With custom clinical relevance filtering and output to TXT
python Genlit.py -q "Crohn's Disease" -e "youremail@gmail.com" --clinical-relevance "Pathogenic" "Likely pathogenic" "Uncertain significance" -o crohns_results.txt

# With Perplexity temperature and tokens customized
python Genlit.py -q "Rheumatoid Arthritis" -e "youremail@gmail.com" --ai-provider perplexity --temperature 0.3 --max-tokens 1024 -o ra_results.csv

# With custom .env file path and a retrieval max of 50 articles
python Genlit.py -q "Type 1 Diabetes" -e "youremail@gmail.com" --ai-provider perplexity --env-path /path/to/custom/.env --retmax 50 -o t1d_results.csv

# High determinism (low temperature)
python Genlit.py -q "Celiac Disease" -e "youremail@gmail.com" --ai-provider perplexity --temperature 0.0 -o celiac_results.csv

# With OpenAI cross-check
python Genlit.py -q "Type 2 Diabetes" -e "youremail@gmail.com" --ai-provider openai --model gpt-4o-mini -o diabetes_results.csv

# With Gemini cross-check
python Genlit.py -q "Breast Cancer" -e "youremail@gmail.com" --ai-provider gemini --model gemini-1.5-flash -o breast_results.csv

# With Claude cross-check
python Genlit.py -q "Lupus" -e "youremail@gmail.com" --ai-provider claude --model claude-sonnet-4-5 -o lupus_results.csv

# With Perplexity
python Genlit.py -q "Pancreatic Cancer" -e "youremail@gmail.com" --ai-provider perplexity --model sonar-pro -o pancreas_results.csv

# With chunking control (cost-conscious mode)
python Genlit.py -q "Pancreatic Cancer" -e "youremail@gmail.com" --ai-provider perplexity --enable-chunking False -o pancreas_cheap.csv

# With custom concurrent API calls (faster processing)
python Genlit.py -q "Breast Cancer" -e "youremail@gmail.com" --ai-provider openai --max-concurrent 5 -o breast_fast.csv

# With custom max retries for rate limit handling
python Genlit.py -q "Lung Cancer" -e "youremail@gmail.com" --ai-provider gemini --max-retries 5 -o lung_robust.csv

# Display help
python Genlit.py -h
python Genlit.py --help
""",
    )

    parser.add_argument(
        "-q",
        "--query",
        type=str,
        help="NCBI inspired query for the Disease term to search for (e.g., 'Hidradenitis Suppurativa')",
        required=True,
    )

    parser.add_argument(
        "-e",
        "--email",
        type=str,
        default=None,
        help="Email for NCBI Entrez API (default: None)",
        required=True,
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help="Output file (CSV or TXT). If not specified, prints to stdout.",
    )

    parser.add_argument(
        "-s",
        "--score",
        type=float,
        default=0.75,
        help="Tagging minimum score for labeling text entities as genes/variants in the NER model. It is recommended a high score to improve accuracy. (default: 0.75)",
    )

    parser.add_argument(
        "--ai-provider",
        type=str,
        choices=["perplexity", "openai", "gemini", "claude"],
        default=None,
        help="""(Optional) AI provider to use for cross-validation of gene/variant-disease associations.
If not specified, no AI cross-check is performed.
- perplexity: Uses Perplexity Sonar API (requires PERPLEXITY_API_KEY in .env)
- openai: Uses OpenAI Chat Completions API (requires OPENAI_API_KEY in .env)
- gemini: Uses Google Gemini API (requires GEMINI_API_KEY in .env)
- claude: Uses Anthropic Claude API (requires ANTHROPIC_API_KEY in .env)
Example: --ai-provider openai""",
    )

    parser.add_argument(
        "--temperature",
        type=float,
        default=0.2,
        help="""Temperature for AI provider.
Range: 0.0 (deterministic) to 1.0 (creative).
Default: 0.2 (recommended for genetics curation)
A temperature above 0.5 will yield more creative but less reliable answers, and will be flagged""",
    )

    parser.add_argument(
        "--model",
        type=str,
        default=None,
        help="""Model name to use for the selected AI provider. If not specified, uses the default for each provider.
Perplexity defaults: sonar-pro. Options: sonar, sonar-pro
OpenAI defaults: gpt-4o. Options: gpt-4o, gpt-4o-mini, gpt-4-turbo
Gemini defaults: gemini-1.5-pro. Options: gemini-1.5-pro, gemini-1.5-flash, gemini-2.0-flash
Claude defaults: claude-opus-4-5. Options: claude-opus-4-5, claude-sonnet-4-5, claude-haiku-3-5
Example: --model gpt-4o-mini""",
    )

    parser.add_argument(
        "--env-path",
        type=str,
        default=None,
        help="""Custom path to .env file containing the API key for the selected provider.
If not specified, looks for .env in current directory.
Example: /path/to/custom/.env or /home/user/configs/.env.prod""",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Verbose output (shows API calls and reasoning)",
    )

    parser.add_argument(
        "-r",
        "--retmax",
        type=int,
        default=10,
        help="Max number of PubMed articles to fetch (default: 10)",
    )

    parser.add_argument(
        "-c",
        "--clinical-relevance",
        nargs="+",
        type=str,
        default=["Pathogenic", "Likely pathogenic"],
        metavar="LEVEL",
        help="""Clinical significance levels to include in ClinVar filtering
(default: 'Pathogenic' 'Likely pathogenic').
Options: 'Pathogenic', 'Likely pathogenic', 'Uncertain significance',
'Benign', 'Likely benign'. Specify multiple values separated by spaces.""",
    )

    parser.add_argument(
        "--min-review-stars",
        type=int,
        choices=[0, 1, 2, 3, 4],
        default=0,
        help=(
            "Minimum ClinVar review status (select numerically 0–4 stars) for variants to be included.\n"
            "Stars reflect ClinVar evaluation criteria (germline & somatic):\n"
            "  ★★★★ practice guideline – submitted record with classification from a practice guideline.\n"
            "  ★★★  reviewed by expert panel – submitted record with classification by an expert panel.\n"
            "  ★★   criteria provided, multiple submitters – multiple submitters with assertion criteria\n"
            "           and evidence; for germline, classifications agree (no conflicts);\n"
            "           for somatic, multiple submitters of clinical impact.\n"
            "  ★    criteria provided – single submitter OR conflicting classifications across submitters,\n"
            "           but assertion criteria and evidence (or public contact) were provided.\n"
            "  0    no assertion criteria or no classification – one or more records without assertion\n"
            "           criteria/evidence, no classification provided, or variant only classified as part of\n"
            "           a haplotype/genotype."
        ),
    )

    parser.add_argument(
        "--enable-chunking",
        type=lambda x: x.lower() in ("true", "1", "yes", "on"),
        default=True,
        help="""Enable smart chunking for large candidate lists (default: True).
Set to False for cost-conscious mode (single API call, may hit token limits on large datasets or classify as Not assesed most of candidates).
Examples: --enable-chunking True or --enable-chunking False""",
    )

    parser.add_argument(
        "--max-concurrent",
        type=int,
        default=3,
        help="""Maximum concurrent API calls when chunking is enabled (default: 3).
Increase for faster processing, decrease if hitting rate limits.
Example: --max-concurrent 5""",
    )

    parser.add_argument(
        "--max-retries",
        type=int,
        default=3,
        help="""Maximum retry attempts per chunk for transient errors (default: 3).
Increase for better resilience to rate limits and timeouts.
Example: --max-retries 5""",
    )

    parser.add_argument(
        "--max-tokens",
        type=int,
        default=2048,
        help="""Token budget for chunking strategy (default: 2048).
Increase to process more items per chunk (fewer API calls) but check for maximum token budget.
Example: --token-limit 4096""",
    )

    parser.add_argument(
        "--token-limit",
        type=int,
        default=4096,
        help="""Token budget limit (default: 4096).
Example: --token-limit 6000""",
    )

    args = parser.parse_args()

    # Validate clinical relevance values based on ClinVar standards
    valid_levels = [
        "Pathogenic",
        "Likely pathogenic",
        "Uncertain significance",
        "Benign",
        "Likely benign",
    ]

    for level in args.clinical_relevance:
        if level not in valid_levels:
            parser.error(
                f"Invalid clinical relevance level: '{level}'. "
                f"Valid options are: {', '.join(valid_levels)}"
            )

    # Validate AI provider parameters
    if args.ai_provider and args.temperature > 0.5:
        logger.warning(
            "High temperature (>0.5) may produce inconsistent gene categorizations. "
            "Recommend using 0.0-0.3 for genetics curation."
        )

    splitter, tagger, gene_linker, disease_linker, species_linker = (
        initialize_flair_models(verbose=args.verbose)
    )

    # Store Flair splitter, tagger and linkers in globals so extract_entities() can access them
    globals()["splitter"] = splitter
    globals()["tagger"] = tagger
    globals()["gene_linker"] = gene_linker
    globals()["disease_linker"] = disease_linker
    globals()["species_linker"] = species_linker

    # Configure and run
    configure_entrez(email=args.email)

    logger.info("=" * 80)
    logger.info("GenLit Pipeline Started")
    logger.info("=" * 80)
    logger.info(f"Disease Query: {args.query}")
    logger.info(f"PubMed retmax: {args.retmax}")
    logger.info(f"Clinical relevance filter: {', '.join(args.clinical_relevance)}")
    logger.info(f"Min ClinVar review stars: {args.min_review_stars}")
    logger.info(f"Tagger minimum score (NLP): {args.score}")
    logger.info(f"NCBI Email: {args.email}")
    logger.info(f"Output file: {args.output if args.output else 'STDOUT'}")

    if args.ai_provider:
        logger.info(f"AI cross-check: ENABLED (provider: {args.ai_provider})")
        logger.info(f"  - Model: {args.model if args.model else 'default for provider'}")
        logger.info(f"  - Temperature: {args.temperature}")
        logger.info(f"  - Token limit: {args.token_limit}")
        logger.info(f"  - Enable chunking: {args.enable_chunking}")
        if args.enable_chunking:
            logger.info(f"  - Max concurrent API calls: {args.max_concurrent}")
            logger.info(f"  - Max retries per chunk: {args.max_retries}")
            logger.info(f"  - Max tokens per chunk: {args.max_tokens}")
        logger.info(f"  - Env file: {args.env_path if args.env_path else 'Default (.env)'}")
    else:
        logger.info("AI cross-check: DISABLED")

    logger.info("=" * 80)

    # ===== RUN MAIN SEARCH =====
    hs_results = run_search(
        disease_term=args.query,
        search_terms=" AND (genetics OR variant)",
        clinical_relevance=args.clinical_relevance,
        retnumber=args.retmax,
        min_score=args.score,
        min_review_stars=args.min_review_stars,
    )

    # ===== OPTIONAL AI CROSS-CHECK =====
    ai_verdict = None

    if args.ai_provider:
        logger.info(f"Running {args.ai_provider} cross-check validation...")

        if args.env_path:
            from dotenv import load_dotenv
            load_dotenv(args.env_path)
            logger.info(f"Loaded custom .env from: {args.env_path}")
        else:
            from dotenv import load_dotenv
            load_dotenv()
            logger.info("Loaded .env from current directory")

        PROVIDER_DEFAULTS = {
            "perplexity": {"key_env": "PERPLEXITY_API_KEY", "default_model": "sonar-pro"},
            "openai": {"key_env": "OPENAI_API_KEY", "default_model": "gpt-4o"},
            "gemini": {"key_env": "GEMINI_API_KEY", "default_model": "gemini-1.5-pro"},
            "claude": {"key_env": "ANTHROPIC_API_KEY", "default_model": "claude-opus-4-5"},
        }

        provider_config = PROVIDER_DEFAULTS[args.ai_provider]
        api_key = os.getenv(provider_config["key_env"])

        if not api_key:
            logger.error(f"{provider_config['key_env']} not found in environment")
            exit(1)

        selected_model = args.model if args.model else provider_config["default_model"]

        common_kwargs = dict(
            query=args.query,
            candidates=hs_results,
            api_key=api_key,
            verbose=args.verbose,
            temperature=args.temperature,
            token_chunk=args.max_tokens,
            model=selected_model,
            max_retries=args.max_retries,
            max_concurrent=args.max_concurrent,
            token_limit=args.token_limit,
            enable_chunking=args.enable_chunking,
        )

        if args.ai_provider == "perplexity":
            ai_verdict = perplexity_API_search(**common_kwargs)
        elif args.ai_provider == "openai":
            ai_verdict = openai_API_search(**common_kwargs)
        elif args.ai_provider == "gemini":
            ai_verdict = gemini_API_search(**common_kwargs)
        elif args.ai_provider == "claude":
            ai_verdict = claude_API_search(**common_kwargs)

    # Output results
    output_results(
        disease=args.query,
        pipeline_results=hs_results,
        ai_verdict=ai_verdict,
        output_file=args.output,
        verbose=args.verbose,
    )

    elapsed = time.time() - start_time

    logger.info("=" * 80)
    logger.info(f"GenLit completed in {elapsed:.2f} seconds!")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
