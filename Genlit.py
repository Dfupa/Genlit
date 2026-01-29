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
    ACRONYM_PATTERN = re.compile(r"^[A-Z0-9\-]{2,8}$")
    text = entity.data_point.text
    # keep acronyms only
    if not ACRONYM_PATTERN.match(text):
        return False
    return True


def has_strong_gene_link(entity, min_score: float = 0.75) -> bool:
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
    search_terms: str = " AND (genetics OR variant)",
    clinical_relevance: list[str] = ["Pathogenic", "Likely pathogenic"],
    retnumber: int = 10,
    min_score: float = 0.75,
) -> dict:
    """
    This function performs the search based on the disease term.
    Uses a NCBI query by default, as well as filtering by clinical relevance
    And incorporating the tagger NLP min scoring for strong evidence
    Returns a dictionary of xlinked terms
    """
    logger.info(f"Searching for: {disease_term}")

    results = search_and_fetch_pubmed(
        str(disease_term + search_terms), retmax=retnumber, fetch_pmc=True
    )

    pmids = [r["pmid"] for r in results]

    abstract_genes = []
    full_text_genes = []
    abstract_variants = []
    full_text_variants = []

    for r in results:
        text = r["abstract"]
        if not text:
            continue

        # Extract entities from abstracts
        genes = extract_entities(text, min_score)
        abstract_genes.extend(genes)

        abstract_VAR = extract_variants(text)
        abstract_variants.extend(abstract_VAR)

        if r["pmc_xml"]:
            full_text = ET.tostring(r["pmc_xml"], encoding="unicode")

            # Extract entities from full text
            full_text_g = extract_entities(full_text, min_score)
            full_text_genes.extend(full_text_g)

            full_VAR = extract_variants(full_text)
            full_text_variants.extend(full_VAR)

    # Fetch ClinVar variants
    ClinVar = search_and_fetch_clinvar(disease_term, retmax=20)

    df = pd.concat(ClinVar, ignore_index=True)

    # Filter by clinical relevance:
    fdf = df[df["clinical_significance"].isin(clinical_relevance)]

    clinvar_genes = (
        fdf["genes"]
        .dropna()
        .apply(lambda x: [g.strip() for g in x.split(",")])
        .explode()
        .unique()
    )

    clinvar_variants = fdf["variation_name"].dropna().unique().tolist()

    # Dedup all variant lists just in case
    abstract_variants = list(set(abstract_variants))
    full_text_variants = list(set(full_text_variants))
    clinvar_variants = list(set(clinvar_variants))

    return {
        "genes_in_both_abstract_clinvar": sorted(
            set(abstract_genes) & set(clinvar_genes)
        ),
        "genes_only_abstract": sorted(set(abstract_genes) - set(clinvar_genes)),
        "genes_only_clinvar": sorted(set(clinvar_genes) - set(abstract_genes)),
        "genes_only_full_text": sorted(set(full_text_genes) - set(abstract_genes)),
        "variants_in_abstract_clinvar": sorted(
            set(abstract_variants) & set(clinvar_variants)
        ),
        "variants_only_text": sorted(set(abstract_variants) - set(clinvar_variants)),
        "variants_only_clinvar": sorted(set(clinvar_variants) - set(abstract_variants)),
        "variants_only_full_text": sorted(
            set(full_text_variants) - set(abstract_variants)
        ),
        "pmids": pmids,
    }


# -------------------------
# Output formatting func.
# -------------------------


def output_results(
    disease: str,
    pipeline_results: dict,
    perplexity_verdict: dict = None,
    output_file: str = None,
    verbose: bool = False,
) -> None:
    """
    Format and output results as CSV or TXT table.

    Args:
        disease: The disease term searched
        pipeline_results: Dictionary from run_search()
        perplexity_verdict: Dictionary from validate_with_perplexity()
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

        # If no pattern match, return original (could be rsID or other format)
        return variant_str

    # 1) Genes section
    if perplexity_verdict and "error" not in perplexity_verdict:

        # Extract gene lists from Perplexity verdict for quick lookup
        confirmed_genes = set(perplexity_verdict.get("confirmed", []))
        uncertain_genes = set(perplexity_verdict.get("uncertain", []))
        rejected_genes = set(perplexity_verdict.get("rejected", []))

        def get_perplexity_status(entity_name: str, is_variant: bool = False) -> str:
            """
            Determine Perplexity status for an entity (gene or variant).

            Args:
                entity_name: Gene name or variant string
                is_variant: If True, extract gene name from variant string first

            Returns:
                Status string: "Confirmed", "Uncertain", "Rejected", or "Not assessed"
            """
            # For variants, extract the gene name first
            lookup_name = (
                extract_gene_from_variant(entity_name) if is_variant else entity_name
            )

            if lookup_name in confirmed_genes:
                return "Confirmed"
            elif lookup_name in uncertain_genes:
                return "Uncertain"
            elif lookup_name in rejected_genes:
                return "Rejected"
            else:
                return "Not assessed"

        # 1) Genes section
        for gene in pipeline_results.get("genes_in_both_abstract_clinvar", []):
            status = get_perplexity_status(gene, is_variant=False)
            rows.append(
                {
                    "Entity Type": "Gene",
                    "Entity Name": gene,
                    "Associated Gene": "",
                    "Source": "Both",
                    "Perplexity Status": status,
                    "Notes": "Found in literature (PubMed|PMC) + ClinVar",
                }
            )

        for gene in pipeline_results.get("genes_only_abstract", []):
            status = get_perplexity_status(gene, is_variant=False)
            rows.append(
                {
                    "Entity Type": "Gene",
                    "Entity Name": gene,
                    "Associated Gene": "",
                    "Source": "Abstract Only",
                    "Perplexity Status": status,
                    "Notes": "Found only in literature (abstracts) but not in ClinVar",
                }
            )

        for gene in pipeline_results.get("genes_only_clinvar", []):
            status = get_perplexity_status(gene, is_variant=False)
            rows.append(
                {
                    "Entity Type": "Gene",
                    "Entity Name": gene,
                    "Associated Gene": "",
                    "Source": "ClinVar Only",
                    "Perplexity Status": status,
                    "Notes": "Found only in ClinVar but not in fetched abstracts",
                }
            )

        for gene in pipeline_results.get("genes_only_full_text", []):
            status = get_perplexity_status(gene, is_variant=False)
            rows.append(
                {
                    "Entity Type": "Gene",
                    "Entity Name": gene,
                    "Associated Gene": "",
                    "Source": "Full Text PMC",
                    "Perplexity Status": status,
                    "Notes": "Found only in full-text articles but not in ClinVar",
                }
            )

        # 2) Variants section
        for variant in pipeline_results.get("variants_in_abstract_clinvar", []):
            variant_str = str(variant)
            gene_from_variant = extract_gene_from_variant(variant_str)
            status = get_perplexity_status(variant_str, is_variant=True)
            rows.append(
                {
                    "Entity Type": "Variant",
                    "Entity Name": variant_str,
                    "Associated Gene": gene_from_variant,
                    "Source": "Both",
                    "Perplexity Status": status,
                    "Notes": "Found in literature (PubMed|PMC) + ClinVar",
                }
            )

        for variant in pipeline_results.get("variants_only_text", []):
            variant_str = str(variant)
            gene_from_variant = extract_gene_from_variant(variant_str)
            status = get_perplexity_status(variant_str, is_variant=True)
            rows.append(
                {
                    "Entity Type": "Variant",
                    "Entity Name": variant_str,
                    "Associated Gene": gene_from_variant,
                    "Source": "Abstract Only",
                    "Perplexity Status": status,
                    "Notes": "Found only in literature (abstracts) but not in ClinVar",
                }
            )

        for variant in pipeline_results.get("variants_only_full_text", []):
            variant_str = str(variant)
            gene_from_variant = extract_gene_from_variant(variant_str)
            status = get_perplexity_status(variant_str, is_variant=True)
            rows.append(
                {
                    "Entity Type": "Variant",
                    "Entity Name": variant_str,
                    "Associated Gene": gene_from_variant,
                    "Source": "Full Text Only",
                    "Perplexity Status": status,
                    "Notes": "Found only in full-text articles but not in ClinVar",
                }
            )

        for variant in pipeline_results.get("variants_only_clinvar", []):
            variant_str = str(variant)
            gene_from_variant = extract_gene_from_variant(variant_str)
            status = get_perplexity_status(variant_str, is_variant=True)
            rows.append(
                {
                    "Entity Type": "Variant",
                    "Entity Name": variant_str,
                    "Associated Gene": gene_from_variant,
                    "Source": "ClinVar Only",
                    "Perplexity Status": status,
                    "Notes": "Found only in ClinVar but not in abstracts or full text",
                }
            )

    else:

        # Fallback when Perplexity verdict is unavailable
        for gene in pipeline_results.get("genes_in_both_abstract_clinvar", []):
            rows.append(
                {
                    "Entity Type": "Gene",
                    "Entity Name": gene,
                    "Associated Gene": "",
                    "Source": "Both",
                    "Notes": "Found in literature (PubMed|PMC) + ClinVar",
                }
            )

        for gene in pipeline_results.get("genes_only_abstract", []):
            rows.append(
                {
                    "Entity Type": "Gene",
                    "Entity Name": gene,
                    "Associated Gene": "",
                    "Source": "Abstract Only",
                    "Notes": "Found only in literature (abstracts) but not in ClinVar",
                }
            )

        for gene in pipeline_results.get("genes_only_clinvar", []):
            rows.append(
                {
                    "Entity Type": "Gene",
                    "Entity Name": gene,
                    "Associated Gene": "",
                    "Source": "ClinVar Only",
                    "Notes": "Found only in ClinVar but not in fetched abstracts",
                }
            )

        for gene in pipeline_results.get("genes_only_full_text", []):
            rows.append(
                {
                    "Entity Type": "Gene",
                    "Entity Name": gene,
                    "Associated Gene": "",
                    "Source": "Full Text PMC",
                    "Notes": "Found only in full-text articles but not in ClinVar",
                }
            )

        # 2) Variants section (without Perplexity status)
        for variant in pipeline_results.get("variants_in_abstract_clinvar", []):
            variant_str = str(variant)
            gene_from_variant = extract_gene_from_variant(variant_str)
            rows.append(
                {
                    "Entity Type": "Variant",
                    "Entity Name": variant_str,
                    "Associated Gene": gene_from_variant,
                    "Source": "Both",
                    "Notes": "Found in literature (PubMed|PMC) + ClinVar",
                }
            )

        for variant in pipeline_results.get("variants_only_text", []):
            variant_str = str(variant)
            gene_from_variant = extract_gene_from_variant(variant_str)
            rows.append(
                {
                    "Entity Type": "Variant",
                    "Entity Name": variant_str,
                    "Associated Gene": gene_from_variant,
                    "Source": "Abstract Only",
                    "Notes": "Found only in literature (abstracts) but not in ClinVar",
                }
            )

        for variant in pipeline_results.get("variants_only_full_text", []):
            variant_str = str(variant)
            gene_from_variant = extract_gene_from_variant(variant_str)
            rows.append(
                {
                    "Entity Type": "Variant",
                    "Entity Name": variant_str,
                    "Associated Gene": gene_from_variant,
                    "Source": "Full Text Only",
                    "Notes": "Found only in full-text articles but not in ClinVar",
                }
            )

        for variant in pipeline_results.get("variants_only_clinvar", []):
            variant_str = str(variant)
            gene_from_variant = extract_gene_from_variant(variant_str)
            rows.append(
                {
                    "Entity Type": "Variant",
                    "Entity Name": variant_str,
                    "Associated Gene": gene_from_variant,
                    "Source": "ClinVar Only",
                    "Notes": "Found only in ClinVar but not in abstracts or full text",
                }
            )

    # Create DataFrame
    df = pd.DataFrame(rows)

    # Prepare output
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    header_text = f"Disease: {disease}\nGenerated: {timestamp}\n{'='*80}\n"

    if output_file:
        if output_file.endswith(".csv"):
            df.to_csv(output_file, index=False)
            logger.info(f"Results saved to CSV: {output_file}")

        else:  # Default to TXT (tab-separated)
            with open(output_file, "w", encoding="utf-8") as f:
                f.write(header_text)
                f.write(df.to_string(index=False))
                f.write(f"\n{'='*80}\n")
                f.write(f"Total entities: {len(df)}\n")
                if perplexity_verdict:
                    f.write(
                        f"Perplexity verdict: {json_lib.dumps(perplexity_verdict, indent=2)}\n"
                    )

            logger.info(f"Results saved to TXT: {output_file}")

    # Also print to console
    if verbose or not output_file:
        logger.info(header_text)
        logger.info(df.to_string(index=False))
        logger.info(f"{'='*80}")
        logger.info(f"Total entities: {len(df)}")
        if perplexity_verdict and "error" not in perplexity_verdict:
            logger.info(
                f"Perplexity Verdict:\n{json_lib.dumps(perplexity_verdict, indent=2)}"
            )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="GeneLit: Extract genes and variants for a disease from literature (Pubmed/PMC) and ClinVar, using NLP from Flair with a BioNER model to extract the entities, with optional Perplexity cross-check.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Author:
Diego Fuentes
diegofupa@gmail.com
Barcelona, 2026-01-29

BASIC USAGE:
python Genlit.py -q "Disease Name" -e "youremail@gmail.com" [options]

EXAMPLES:

# Basic run (PubMed/PMC + ClinVar, no Perplexity)
python Genlit.py -q "Hidradenitis Suppurativa" -e "youremail@gmail.com"

# With Perplexity cross-check and verbose output to CSV
python Genlit.py -q "Type 2 Diabetes" -e "youremail@gmail.com" --perplexity --verbose -o diabetes_results.csv

# With custom clinical relevance filtering and output to TXT
python Genlit.py -q "Crohn's Disease" -e "youremail@gmail.com" --clinical-relevance "Pathogenic" "Likely pathogenic" "Uncertain significance" -o crohns_results.txt

# With Perplexity temperature and tokens customized
python Genlit.py -q "Rheumatoid Arthritis" -e "youremail@gmail.com" --perplexity --temperature 0.3 --max-tokens 1024 -o ra_results.csv

# With custom .env file path and a retrieval max of 50 articles
python Genlit.py -q "Type 1 Diabetes" -e "youremail@gmail.com" --perplexity --env-path /path/to/custom/.env --retmax 50 -o t1d_results.csv

# High determinism (low temperature)
python Genlit.py -q "Celiac Disease" -e "youremail@gmail.com" --perplexity --temperature 0.0 -o celiac_results.csv

# Custom model selection
python Genlit.py -q "Lupus" -e "youremail@gmail.com" --perplexity --model sonar -o lupus_results.csv

# With chunking control (cost-conscious mode)
python Genlit.py -q "Pancreatic Cancer" -e "youremail@gmail.com" --perplexity --enable-chunking False -o pancreas_cheap.csv

# With custom concurrent API calls (faster processing)
python Genlit.py -q "Breast Cancer" -e "youremail@gmail.com" --perplexity --max-concurrent 5 -o breast_fast.csv

# With custom max retries for rate limit handling
python Genlit.py -q "Lung Cancer" -e "youremail@gmail.com" --perplexity --max-retries 5 -o lung_robust.csv

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
        "-p",
        "--perplexity",
        action="store_true",
        help="(Optional) Enable Perplexity API cross-check (requires PERPLEXITY_API_KEY stated as such in a .env)",
    )

    parser.add_argument(
        "--temperature",
        type=float,
        default=0.2,
        help="""Temperature for Perplexity.
Range: 0.0 (deterministic) to 1.0 (creative).
Default: 0.2 (recommended for genetics curation)
A temperature above 0.5 will yield more creative but less reliable answers, and will be flagged""",
    )

    parser.add_argument(
        "--model",
        type=str,
        choices=["sonar", "sonar-pro"],
        default="sonar-pro",
        help="""Perplexity model to use (control mode only).
- sonar: Faster, cheaper
- sonar-pro: Better reasoning (default)""",
    )

    parser.add_argument(
        "--env-path",
        type=str,
        default=None,
        help="""Custom path to .env file containing PERPLEXITY_API_KEY.
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

    # Validate Perplexity control mode parameters
    if args.perplexity and args.temperature > 0.5:
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
    logger.info(f"Tagger minimum score (NLP): {args.score}")
    logger.info(f"NCBI Email: {args.email}")
    logger.info(f"Output file: {args.output if args.output else 'STDOUT'}")

    if args.perplexity:
        logger.info("Perplexity cross-check: ENABLED")
        logger.info(f"  - Model: {args.model}")
        logger.info(f"  - Temperature: {args.temperature}")
        logger.info(f"  - Token limit: {args.token_limit}")
        logger.info(f"  - Enable chunking: {args.enable_chunking}")
        if args.enable_chunking:
            logger.info(f"  - Max concurrent API calls: {args.max_concurrent}")
            logger.info(f"  - Max retries per chunk: {args.max_retries}")
            logger.info(f"  - Max tokens per chunk: {args.max_tokens}")
        logger.info(
            f"  - Env file: {args.env_path if args.env_path else 'Default (.env)'}"
        )
    else:
        logger.info("Perplexity cross-check: DISABLED")

    logger.info("=" * 80)

    # ===== RUN MAIN SEARCH =====
    hs_results = run_search(
        disease_term=args.query,
        search_terms=" AND (genetics OR variant)",
        clinical_relevance=args.clinical_relevance,
        retnumber=args.retmax,
        min_score=args.score,
    )

    # ===== OPTIONAL PERPLEXITY CROSS-CHECK =====
    perplexity_verdict = None

    if args.perplexity:
        logger.info("Running Perplexity cross-check validation...")

        # Load custom .env if provided
        if args.env_path:
            from dotenv import load_dotenv

            load_dotenv(args.env_path)  # Load from custom path
            logger.info(f"Loaded custom .env from: {args.env_path}")

        else:
            from dotenv import load_dotenv

            load_dotenv()  # Load from current directory .env
            logger.info("Loaded .env from current directory")

        api_key = os.getenv("PERPLEXITY_API_KEY")

        if not api_key:
            logger.error("PERPLEXITY_API_KEY not found in environment")
            exit(1)

        # Import and call with NEW parameters

        perplexity_verdict = perplexity_API_search(
            query=args.query,
            candidates=hs_results,
            api_key=api_key,
            verbose=args.verbose,
            temperature=args.temperature,
            token_chunk=args.max_tokens,
            model=args.model,
            max_retries=args.max_retries,
            max_concurrent=args.max_concurrent,
            token_limit=args.token_limit,
            enable_chunking=args.enable_chunking,
        )

    # Output results
    output_results(
        disease=args.query,
        pipeline_results=hs_results,
        perplexity_verdict=perplexity_verdict,
        output_file=args.output,
        verbose=args.verbose,
    )

    elapsed = time.time() - start_time

    logger.info("=" * 80)
    logger.info(f"GenLit completed in {elapsed:.2f} seconds!")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
