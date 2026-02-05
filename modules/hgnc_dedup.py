"""
Author: Diego Fuentes
Contact email: diegofupa@gmail.com
Barcelona
Date:2026-02-04

Module for HGNC database management, gene normalization, and deduplication
The main wrapper is verify_genes_with_hgnc()
"""

import json
import urllib.request
import re
import logging
from pathlib import Path

# Configure logger
logging.basicConfig(level=logging.DEBUG, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


# ============================================================================
# HGNC DATABASE FUNCTIONS
# ============================================================================


def load_hgnc_database(
    use_cache: bool = True, cache_file: str = "hgnc_cache.json"
) -> dict:
    """
    Load HGNC gene nomenclature database from official Google Cloud Storage.
    Caches locally to avoid repeated downloads.

    Parameters
    ----------
    use_cache : bool
        Use cached version if available (recommended)
    cache_file : str
        Path to cache file

    Returns
    -------
    dict
        Dictionary with gene symbol as key and {approved_symbol, aliases, previous} as value
    """

    cache_path = Path(cache_file)

    # Try to load from cache first
    if use_cache and cache_path.exists():
        logger.debug(f"Loading HGNC from cache: {cache_file}")
        with open(cache_path, "r", encoding="utf-8") as f:
            return json.load(f)

    logger.debug("Downloading HGNC database from Google Cloud Storage...")

    hgnc_url = "https://storage.googleapis.com/public-download-files/hgnc/json/json/hgnc_complete_set.json"

    try:
        with urllib.request.urlopen(hgnc_url, timeout=30) as response:
            raw_data = json.loads(response.read().decode("utf-8"))
    except Exception as e:
        logger.error(f"Error downloading HGNC: {e}")
        logger.warning("Falling back to basic synonym mapping")
        return build_basic_synonym_dict()

    # Build efficient lookup dictionary
    hgnc_dict = {}

    for gene in raw_data.get("response", {}).get("docs", []):
        approved_symbol = gene.get("symbol", "")

        if not approved_symbol:
            continue

        # Collect all names for this gene
        aliases = set()

        # Add approved symbol
        aliases.add(approved_symbol)

        # Add alias symbols
        if "alias_symbol" in gene:
            alias_list = gene["alias_symbol"]
            if isinstance(alias_list, list):
                aliases.update(alias_list)
            else:
                aliases.add(alias_list)

        # Add previous symbols
        if "prev_symbol" in gene:
            prev_list = gene["prev_symbol"]
            if isinstance(prev_list, list):
                aliases.update(prev_list)
            else:
                aliases.add(prev_list)

        # Remove dashes and add normalized versions
        aliases_normalized = set()
        for alias in aliases:
            aliases_normalized.add(alias)
            # Add version without dashes for matching
            if "-" in alias:
                aliases_normalized.add(alias.replace("-", ""))

        # Create entry: each alias points to approved symbol
        for alias in aliases_normalized:
            hgnc_dict[alias.upper()] = {
                "approved_symbol": approved_symbol,
                "all_symbols": list(aliases),
            }

    # Cache the result
    with open(cache_path, "w", encoding="utf-8") as f:
        json.dump(hgnc_dict, f, indent=2)

    logger.debug(f"HGNC database loaded ({len(hgnc_dict)} entries) and cached")
    return hgnc_dict


def build_basic_synonym_dict() -> dict:
    """Fallback: Basic common gene synonyms (if download fails)"""
    basic_dict = {
        "IL8": {"approved_symbol": "IL8", "all_symbols": ["IL8", "IL-8", "CXCL8"]},
        "IL-8": {"approved_symbol": "IL8", "all_symbols": ["IL8", "IL-8", "CXCL8"]},
        "CXCL8": {"approved_symbol": "IL8", "all_symbols": ["IL8", "IL-8", "CXCL8"]},
        "TNF": {"approved_symbol": "TNF", "all_symbols": ["TNF", "TNF-ALPHA"]},
        "TNF-ALPHA": {
            "approved_symbol": "TNF",
            "all_symbols": ["TNF", "TNF-ALPHA"],
        },
    }
    return basic_dict


def normalize_gene_name(gene_name: str, hgnc_dict: dict) -> str | None:
    """
    Normalize gene name using HGNC database.
    Handles dash variants and synonyms.

    Parameters
    ----------
    gene_name : str
        Gene name to normalize
    hgnc_dict : dict
        HGNC database dictionary

    Returns
    -------
    str or None
        Approved gene symbol or None if not found
    """

    gene_upper = gene_name.upper()

    # Direct match
    if gene_upper in hgnc_dict:
        return hgnc_dict[gene_upper]["approved_symbol"]

    # Try without dashes
    gene_no_dash = gene_upper.replace("-", "")

    if gene_no_dash in hgnc_dict:
        return hgnc_dict[gene_no_dash]["approved_symbol"]

    # Try with dashes in different positions (for genes like IL17 -> IL-17)
    for i in range(1, len(gene_upper)):
        with_dash = gene_upper[:i] + "-" + gene_upper[i:]
        if with_dash in hgnc_dict:
            return hgnc_dict[with_dash]["approved_symbol"]

    return None


# ============================================================================
# GENE/VARIANT CLASSIFICATION
# ============================================================================


def classify_and_filter_genes(gene_name: str, hgnc_dict: dict) -> tuple:
    """
    Classify gene/variant and filter out invalid/truncated entries.

    Parameters
    ----------
    gene_name : str
        Gene name or variant identifier
    hgnc_dict : dict
        HGNC database dictionary

    Returns
    -------
    tuple
        (is_valid: bool, normalized_name: str, variant_type: str, message: str)
    """

    gene_name = str(gene_name).strip()

    if not gene_name:
        return False, None, None, "Empty gene name"

    # Check for truncated/incomplete genes
    if gene_name.endswith("-") or gene_name.startswith("-"):
        return False, None, "TRUNCATED", f"Truncated gene name: {gene_name}"

    if gene_name in ["HLA", "IL", "TNF", "MHC"]:  # Common truncations
        return False, None, "TRUNCATED", f"Incomplete gene: {gene_name}"

    # Classify by pattern
    rs_pattern = r"^rs\d+$"
    protein_pattern = r"^p\.[A-Za-z0-9\*]+$"
    cdna_pattern = r"^c\.[A-Za-z0-9\+\-_]+$"
    refseq_pattern = r"^NM_\d+.*"

    if re.match(rs_pattern, gene_name):
        logger.debug(f"SNP identifier detected: {gene_name} (not processable)")
        return False, None, "SNP_ID", f"SNP identifier (not processable): {gene_name}"

    if re.match(protein_pattern, gene_name) or re.match(cdna_pattern, gene_name):
        return (
            False,
            None,
            "PROTEIN_VAR",
            f"Protein/cDNA variant (not processable): {gene_name}",
        )

    if re.match(refseq_pattern, gene_name):
        return (
            False,
            None,
            "REFSEQ_ID",
            f"RefSeq identifier (not processable): {gene_name}",
        )

    # Try to normalize the gene name
    normalized = normalize_gene_name(gene_name, hgnc_dict)

    if normalized:
        logger.debug(f"Gene normalized: {gene_name} → {normalized}")
        return True, normalized, "GENE", "Valid gene"
    else:
        # Check if it looks like an attempt at a gene name
        if len(gene_name) > 2 and gene_name.isalpha():
            logger.debug(f"Gene not in HGNC (unmapped): {gene_name}")
            return (
                True,
                gene_name.upper(),
                "GENE_UNMAPPED",
                f"Gene not found in HGNC: {gene_name}",
            )
        else:
            return (
                False,
                None,
                "INVALID",
                f"Invalid gene identifier: {gene_name}",
            )


# ============================================================================
# DEDUPLICATION FUNCTION
# ============================================================================


def deduplicate_search_results(
    genes_list: list,
    variants_list: list,
    hgnc_dict: dict = None,
    variants_with_genes: dict = None,
) -> dict:
    """
    Deduplicate and verify genes/variants from search results.

    Internal function used by verify_genes_with_hgnc wrapper.

    Parameters
    ----------
    genes_list : list
        List of gene names extracted from literature
    variants_list : list
        List of variant identifiers
    hgnc_dict : dict, optional
        HGNC database dictionary. If None, will load it.
    variants_with_genes : dict, optional
        Mapping of variant_id -> gene_name (for variants with associated genes)

    Returns
    -------
    dict
        Complete deduplication results
    """

    if hgnc_dict is None:
        hgnc_dict = load_hgnc_database(use_cache=True)

    if variants_with_genes is None:
        variants_with_genes = {}

    logger.debug(
        f"Deduplicating search results: {len(genes_list)} genes, {len(variants_list)} variants"
    )

    # =========================================================================
    # PROCESS GENES
    # =========================================================================

    genes_verified = []
    genes_unmapped = []
    genes_invalid = []

    # Deduplicate gene list first
    unique_genes = list(set(genes_list))

    for gene_name in unique_genes:
        is_valid, canonical, gene_type, message = classify_and_filter_genes(
            gene_name, hgnc_dict
        )

        if not is_valid:
            logger.debug(f"Gene filtered: {gene_name} - {message}")
            genes_invalid.append(
                {"original": gene_name, "type": gene_type, "reason": message}
            )
        elif gene_type == "GENE_UNMAPPED":
            logger.debug(f"Gene unmapped (not in HGNC): {gene_name}")
            genes_unmapped.append({"original": gene_name, "canonical": canonical})
        else:  # GENE_UNMAPPED or GENE
            genes_verified.append(
                {"original": gene_name, "canonical": canonical, "type": gene_type}
            )

    # Deduplicate by canonical name
    verified_canonicals = list({g["canonical"] for g in genes_verified})

    logger.debug(
        f"Genes verified: {len(verified_canonicals)} canonical (from {len(unique_genes)} unique)"
    )

    # =========================================================================
    # PROCESS VARIANTS
    # =========================================================================

    variants_with_genes_verified = []
    variants_invalid = []

    unique_variants = list(set(variants_list))

    for variant in unique_variants:
        # Check if variant has associated gene
        if variant in variants_with_genes:
            associated_gene = variants_with_genes[variant]

            # Verify the associated gene
            is_valid, canonical, gene_type, message = classify_and_filter_genes(
                associated_gene, hgnc_dict
            )

            if is_valid:
                logger.debug(f"Variant with gene verified: {variant} → {canonical}")
                variants_with_genes_verified.append(
                    {
                        "variant": variant,
                        "associated_gene": associated_gene,
                        "canonical_gene": canonical,
                        "gene_type": gene_type,
                    }
                )
            else:
                logger.debug(
                    f"Variant with invalid gene: {variant} → {associated_gene} ({message})"
                )
                variants_invalid.append(
                    {
                        "variant": variant,
                        "associated_gene": associated_gene,
                        "reason": message,
                    }
                )
        else:
            # Variant without associated gene
            logger.debug(
                f"Variant without associated gene (not processable): {variant}"
            )
            variants_invalid.append(
                {
                    "variant": variant,
                    "associated_gene": None,
                    "reason": "No associated gene",
                }
            )

    # =========================================================================
    # COMPILE DEDUP REPORT
    # =========================================================================

    dedup_report = {
        "input_genes": len(genes_list),
        "unique_genes": len(unique_genes),
        "genes_verified_canonical": len(verified_canonicals),
        "genes_unmapped": len(genes_unmapped),
        "genes_invalid": len(genes_invalid),
        "input_variants": len(variants_list),
        "unique_variants": len(unique_variants),
        "variants_with_genes_verified": len(variants_with_genes_verified),
        "variants_invalid": len(variants_invalid),
    }

    logger.debug(f"Deduplication complete: {dedup_report}")

    return {
        "genes_verified": genes_verified,
        "genes_canonical": verified_canonicals,
        "genes_unmapped": genes_unmapped,
        "genes_invalid": genes_invalid,
        "variants_with_genes": variants_with_genes_verified,
        "variants_invalid": variants_invalid,
        "dedup_report": dedup_report,
        "hgnc_dict": hgnc_dict,
    }


# ============================================================================
# MAIN WRAPPER FUNCTION - SINGLE ENTRY POINT FOR GENLIT
# ============================================================================


def verify_genes_with_hgnc(
    genes_list: list,
    variants_list: list = None,
    use_cache: bool = True,
    hgnc_cache_file: str = "hgnc_cache.json",
) -> dict:
    """
    Main wrapper function for HGNC gene verification and deduplication.

    This is the SINGLE ENTRY POINT for all HGNC functionality in Genlit.
    Handles gene normalization, classification, deduplication, and variant
    verification against HGNC database.

    Parameters
    ----------
    genes_list : list
        List of gene names/identifiers to verify and deduplicate
    variants_list : list, optional
        List of variant identifiers (default: empty list)
    use_cache : bool
        Whether to use cached HGNC database (default: True)
    hgnc_cache_file : str
        Path to HGNC cache file (default: "hgnc_cache.json")

    Returns
    -------
    dict
        Comprehensive results dictionary containing:

        - 'genes_verified': list of verified genes with canonical names
        - 'genes_canonical': list of canonical gene symbols (deduplicated)
        - 'genes_unmapped': list of genes not found in HGNC
        - 'genes_invalid': list of invalid/filtered gene identifiers
        - 'variants_with_genes': list of variants with verified gene associations
        - 'variants_invalid': list of unverified variants
        - 'dedup_report': detailed deduplication statistics
        - 'hgnc_dict': loaded HGNC dictionary for downstream use

    Examples
    --------
    >>> genes = ["IL-8", "IL8", "CXCL8", "TNF"]
    >>> results = verify_genes_with_hgnc(genes)
    >>> print(results['genes_canonical'])
    ['IL8', 'TNF']

    >>> genes = ["IL8", "TNF"]
    >>> variants = ["rs123456"]
    >>> results = verify_genes_with_hgnc(genes, variants)
    >>> print(f"Verified: {len(results['genes_canonical'])}")
    Verified: 2

    Notes
    -----
    - HGNC database is cached locally (~40MB) on first run
    - Subsequent calls use cached version (instant)
    - All logging output is at DEBUG level
    - Fallback to basic synonym dictionary if HGNC unavailable
    """

    if variants_list is None:
        variants_list = []

    logger.debug(
        f"Gene verification started: {len(genes_list)} genes, "
        f"{len(variants_list)} variants"
    )

    # Load HGNC database
    hgnc_dict = load_hgnc_database(use_cache=use_cache, cache_file=hgnc_cache_file)

    # Deduplicate and verify
    results = deduplicate_search_results(
        genes_list=genes_list,
        variants_list=variants_list,
        hgnc_dict=hgnc_dict,
        variants_with_genes=None,
    )

    logger.debug(
        f"Gene verification complete: {len(results['genes_canonical'])} verified, "
        f"{len(results['genes_unmapped'])} unmapped, "
        f"{len(results['genes_invalid'])} invalid"
    )

    return results
