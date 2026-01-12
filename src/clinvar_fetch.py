#!/usr/bin/env python3
# Author: Diego Fuentes
# Contact email: diegofupa@gmail.com
# Barcelona
# Date:2025-11-27

from Bio import Entrez
from xml.etree import ElementTree as ET
import pandas as pd


# -------------------------
# Entrez configuration
# -------------------------


def configure_entrez(email, api_key=None, tool="gene_textminer"):
    """
    Configure Entrez globally. Must be called once.
    Select the desired tool.
    Api key if necessary (shouldn't be for basic fetching)
    """
    Entrez.email = email
    Entrez.tool = tool

    if api_key:
        Entrez.api_key = api_key


# --------------------------------------
# ClinVar Search
# --------------------------------------


def clinvar_search(query, retmax=200):
    """
    Search ClinVar using Entrez.esearch. Retrival by default set to 200.
    Query can be: gene symbol, rsID, HGVS string, condition, etc.
    Based on ClinVar query synthax
    Returns a list of ClinVar IDs (accessions).
    """
    handle = Entrez.esearch(db="clinvar", term=query, retmax=retmax, retmode="xml")
    record = Entrez.read(handle, validate=False)
    handle.close()
    return record.get("IdList", [])


# --------------------------------------
# ClinVar Summary Fetch
# --------------------------------------


def clinvar_summaries(ids):
    """
    Fetch ClinVar summary records (JSON-like dictionary).
    Summaries provide clinical significance, gene names, variants, etc.
    Returns a list of summary dicts.
    """
    if isinstance(ids, list):
        ids = ",".join(ids)

    handle = Entrez.esummary(db="clinvar", id=ids, retmode="xml")
    summary = Entrez.read(handle, validate=False)
    handle.close()
    # No need to further parse, already parsed into Python dicts
    return summary["DocumentSummarySet"]["DocumentSummary"]


# --------------------------------------
# ClinVar Full-Record Fetch
# --------------------------------------


def clinvar_fetch_full(ids):
    """
    Fetch full XML ClinVar records.
    Useful for detailed parsing: submissions, alleles, review status, etc.
    Needs to be parsed downstream, returns XML ElementTree root.
    """
    if isinstance(ids, list):
        ids = ",".join(ids)

    handle = Entrez.efetch(db="clinvar", id=ids, retmode="xml")
    xml_data = handle.read()
    handle.close()

    root = ET.fromstring(xml_data)
    return root


# --------------------------------------
# Helper: Extract useful info from summary dicts
# --------------------------------------


def parse_clinvar_summary(summary):
    """
    Convert a raw ClinVar DictElement into a clean Python dict
    with all relevant fields extracted.
    """

    def get(obj, key, default=None):
        try:
            if hasattr(obj, "get"):
                return obj.get(key, default)
        except Exception:
            pass
        return default

    def norm_list(x):
        if x is None:
            return []
        return list(x) if isinstance(x, (list, tuple)) else [x]

    clean = {}

    # --- Basic Identifiers ---
    clean["uid"] = summary.attributes.get("uid")
    clean["accession"] = get(summary, "accession")
    clean["accession_version"] = get(summary, "accession_version")
    clean["title"] = get(summary, "title")
    clean["variant_type"] = get(summary, "obj_type")

    # --- Variation Set (primary variant block) ---
    vset = get(summary, "variation_set", [])
    v = vset[0] if vset else {}

    clean["measure_id"] = get(v, "measure_id")
    clean["variation_name"] = get(v, "variation_name")
    clean["cdna_change"] = get(v, "cdna_change")
    clean["aliases"] = norm_list(get(v, "aliases"))
    clean["canonical_spdi"] = get(v, "canonical_spdi")

    # --- Molecular Consequences ---
    clean["molecular_consequences"] = norm_list(
        get(summary, "molecular_consequence_list", [])
    )

    # --- Protein Change ---
    clean["protein_change"] = get(summary, "protein_change")

    # --- Gene Information ---
    genes = get(summary, "genes", [])
    clean["genes"] = [
        {
            "symbol": get(g, "symbol"),
            "gene_id": get(g, "GeneID"),
            "strand": get(g, "strand"),
        }
        for g in genes
    ]

    # --- Genomic Coordinates (all assemblies) ---
    coords = []
    for loc in norm_list(get(v, "variation_loc", [])):
        coords.append(
            {
                "assembly": get(loc, "assembly_name"),
                "status": get(loc, "status"),
                "chr": get(loc, "chr"),
                "band": get(loc, "band"),
                "start": get(loc, "start"),
                "stop": get(loc, "stop"),
                "assembly_accession": get(loc, "assembly_acc_ver"),
            }
        )
    clean["coordinates"] = coords

    # --- Supporting Submissions ---
    subs = get(summary, "supporting_submissions", {})
    clean["scv"] = norm_list(get(subs, "scv"))
    clean["rcv"] = norm_list(get(subs, "rcv"))

    # --- Germline Classification ---
    germ = get(summary, "germline_classification", {})
    clean["clinical_significance"] = {
        "description": get(germ, "description"),
        "last_evaluated": get(germ, "last_evaluated"),
        "review_status": get(germ, "review_status"),
        "traits": [
            {
                "trait_name": get(t, "trait_name"),
                "trait_xrefs": [
                    {"db": get(x, "db_source"), "id": get(x, "db_id")}
                    for x in norm_list(get(t, "trait_xrefs", []))
                ],
            }
            for t in norm_list(get(germ, "trait_set", []))
        ],
    }

    # --- Oncogenicity & Impact (rarely used but included) ---
    clean["oncogenicity"] = get(summary, "oncogenicity_classification")
    clean["clinical_impact"] = get(summary, "clinical_impact_classification")

    # --- Sorting & Misc ---
    clean["gene_sort"] = get(summary, "gene_sort")
    clean["chr_sort"] = get(summary, "chr_sort")
    clean["location_sort"] = get(summary, "location_sort")

    return clean


def flatten_clinvar_summary(clean):
    """
    Flatten one parsed ClinVar summary (a dict) into
    a single-row pandas DataFrame.
    """

    def join_list(x):
        if not x:
            return None
        if isinstance(x, list):
            return ", ".join(str(i) for i in x)
        return str(x)

    # --- Coordinate extraction ---
    coords38 = next(
        (c for c in clean.get("coordinates", []) if c.get("assembly") == "GRCh38"), {}
    )
    coords37 = next(
        (c for c in clean.get("coordinates", []) if c.get("assembly") == "GRCh37"), {}
    )

    # --- Trait extraction ---
    traits = clean.get("clinical_significance", {}).get("traits", [])
    trait_names = join_list([t.get("trait_name") for t in traits])
    trait_ids = join_list(
        [
            f"{xref.get('db')}:{xref.get('id')}"
            for t in traits
            for xref in t.get("trait_xrefs", [])
        ]
    )

    # --- Gene extraction ---
    genes = clean.get("genes", [])
    gene_symbols = join_list([g.get("symbol") for g in genes])
    gene_ids = join_list([g.get("gene_id") for g in genes])

    # --- Molecular consequences ---
    consequences = join_list(clean.get("molecular_consequences"))

    # --- Build flat dict ---
    flat = {
        "uid": clean.get("uid"),
        "accession": clean.get("accession"),
        "accession_version": clean.get("accession_version"),
        "title": clean.get("title"),
        "variant_type": clean.get("variant_type"),
        "variation_name": clean.get("variation_name"),
        "measure_id": clean.get("measure_id"),
        "cdna_change": clean.get("cdna_change"),
        "protein_change": clean.get("protein_change"),
        "aliases": join_list(clean.get("aliases")),
        "canonical_spdi": clean.get("canonical_spdi"),
        # Genes
        "genes": gene_symbols,
        "gene_ids": gene_ids,
        # Coordinates GRCh38
        "chr38": coords38.get("chr"),
        "start38": coords38.get("start"),
        "stop38": coords38.get("stop"),
        # Coordinates GRCh37
        "chr37": coords37.get("chr"),
        "start37": coords37.get("start"),
        "stop37": coords37.get("stop"),
        # Clinical significance
        "clinical_significance": clean.get("clinical_significance", {}).get(
            "description"
        ),
        "review_status": clean.get("clinical_significance", {}).get("review_status"),
        "last_evaluated": clean.get("clinical_significance", {}).get("last_evaluated"),
        "trait_names": trait_names,
        "trait_ids": trait_ids,
        # Molecular consequences
        "molecular_consequences": consequences,
        # Submissions
        "scv": join_list(clean.get("scv")),
        "rcv": join_list(clean.get("rcv")),
    }

    return pd.DataFrame([flat])


# --------------------------------------
# Composite: Search → Summary → Parsed
# --------------------------------------


def search_and_fetch_clinvar(query, retmax=50, parse=True):
    """
    Convenience wrapper:
      - search ClinVar
      - fetch summary metadata
      - optionally parse into structured dicts
    """
    ids = clinvar_search(query, retmax=retmax)
    if not ids:
        return []
    summaries = clinvar_summaries(ids)

    if parse:
        # return [parse_clinvar_summary(s) for s in summaries]
        return [flatten_clinvar_summary(parse_clinvar_summary(s)) for s in summaries]

    return summaries


# --------------------------------------
# Example usage
# --------------------------------------

if __name__ == "__main__":
    configure_entrez(email="diegofupa@gmail.com")

    # example: variants associated to Hidradenitis Suppurativa
    results = search_and_fetch_clinvar("Hidradenitis Suppurativa", retmax=20)

    df = pd.concat(results, ignore_index=True)

    print(df.head())
