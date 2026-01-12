#!/usr/bin/env python3
# Author: Diego Fuentes
# Contact email: diegofupa@gmail.com
# Barcelona
# Date:2025-12-03

import spacy
from scispacy.umls_linking import UmlsEntityLinker
from literature_fetch import search_and_fetch
from clinvar_fetch import search_and_fetch_clinvar
import warnings
import re

# Handle warnings
warnings.filterwarnings("ignore", category=FutureWarning)


# -------------------------
# Set up ther NER model
# -------------------------
# The model en_ner_bionlp13cg_md is specialized in biomedical NER for cancer geonomics.
# It can recognize genes and proteins among other entities
# Was downloaded from: https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.3/en_ner_bionlp13cg_md-0.5.3.tar.gz
# Make sure to install through pip scispacy and the model before running this code

nlp = spacy.load("en_ner_bionlp13cg_md")
linker = UmlsEntityLinker(resolve_abbreviations=True, name="umls")
nlp.add_pipe("scispacy_linker", config={"resolve_abbreviations": True})


def extract_entities(text):
    """
    Return gene and protein mentions from text.
    Uses the NLP pipeline previously defined
    Returns a list of entities
    """
    doc = nlp(text)
    genes = []
    for ent in doc.ents:
        if ent.label_ in ("GENE_OR_GENE_PRODUCT", "GENE", "PROTEIN"):
            genes.append(ent.text)
    return genes


# -------------------------
# Helper regexpresion func.
# -------------------------
def extract_variants(text):
    """
    This helper function extracts rs and p.
    """
    # HGVS protein-level
    p_hgvs = re.compile(r"p\.[A-Za-z]+[0-9]+[A-Za-z]+")
    # dbSNP IDs
    rs_regex = re.compile(r"rs[0-9]+")

    prots = p_hgvs.findall(text)
    rsids = rs_regex.findall(text)
    return {"hgvs": prots, "rsids": rsids}


# -------------------------
# Interlink NLP + Clinvar
# -------------------------
def run_search(
    disease_term,
    search_terms=" AND (genetics OR variant)",
    clinical_relevance=["Pathogenic", "Likely pathogenic"],
):
    """
    This function performs the search based on the disease term.
    Uses a NCBI query by default, as well as filtering by clinical relevance
    Returns a dictionary of xlinked terms
    """
    pmids = [r["pmid"] for r in results]
    print(f"Searching for: {disease_term}")
    results = search_and_fetch(
        str(disease_term + search_terms), retmax=10, fetch_pmc=True
    )

    abstract_genes = []
    full_text_genes = []
    abstract_variants = []
    full_text_variants = []

    for r in results:
        text = r["abstract"]
        if not text:
            continue

        # NER for abstracts
        genes = extract_entities(text)
        abstract_genes.extend(genes)
        abstract_VAR = extract_variants(text)
        abstract_variants.extend(abstract_VAR)

        if r["pmc_xml"]:
            full_text = ET.tostring(a, encoding="unicode")
            # NER for full text
            full_text_g = extract_entities(full_text)
            full_text_genes.extend(full_text_g)
            full_VAR = extract_variants(full_text)
            full_text_variants.extend(full_VAR)

        # ClinVar variants
        ClinVar = search_and_fetch_clinvar(disease_term, retmax=20)
        df = pd.concat(ClinVar, ignore_index=True)
        # Filter by clnical relevance:
        fdf = df[df["clinical_significance"].isin(clinical_relevance)]

        clinvar_genes = (
            fdf["genes"]
            .dropna()
            .apply(lambda x: [g.strip() for g in x.split(",")])
            .explode()
            .unique()
        )

        clinvar_variants = fdf["cdna_change"].dropna().unique()
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


# Working example
if __name__ == "__main__":
    query_hs = "Hidradenitis Suppurativa"
    hs_results = run_search(query_hs)
    print("HS genes:", hs_results)
