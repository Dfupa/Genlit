#!/usr/bin/env python3
# Author: Diego Fuentes
# Contact email: diegofupa@gmail.com
# Barcelona
# Date:2025-12-03

import warnings
import re
import flair
import time
from flair.data import Sentence
from flair.models import EntityMentionLinker
from flair.tokenization import SciSpacyTokenizer
from flair.splitter import SciSpacySentenceSplitter
from flair.nn import Classifier
from literature_fetch import search_and_fetch
from clinvar_fetch import search_and_fetch_clinvar

start_time = time.time()

# Handle warnings
warnings.filterwarnings("ignore", category=FutureWarning)

flair.device = "cpu"  # Use CPU for compatibility

# ---------------------------
# Set up ther Flair NER model
# ---------------------------
# The model in question is based on hunflair2, accessible from Flair
# It is a biomedical NER tagger trained on the BioNLP13CG corpus
# and other datasets for better performance in biomedical text
# Make sure to install through pip flair, scispacy and the model before running this code:
# pip install -U spacy==3.7.4
# pip install -U scispacy
# pip install -U en_core_sci_sm
# pip install flair

sentence = Sentence(abstract, use_tokenizer=SciSpacyTokenizer())

# initialize the sentence splitter
splitter = SciSpacySentenceSplitter()

# load biomedical NER tagger hunflair2
tagger = Classifier.load("hunflair2")

# You can also use other domain-specific linkers available in Flair:
# Classifier.load("hunflair2-gene")
# Classifier.load("hunflair2-disease")
# Classifier.load("hunflair2-species")
# Classifier.load("hunflair2-chemical")


# load gene linker and perform normalization
gene_linker = EntityMentionLinker.load("gene-linker")


# load disease linker and perform normalization
disease_linker = EntityMentionLinker.load("disease-linker")


# load species linker and perform normalization
species_linker = EntityMentionLinker.load("species-linker")


def extract_entities(text):
    """
    Return gene and protein mentions from text.
    Uses the NLP pipeline previously defined
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
            if entity.value == "Gene":
                genes.append(entity.data_point.text)
    return genes


# -------------------------
# Helper regexpresion func.
# -------------------------
def extract_variants(text):
    """
    This helper function extracts rs and p.
    Returns a dictionary with lists of HGVS and rsIDs from the full text
    """
    # HGVS protein-level
    p_hgvs = re.compile(r"p\.[A-Za-z]+[0-9]+[A-Za-z]+")
    # dbSNP IDs
    rs_regex = re.compile(r"rs[0-9]+")

    prots = p_hgvs.findall(text)
    rsids = rs_regex.findall(text)
    if len(prots) > 0:
        return {"hgvs": prots, "rsids": rsids}
    elif len(rsids) > 0:
        return {"hgvs": prots, "rsids": rsids}
    else:
        return {}


# -------------------------
# Interlink NLP + Clinvar
# -------------------------
def run_search(
    disease_term,
    search_terms=" AND (genetics OR variant)",
    clinical_relevance=["Pathogenic", "Likely pathogenic"],
    retnumber=10,
):
    """
    This function performs the search based on the disease term.
    Uses a NCBI query by default, as well as filtering by clinical relevance
    Returns a dictionary of xlinked terms
    """
    print(f"Searching for: {disease_term}\n")
    results = search_and_fetch(
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
        genes = extract_entities(text)
        abstract_genes.extend(genes)
        abstract_VAR = extract_variants(text)
        abstract_variants.extend(abstract_VAR)

        if r["pmc_xml"]:
            full_text = ET.tostring(r["pmc_xml"], encoding="unicode")
            # Extract entities from full text
            full_text_g = extract_entities(full_text)
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
    print("\n\n--- %s seconds ---" % (time.time() - start_time))
