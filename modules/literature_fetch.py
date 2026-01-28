#!/usr/bin/env python3
"""
Author: Diego Fuentes
Contact email: diegofupa@gmail.com
Barcelona
Date:2025-11-27

Module for fetching literature data from NCBI Entrez (PubMed and PMC).
Can be used to search PubMed, fetch article metadata and abstracts, convert between PMIDs and PMCIDs
and fetch full-text XML from PMC. Can be executed as is to demonstrate functionality.
"""


from Bio import Entrez
from xml.etree import ElementTree as ET
import time

# -------------------------
# Entrez configuration
# -------------------------


def configure_entrez(
    email: str, api_key: str = None, tool: str = "gene_textminer"
) -> None:
    """
    Configure Entrez globally. Must be called once.
    Select the desired tool.
    Api key if necessary (shouldn't be necessary for basic fetching)
    """
    Entrez.email = email
    Entrez.tool = tool

    if api_key:
        Entrez.api_key = api_key


def entrez_with_retry(func: callable, max_retries: int = 3, delay: float = 1.0) -> None:
    for attempt in range(max_retries):
        try:
            return func()
        except RuntimeError as e:
            if attempt == max_retries - 1:
                raise e
            time.sleep(delay * (attempt + 1))


# -------------------------
# PubMed Search
# -------------------------


def pubmed_search(query: str, retmax: int = 500) -> list[str]:
    """
    Run a PubMed search using Entrez.esearch.
    Retrieve a maximum of 500 papers by default.
    Query is the most important parameter: uses PubMed synthax
    Returns a list of PMIDs.
    """
    handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax, retmode="xml")
    record = entrez_with_retry(lambda: Entrez.read(handle))
    handle.close()
    return record.get("IdList", [])


# -------------------------
# PubMed Fetch
# -------------------------


def pubmed_fetch_details(pmids: list[str]) -> list[dict]:
    """
    Fetch article metadata and abstracts from PubMed.
    Returns a list of dicts.
    """
    if isinstance(pmids, list):
        pmids = ",".join(pmids)

    handle = Entrez.efetch(db="pubmed", id=pmids, rettype="xml", retmode="xml")
    data = entrez_with_retry(lambda: handle.read())
    handle.close()

    root = ET.fromstring(data)
    out = []

    for art in root.findall(".//PubmedArticle"):
        try:
            pmid = art.find("./MedlineCitation/PMID").text
            title = art.findtext(".//ArticleTitle", default="")

            abstract_parts = [
                elem.text
                for elem in art.findall(".//Abstract/AbstractText")
                if elem is not None and elem.text
            ]
            abstract = " ".join(abstract_parts).strip()

            out.append({"pmid": pmid, "title": title, "abstract": abstract or None})
        except Exception:
            continue

    return out


# -------------------------
# PMCID to PMID Conversion
# -------------------------


def convert_pmid_to_pmcid(pmid: str) -> str | None:
    """
    Convert a PMID to PMCID using Entrez.elink.
    Just in case the article has a full text available
    Returns PMCID or None.
    """
    handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid, linkname="pubmed_pmc")

    try:
        record = entrez_with_retry(lambda: Entrez.read(handle))
        handle.close()
        link = record[0]["LinkSetDb"][0]["Link"][0]["Id"]
        return f"PMC{link}"
    except Exception:
        return None


# -------------------------
# PMC Full Text Fetch (XML)
# -------------------------


def pmc_fetch_fulltext(pmcid: str) -> ET.Element | None:
    """
    Fetch the full-text XML for a given PMCID.
    Needs to be parsed downstream, returns an XML ElementTree root.
    """
    pmcid_clean = pmcid.replace("PMC", "")

    handle = Entrez.efetch(
        db="pmc", id=pmcid_clean, rettype="full", retmode="xml", timeout=60
    )
    xml_content = entrez_with_retry(lambda: handle.read())
    handle.close()

    try:
        root = ET.fromstring(xml_content)
        return root
    except Exception:
        return None


# -------------------------
# Convenience composite
# -------------------------


def search_and_fetch_pubmed(
    query: str, retmax: int = 200, fetch_pmc: bool = False
) -> list[dict]:
    """
    Searches PubMed and fetches details.
    Fetches titles + abstracts (and optionally full PMC text).
    Full PMC is available only if the article has a PMCID.
    Returns a list of dicts.
    """
    pmids = pubmed_search(query, retmax=retmax)
    papers = pubmed_fetch_details(pmids)

    if fetch_pmc:
        for paper in papers:
            pmcid = convert_pmid_to_pmcid(paper["pmid"])
            if pmcid:
                paper["pmcid"] = pmcid
                paper["pmc_xml"] = pmc_fetch_fulltext(pmcid)
            else:
                paper["pmcid"] = None
                paper["pmc_xml"] = None

    return papers


# -------------------------
# Example usage
# -------------------------

if __name__ == "__main__":
    configure_entrez(email="examplemail@gmail.com")

    query = "Hidradenitis Suppurativa AND (genetics OR variant)"
    results = search_and_fetch_pubmed(query, retmax=10, fetch_pmc=True)

    for r in results:
        print(f"{r['pmid']} — {r['title']}")
        if r["abstract"]:
            print(r["abstract"][:200], "...\n")
        # Optional, but necessary to decode the PMC fulltext
        if r["pmc_xml"]:
            a = r["pmc_xml"]
            print(ET.tostring(a, encoding="unicode")[:200], "...\n")
