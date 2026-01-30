#!/usr/bin/env python3

"""
Author: Diego Fuentes
Contact email: diegofupa@gmail.com
Barcelona
Date: 2026-01-30

Module for fetching literature data from NCBI Entrez (PubMed and PMC).
Can be used to search PubMed, fetch article metadata and abstracts, convert between PMIDs and PMCIDs
and fetch full-text XML from PMC.
It uses retry logic with exponential backoff for all Entrez API calls.
Handles HTTP errors, timeouts, connection errors, and malformed responses.

Can be executed as is to demonstrate functionality.
"""

from Bio import Entrez
from Bio.Entrez.Parser import NotXMLError
from xml.etree import ElementTree as ET
import time
import random
import logging
from typing import Optional, List, Dict, Callable, Any
from urllib.error import HTTPError, URLError
import socket

# Configure logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

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
    logger.info(f"Entrez configured with email: {email}, tool: {tool}")


# -------------------------
# Robust Retry Logic with Backoff
# -------------------------


def entrez_with_retry(
    func: Callable,
    max_retries: int = 5,
    base_delay: float = 1.0,
    function_name: str = "Entrez API",
    verbose: bool = False,
) -> Any:
    """
    Execute a function with exponential backoff retry logic.
    Handles multiple error types: HTTP errors, connection errors, timeouts, and parsing errors.

    Args:
        func: Callable function to execute
        max_retries: Maximum retry attempts (default 5)
        base_delay: Base delay in seconds for exponential backoff (default 1.0)
        function_name: Name of function for logging purposes
        verbose: Print retry details if True

    Returns:
        Result of the function call

    Raises:
        Exception: If all retries fail
    """

    retriable_errors = (
        HTTPError,  # HTTP errors (500, 503, 429, etc.)
        URLError,  # Connection/network errors
        socket.timeout,  # Socket timeout
        TimeoutError,  # Timeout errors
        ConnectionError,  # Connection errors
        RuntimeError,  # Generic runtime errors
        BrokenPipeError,  # Broken pipe
        OSError,  # OS-level errors
        NotXMLError,  # NOTXmL error
    )

    for attempt in range(max_retries):
        try:
            return func()

        except retriable_errors as e:
            is_final_attempt = attempt == max_retries - 1
            error_type = type(e).__name__

            # Extract HTTP error code if available
            http_code = None
            if isinstance(e, HTTPError):
                http_code = e.code

            if is_final_attempt:
                error_msg = (
                    f"[ERROR] Max retries ({max_retries}) exceeded for {function_name}"
                )
                if http_code:
                    error_msg += f" (HTTP {http_code})"
                logger.error(error_msg)
                raise e

            # Calculate backoff with jitter
            delay = (base_delay * (2**attempt)) + random.uniform(0, 1)

            if verbose or logger.level <= logging.INFO:
                msg = f"[RETRY {attempt + 1}/{max_retries}] {function_name} failed: {error_type}"
                if http_code:
                    msg += f" (HTTP {http_code})"
                msg += f" — Retrying in {delay:.1f}s..."
                logger.debug(msg)

            time.sleep(delay)

        except Exception as e:
            # Non-retriable errors (parsing errors, etc.)
            error_type = type(e).__name__
            logger.debug(
                f"[ERROR] Non-retriable error in {function_name}: {error_type}: {str(e)[:-1]}"
            )
            raise e


# -------------------------
# PubMed Search
# -------------------------


def pubmed_search(
    query: str, retmax: int = 500, max_retries: int = 5, verbose: bool = False
) -> List[str]:
    """
    Run a PubMed search using Entrez.esearch.
    Retrieve a maximum of 500 papers by default.
    Query uses PubMed syntax.

    Args:
        query: PubMed search query string
        retmax: Maximum number of results to return (default 500)
        max_retries: Maximum retry attempts (default 5)
        verbose: Print verbose output if True

    Returns:
        List of PMIDs
    """

    def _search():
        handle = Entrez.esearch(
            db="pubmed", term=query, retmax=retmax, retmode="xml", timeout=30
        )
        return handle

    handle = entrez_with_retry(
        _search,
        max_retries=max_retries,
        base_delay=2.0,
        function_name=f"pubmed_search({query[:50]})",
        verbose=verbose,
    )

    record = entrez_with_retry(
        lambda: Entrez.read(handle),
        max_retries=max_retries,
        base_delay=4.0,
        function_name="PubMed search result parsing",
        verbose=verbose,
    )
    handle.close()

    pmid_list = record.get("IdList", [])
    logger.info(f"PubMed search found {len(pmid_list)} articles")
    return pmid_list


# -------------------------
# PubMed Fetch
# -------------------------


def pubmed_fetch_details(
    pmids: List[str], max_retries: int = 5, verbose: bool = False
) -> List[Dict]:
    """
    Fetch article metadata and abstracts from PubMed.
    Handles batching for large PMID lists.

    Args:
        pmids: List of PMIDs to fetch
        max_retries: Maximum retry attempts (default 5)
        verbose: Print verbose output if True

    Returns:
        List of dictionaries with pmid, title, abstract
    """
    if isinstance(pmids, list):
        pmids_str = ",".join(pmids)
    else:
        pmids_str = pmids

    def _fetch():
        handle = Entrez.efetch(
            db="pubmed", id=pmids_str, rettype="xml", retmode="xml", timeout=60
        )
        return handle.read()

    data = entrez_with_retry(
        _fetch,
        max_retries=max_retries,
        base_delay=2.0,
        function_name=f"pubmed_fetch_details({len(pmids)} articles)",
        verbose=verbose,
    )

    out = []
    try:
        root = ET.fromstring(data)
    except ET.ParseError as e:
        logger.error(f"Failed to parse PubMed XML response: {e}")
        return out

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
        except Exception as e:
            # Skip malformed articles
            logger.debug(f"Skipping malformed article: {e}")
            continue

    logger.info(f"Successfully fetched details for {len(out)}/{len(pmids)} articles")
    return out


# -------------------------
# PMCID to PMID Conversion
# -------------------------


def convert_pmid_to_pmcid(
    pmid: str, max_retries: int = 5, verbose: bool = False
) -> Optional[str]:
    """
    Convert a PMID to PMCID using Entrez.elink.
    Returns PMCID if article has full text available in PMC.

    Args:
        pmid: PubMed ID to convert
        max_retries: Maximum retry attempts (default 5)
        verbose: Print verbose output if True

    Returns:
        PMCID string (e.g., "PMC123456") or None if not found
    """

    def _elink():
        handle = Entrez.elink(
            dbfrom="pubmed", db="pmc", id=pmid, linkname="pubmed_pmc", timeout=60
        )
        return handle

    try:
        handle = entrez_with_retry(
            _elink,
            max_retries=max_retries,
            base_delay=2.0,
            function_name=f"convert_pmid_to_pmcid({pmid})",
            verbose=verbose,
        )

        record = entrez_with_retry(
            lambda: Entrez.read(handle),
            max_retries=max_retries,
            base_delay=3.0,
            function_name=f"PMCID lookup parsing",
            verbose=verbose,
        )
        handle.close()

        # Try to extract PMCID from response
        try:
            link = record[0]["LinkSetDb"][0]["Link"][0]["Id"]
            pmcid = f"PMC{link}"
            logger.debug(f"Found PMCID {pmcid} for PMID {pmid}")
            return pmcid
        except (IndexError, KeyError, TypeError):
            # Article doesn't have PMC full text
            logger.debug(f"No PMC entry found for PMID {pmid}")
            return None

    except Exception as e:
        logger.warning(
            f"Failed to convert PMID {pmid} to PMCID: {type(e).__name__}: {str(e)[:-1]}"
        )
        return None


# -------------------------
# PMC Full Text Fetch (XML)
# -------------------------


def pmc_fetch_fulltext(
    pmcid: str, max_retries: int = 5, verbose: bool = False
) -> Optional[ET.Element]:
    """
    Fetch the full-text XML for a given PMCID.
    Returns an XML ElementTree root element.

    Args:
        pmcid: PMC ID to fetch (with or without "PMC" prefix)
        max_retries: Maximum retry attempts (default 5)
        verbose: Print verbose output if True

    Returns:
        XML ElementTree root or None if fetch failed
    """
    pmcid_clean = pmcid.replace("PMC", "")

    def _fetch_xml():
        handle = Entrez.efetch(
            db="pmc", id=pmcid_clean, rettype="full", retmode="xml", timeout=60
        )
        return handle.read()

    try:
        xml_content = entrez_with_retry(
            _fetch_xml,
            max_retries=max_retries,
            base_delay=3.0,
            function_name=f"pmc_fetch_fulltext({pmcid})",
            verbose=verbose,
        )

        try:
            root = ET.fromstring(xml_content)
            logger.debug(f"Successfully fetched PMC full text for available id {pmcid}")
            return root
        except ET.ParseError as e:
            logger.warning(f"Failed to parse PMC XML for id {pmcid}: {e}")
            return None

    except Exception as e:
        logger.warning(
            f"Failed to fetch PMC full text for id {pmcid}: {type(e).__name__}: {str(e)[:-1]}"
        )
        return None


# -------------------------
# Convenience composite
# -------------------------


def search_and_fetch_pubmed(
    query: str,
    retmax: int = 200,
    fetch_pmc: bool = True,
    max_retries: int = 5,
    verbose: bool = False,
) -> List[Dict]:
    """
    Searches PubMed and fetches details with full retry logic.
    Optionally fetches full PMC text if available.

    Args:
        query: PubMed search query
        retmax: Maximum number of articles to retrieve (default 200)
        fetch_pmc: Whether to fetch full PMC text (default True)
        max_retries: Maximum retry attempts for API calls (default 5)
        verbose: Print verbose output if True

    Returns:
        List of dictionaries with article data:
        - pmid: PubMed ID
        - title: Article title
        - abstract: Article abstract (or None)
        - pmcid: PMC ID (if available and fetch_pmc=True)
        - pmc_xml: XML ElementTree root (if available and fetch_pmc=True)
    """
    logger.info(f"Starting PubMed search: '{query}' (max {retmax} articles)")

    # Step 1: Search PubMed
    pmids = pubmed_search(
        query, retmax=retmax, max_retries=max_retries, verbose=verbose
    )

    if not pmids:
        logger.warning("No articles found for this query")
        return []

    # Step 2: Fetch article details
    papers = pubmed_fetch_details(pmids, max_retries=max_retries, verbose=verbose)

    # Step 3: Optionally fetch PMC full text
    if fetch_pmc:
        logger.info(f"Fetching PMC full text for {len(papers)} articles...")
        for i, paper in enumerate(papers, 1):
            try:
                pmcid = convert_pmid_to_pmcid(
                    paper["pmid"], max_retries=max_retries, verbose=verbose
                )

                if pmcid:
                    pmc_xml = pmc_fetch_fulltext(
                        pmcid, max_retries=max_retries, verbose=verbose
                    )
                    paper["pmcid"] = pmcid
                    paper["pmc_xml"] = pmc_xml
                else:
                    paper["pmcid"] = None
                    paper["pmc_xml"] = None

                if verbose and i % 5 == 0:
                    logger.info(f"Processed {i}/{len(papers)} articles for PMC data")

            except Exception as e:
                logger.warning(
                    f"Error processing PMC data for PMID {paper['pmid']}: {e}"
                )
                paper["pmcid"] = None
                paper["pmc_xml"] = None

        pmc_success = sum(1 for p in papers if p.get("pmc_xml") is not None)
        logger.info(
            f"Successfully fetched PMC full text for {pmc_success} for {len(papers)} articles"
        )

    return papers


# -------------------------
# Example usage
# -------------------------

if __name__ == "__main__":
    configure_entrez(email="your.email@gmail.com")

    query = "Pancreatic Cancer AND (genetics OR variant)"
    logger.info(f"Fetching articles for: {query}")

    results = search_and_fetch_pubmed(query, retmax=50, fetch_pmc=True, verbose=True)

    logger.info(f"\nResults ({len(results)} articles):")
    for r in results:
        logger.info(f"\nPMID: {r['pmid']}")
        logger.info(f"Title: {r['title']}")

        if r.get("pmcid"):
            logger.info(f"PMCID: {r['pmcid']}")

        if r["abstract"]:
            logger.info(f"Abstract: {r['abstract'][:200]}...")

        if r.get("pmc_xml"):
            logger.info(f"Full text: Available (XML)")
