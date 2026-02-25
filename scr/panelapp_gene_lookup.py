#!/usr/bin/env python3
"""
Author: Diego Fuentes
Contact email: diegofupa@gmail.com
Barcelona
Date: 2026-02-23

Module to query PanelApp API for gene membership, evidence color and genome annotation.
Designed to be run as a script with an input CSV of gene symbols and output an annotated CSV
with one row per gene–panel association, including genome annotation and evidence color.
"""

import argparse
import csv
import sys
import time
from typing import Dict, List, Optional

import requests


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Query PanelApp for gene membership, evidence color and genome annotation."
    )
    parser.add_argument(
        "input_csv", help="Input CSV file with a gene column (default: Gene_Canonical)."
    )
    parser.add_argument(
        "-o",
        "--output-csv",
        required=True,
        help="Output CSV file with PanelApp annotation.",
    )
    parser.add_argument(
        "--gene-column",
        default="Gene_Canonical",
        help="Name of the column containing gene symbols (default: Gene_Canonical).",
    )
    parser.add_argument(
        "--base-url",
        default="https://panelapp.genomicsengland.co.uk/api/v1",
        help="Base URL for PanelApp API (default: https://panelapp.genomicsengland.co.uk/api/v1).",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.2,
        help="Seconds to sleep between API calls (default: 0.2).",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=20.0,
        help="Request timeout in seconds (default: 20).",
    )
    return parser.parse_args(argv)


def read_genes_from_csv(path: str, gene_column: str) -> List[str]:
    """Read gene symbols from a specified column in a CSV file. De-duplicates while preserving order.
    Arguments:
        - path: Path to the input CSV file.
        - gene_column: Name of the column containing gene symbols.
    Returns:     List of unique gene symbols in the order they appear in the file.
     Raises:
        - ValueError if the specified gene column is not found in the CSV header.
    """
    genes: List[str] = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if gene_column not in reader.fieldnames:
            raise ValueError(
                f"Gene column '{gene_column}' not found in CSV header: {reader.fieldnames}"
            )
        for row in reader:
            gene = row.get(gene_column, "").strip()
            if gene:
                genes.append(gene)
    # Keep order but de-duplicate
    seen = set()
    unique_genes = []
    for g in genes:
        if g not in seen:
            seen.add(g)
            unique_genes.append(g)
    return unique_genes


def confidence_to_color(conf: Optional[str]) -> str:
    """
    Map PanelApp confidence_level (0/1/2/3) to color.
    3 -> green, 2 -> amber, 1 -> red, 0 -> grey, other/None -> unknown.
    """
    if conf is None:
        return "unknown"
    try:
        lvl = int(conf)
    except (TypeError, ValueError):
        return "unknown"
    if lvl == 3:
        return "green"
    if lvl == 2:
        return "amber"
    if lvl == 1:
        return "red"
    if lvl == 0:
        return "grey"
    return "unknown"


def build_genome_annotation(entry: Dict) -> str:
    """
    Build a single string summarising Ensembl IDs and locations for GRCh37 and GRCh38.

    Arguments:
        - entry: A single gene entry from PanelApp API response.

    Format examples:
      GRCh37:ENSG00000167207|16:50727514-50766988; GRCh38:ENSG00000167207|16:50693603-50734041
    Multiple versions per build are included as additional '; '-separated fragments.
    """
    gene_data = entry.get("gene_data") or {}
    ens = gene_data.get("ensembl_genes") or {}
    parts: List[str] = []

    # Handle both GRch37/GRCh37 and GRch38/GRCh38 key spellings
    build_key_map = {
        "GRch37": "GRCh37",
        "GRCh37": "GRCh37",
        "GRch38": "GRCh38",
        "GRCh38": "GRCh38",
    }

    for key, label in build_key_map.items():
        build_info = ens.get(key)
        if not isinstance(build_info, dict):
            continue
        for ver, vinfo in build_info.items():
            if not isinstance(vinfo, dict):
                continue
            eid = vinfo.get("ensembl_id") or ""
            loc = vinfo.get("location") or ""
            if not eid and not loc:
                continue
            parts.append(f"{label}:{eid}|{loc}")

    return "; ".join(parts)


def fetch_gene_entries(
    gene_symbol: str, base_url: str, timeout: float = 20.0
) -> Optional[List[Dict]]:
    """
    Query PanelApp gene search endpoint and return a list of entries
    (one per gene–panel association).

    Arguments:
        - gene_symbol: HGNC gene symbol to search for.
        - base_url: Base URL for PanelApp API.
        - timeout: Request timeout in seconds (default: 20.0).

    Returns:
        list of dicts for each panel membership, or None if gene not found.
    """
    url = f"{base_url.rstrip('/')}/genes/"
    try:
        resp = requests.get(
            url,
            params={"entity_name": gene_symbol, "format": "json"},
            timeout=timeout,
        )
    except requests.RequestException as e:
        sys.stderr.write(f"[WARN] Request error for {gene_symbol}: {e}\n")
        return None

    if not resp.ok:
        sys.stderr.write(
            f"[WARN] Non-OK status for {gene_symbol}: {resp.status_code} {resp.text[:200]}\n"
        )
        return None

    data = resp.json()
    count = data.get("count", 0)
    if count == 0:
        # Gene not present on any panel
        return None

    results = data.get("results") or []
    if not results:
        return None

    rows: List[Dict] = []
    for entry in results:
        panel = entry.get("panel") or {}
        confidence_level = entry.get("confidence_level")
        color = confidence_to_color(confidence_level)

        phenotypes = entry.get("phenotypes") or []
        evidence_list = entry.get("evidence") or []
        relevant_disorders = panel.get("relevant_disorders") or []

        genome_annotation = build_genome_annotation(entry)

        rows.append(
            {
                "Gene": gene_symbol,
                "Genome_Annotation": genome_annotation,
                "Panel_ID": panel.get("id"),
                "Panel_Name": panel.get("name"),
                "Disease_Group": panel.get("disease_group"),
                "Disease_Subgroup": panel.get("disease_sub_group"),
                "Confidence_Level": confidence_level,
                "color": color,
                "Phenotypes": "; ".join(phenotypes),
                "mode_of_inheritance": entry.get("mode_of_inheritance"),
                "Evidence": "; ".join(evidence_list),
                "relevant_disorders": "; ".join(relevant_disorders),
            }
        )

    return rows if rows else None


def annotate_genes(
    genes: List[str],
    base_url: str,
    output_csv: str,
    sleep_sec: float = 0.2,
    timeout: float = 20.0,
) -> None:
    """
    For each gene, query PanelApp and write per-panel rows to output CSV.
    If a gene is not registered (no panels), a single row with Is_Registered = False is written.
    Arguments:
    - genes: List of gene symbols to query.
    - base_url: Base URL for PanelApp API.
    - output_csv: Path to output CSV file.
    - sleep_sec: Seconds to sleep between API calls (default: 0.2).
    - timeout: Request timeout in seconds (default: 20.0).
    """
    fieldnames = [
        "Gene",
        "Is_Registered",
        "Genome_Annotation",
        "Panel_ID",
        "Panel_Name",
        "Disease_Group",
        "Disease_Subgroup",
        "Confidence_Level",
        "color",
        "Phenotypes",
        "mode_of_inheritance",
        "Evidence",
        "relevant_disorders",
    ]

    with open(output_csv, "w", newline="", encoding="utf-8") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=fieldnames)
        writer.writeheader()

        for i, gene in enumerate(genes, start=1):
            sys.stderr.write(f"[INFO] ({i}/{len(genes)}) Querying {gene}\n")
            rows = fetch_gene_entries(
                gene_symbol=gene,
                base_url=base_url,
                timeout=timeout,
            )

            if not rows:
                # Gene not registered on any panel
                writer.writerow(
                    {
                        "Gene": gene,
                        "Is_Registered": False,
                        "Genome_Annotation": "",
                        "Panel_ID": "",
                        "Panel_Name": "",
                        "Disease_Group": "",
                        "Disease_Subgroup": "",
                        "Confidence_Level": "",
                        "color": "",
                        "Phenotypes": "",
                        "mode_of_inheritance": "",
                        "Evidence": "",
                        "relevant_disorders": "",
                    }
                )
            else:
                # One row per gene–panel association
                for row in rows:
                    row_out = {
                        **row,
                        "Is_Registered": True,
                    }
                    writer.writerow(row_out)

            if sleep_sec > 0:
                time.sleep(sleep_sec)


def main(argv: Optional[List[str]] = None) -> None:
    args = parse_args(argv)
    genes = read_genes_from_csv(args.input_csv, args.gene_column)
    annotate_genes(
        genes=genes,
        base_url=args.base_url,
        output_csv=args.output_csv,
        sleep_sec=args.sleep,
        timeout=args.timeout,
    )


if __name__ == "__main__":
    main()
