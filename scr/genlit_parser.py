#!/usr/bin/env python3
"""
Genlit Parser - Convert .txt results to .csv

Usage:
    python parse_genlit.py input_file.txt [output_file.csv]

Arguments:
    input_file.txt     Input Genlit results file (.txt)
    output_file.csv    Output CSV file (optional, defaults to input_file.csv)

Example:
    python parse_genlit.py hidradenitis_results.txt
    python parse_genlit.py hs_results.txt hidradenitis.csv
"""

import argparse
import sys
import pandas as pd
import re
from pathlib import Path


def parse_bioinfo_result_file(filepath):
    """
    Generalized parser for Genlit .txt result files.

    Designed to handle the formatted .txt files with:
    - Disease/date metadata header
    - Fixed-width or loosely-aligned tabular columns
    - Variable-length text fields (Notes column)
    - JSON payload at end of file (automatically ignored)

    Parameters
    ----------
    filepath : str or Path
        Path to the .txt results file

    Returns
    -------
    pd.DataFrame
        DataFrame with columns:
        - Entity Type: "Gene" or "Variant"
        - Entity Name: The gene/variant identifier
        - Associated Gene: Gene association (if applicable)
        - Source: Where entity was found (e.g., "Both", "Abstract Only", "ClinVar Only", "Full Text PMC")
        - Perplexity Status: Classification ("Confirmed", "Uncertain", "Rejected", "Not assessed")
        - Notes: Description/explanation text
    """

    with open(filepath, "r", encoding="utf-8") as f:
        content = f.read()

    lines = content.split("\n")

    # Find where data section starts (after header with "Entity Type")
    data_start = 0
    for i, line in enumerate(lines):
        if "Entity Type" in line:
            data_start = i + 1
            break

    # Find where data ends (before JSON object or summary section)
    data_end = len(lines)
    for i in range(data_start, len(lines)):
        if lines[i].strip().startswith("{") or "Total entities:" in lines[i]:
            data_end = i
            break

    # Parse data rows
    rows = []

    for line in lines[data_start:data_end]:
        # Skip separators and empty lines
        if not line.strip() or "====" in line:
            continue

        # Valid data rows start with whitespace followed by entity type
        if not re.match(r"^\s+(Gene|Variant)", line):
            continue

        # Extract entity type
        entity_type_match = re.search(r"(Gene|Variant)", line)
        entity_type = entity_type_match.group(1)

        # Parse remainder of line
        rest_of_line = line[entity_type_match.end() :].lstrip()

        # Find status (anchor point for parsing)
        status_pattern = r"(Confirmed|Rejected|Uncertain|Not\s+assessed)"
        status_match = re.search(status_pattern, rest_of_line)

        if not status_match:
            continue

        status = status_match.group(1).strip()

        # Find source (search backwards from status)
        source_pattern = r"(Both|Abstract\s+Only|Full\s+Text\s+PMC|ClinVar\s+Only|Full\s+Text\s+Only)"
        source_search = rest_of_line[: status_match.start()]
        source_match = re.search(source_pattern, source_search)

        if source_match:
            source = source_match.group(1).strip()
            entity_and_gene = rest_of_line[: source_match.start()].strip()
        else:
            source = "Unknown"
            entity_and_gene = source_search.strip()

        # Split entity name from associated gene
        parts = entity_and_gene.rsplit(maxsplit=1)

        if len(parts) == 2:
            entity_name, associated_gene = parts
        elif len(parts) == 1:
            entity_name = parts[0]
            associated_gene = ""
        else:
            continue

        # Notes = everything after status
        notes = rest_of_line[status_match.end() :].strip()

        rows.append(
            {
                "Entity Type": entity_type,
                "Entity Name": entity_name.strip(),
                "Associated Gene": associated_gene.strip(),
                "Source": source,
                "Perplexity Status": status,
                "Notes": notes,
            }
        )

    df = pd.DataFrame(rows)
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Convert Genlit .txt results to structured CSV format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s hidradenitis_results.txt
  %(prog)s hs_results.txt hidradenitis_results.csv
  %(prog)s spa_results.txt --verbose
        """,
    )

    parser.add_argument("input_file", type=str, help="Input Genlit results file (.txt)")

    parser.add_argument(
        "output_file",
        nargs="?",
        type=str,
        default=None,
        help="Output CSV file (optional, defaults to input_file.csv)",
    )

    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Print detailed parsing information",
    )

    args = parser.parse_args()

    # Input validation
    input_path = Path(args.input_file)

    if not input_path.exists():
        print(f"[ERROR]: Input file '{args.input_file}' not found!")
        sys.exit(1)

    if not input_path.suffix == ".txt":
        print("[WARNING]: Input file should have .txt extension")

    # Determine output path
    if args.output_file:
        output_path = Path(args.output_file)
    else:
        output_path = input_path.with_suffix(".csv")

    if args.verbose:
        print(f"Input file:  {input_path.absolute()}")
        print(f"Output file: {output_path.absolute()}")
        print()

    try:
        # Parse the file
        df = parse_bioinfo_result_file(input_path)

        # Save to CSV
        df.to_csv(output_path, index=False)

        # Summary
        total_entities = len(df)
        status_counts = df["Perplexity Status"].value_counts()
        entity_types = df["Entity Type"].value_counts()

        print(f"Successfully parsed {total_entities} entities!")
        print(f"Breakdown:")
        print(f"   Genes:     {entity_types.get('Gene', 0)}")
        print(f"   Variants:  {entity_types.get('Variant', 0)}")
        print(f"   Confirmed: {status_counts.get('Confirmed', 0)}")
        print(f"   Uncertain: {status_counts.get('Uncertain', 0)}")
        print(f"   Rejected:  {status_counts.get('Rejected', 0)}")
        print(f"   Not assessed: {status_counts.get('Not assessed', 0)}")
        print(f"Results saved to: {output_path.absolute()}")

        if args.verbose:
            print("\nPreview of first 5 rows:")
            print(df.head().to_string())

    except Exception as e:
        print(f"Error during parsing: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
