#!/usr/bin/env python3
"""
Author: Diego Fuentes
Contact email: diegofupa@gmail.com
Barcelona
Date:2026-02-04
"""

import argparse
import logging
import sys
from pathlib import Path
import pandas as pd
import re

# Configure logging
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


# ============================================================================
# FILE PARSING FUNCTIONS
# ============================================================================


def parse_bioinfo_result_file(filepath: str | Path) -> pd.DataFrame:
    """
    Parser for biomedical research result files.
    Handles metadata, variable spacing, and structured entity listings.

    Parameters
    ----------
    filepath : str or Path
        Path to the result file (.txt format)

    Returns
    -------
    pd.DataFrame
        Parsed data with columns:
        - Entity Type (Gene/Variant)
        - Entity Name
        - Associated Gene
        - Source
        - Perplexity Status
        - Notes
    """

    filepath = Path(filepath)

    if not filepath.exists():
        logger.error(f"File not found: {filepath}")
        raise FileNotFoundError(f"File not found: {filepath}")

    logger.info(f"Parsing file: {filepath.name}")

    with open(filepath, "r", encoding="utf-8") as f:
        content = f.read()

    lines = content.split("\n")

    # Find data start marker
    data_start = 0
    for i, line in enumerate(lines):
        if "Entity Type" in line:
            data_start = i + 1
            logger.debug(f"Found data start at line {i}")
            break

    # Find data end marker
    data_end = len(lines)
    for i in range(data_start, len(lines)):
        if lines[i].strip().startswith("{") or "Total entities:" in lines[i]:
            data_end = i
            logger.debug(f"Found data end at line {i}")
            break

    rows = []

    for line_num, line in enumerate(lines[data_start:data_end], start=data_start):
        if not line.strip() or "====" in line:
            continue

        if not re.match(r"^\s+(Gene|Variant)", line):
            continue

        # Extract entity type
        entity_type_match = re.search(r"(Gene|Variant)", line)
        entity_type = entity_type_match.group(1)

        rest_of_line = line[entity_type_match.end() :].lstrip()

        # Extract status
        status_pattern = r"(Confirmed|Rejected|Uncertain|Not\s+assessed)"
        status_match = re.search(status_pattern, rest_of_line)

        if not status_match:
            logger.debug(f"No status found in line {line_num}: {line[:50]}...")
            continue

        status = status_match.group(1).strip()

        # Extract source
        source_pattern = r"(Both|Abstract\s+Only|Full\s+Text\s+PMC|ClinVar\s+Only|Full\s+Text\s+Only)"
        source_search = rest_of_line[: status_match.start()]
        source_match = re.search(source_pattern, source_search)

        if source_match:
            source = source_match.group(1).strip()
            entity_and_gene = rest_of_line[: source_match.start()].strip()
        else:
            source = "Unknown"
            entity_and_gene = source_search.strip()

        # Extract entity name and associated gene
        parts = entity_and_gene.rsplit(maxsplit=1)

        if len(parts) == 2:
            entity_name, associated_gene = parts
        elif len(parts) == 1:
            entity_name = parts[0]
            associated_gene = ""
        else:
            logger.debug(f"Could not parse entity/gene on line {line_num}")
            continue

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
    logger.info(f"Successfully parsed {len(df)} entities from {filepath.name}")

    return df


def normalize_gene_field(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize gene field: for Variants use Associated Gene, for Genes use Entity Name.
    Creates 'Gene_Normalized' column for consistent analysis.

    Parameters
    ----------
    df : pd.DataFrame
        Parsed dataframe from parse_bioinfo_result_file()

    Returns
    -------
    pd.DataFrame
        Dataframe with added 'Gene_Normalized' column
    """

    df = df.copy()

    df["Gene_Normalized"] = df.apply(
        lambda row: (
            row["Associated Gene"]
            if row["Entity Type"] == "Variant" and row["Associated Gene"]
            else row["Entity Name"]
        ),
        axis=1,
    )

    logger.debug(
        f"Gene field normalized: {df['Gene_Normalized'].nunique()} unique genes"
    )

    return df


# ============================================================================
# ANALYSIS FUNCTIONS
# ============================================================================


def analyze_file_statistics(df: pd.DataFrame, filename: str = "") -> dict:
    """
    Generate summary statistics for parsed file.

    Parameters
    ----------
    df : pd.DataFrame
        Parsed dataframe
    filename : str
        Original filename for logging

    Returns
    -------
    dict
        Summary statistics
    """

    stats = {
        "filename": filename,
        "total_entities": len(df),
        "genes": len(df[df["Entity Type"] == "Gene"]),
        "variants": len(df[df["Entity Type"] == "Variant"]),
        "confirmed": len(df[df["Perplexity Status"] == "Confirmed"]),
        "uncertain": len(df[df["Perplexity Status"] == "Uncertain"]),
        "rejected": len(df[df["Perplexity Status"] == "Rejected"]),
        "not_assessed": len(df[df["Perplexity Status"] == "Not assessed"]),
        "unique_genes": df["Gene_Normalized"].nunique(),
        "sources": df["Source"].value_counts().to_dict(),
    }

    return stats


def print_file_summary(stats: dict) -> None:
    """Print nicely formatted summary statistics."""

    print("\n" + "=" * 80)
    print(f"FILE SUMMARY: {Path(stats['filename']).name}")
    print("=" * 80)

    print(f"\nEntity Types:")
    print(f"  Total entities:      {stats['total_entities']}")
    print(f"  Genes:               {stats['genes']}")
    print(f"  Variants:            {stats['variants']}")

    print(f"\nStatus Distribution:")
    print(f"  Confirmed:           {stats['confirmed']}")
    print(f"  Uncertain:           {stats['uncertain']}")
    print(f"  Rejected:            {stats['rejected']}")
    print(f"  Not assessed:        {stats['not_assessed']}")

    print(f"\nGene Analysis:")
    print(f"  Unique genes:        {stats['unique_genes']}")

    print(f"\nSources:")
    for source, count in sorted(
        stats["sources"].items(), key=lambda x: x[1], reverse=True
    ):
        pct = 100 * count / stats["total_entities"]
        print(f"  {source:30s}: {count:4d} ({pct:5.1f}%)")

    print()


# ============================================================================
# MAIN ARGUMENT PARSER & EXECUTION
# ============================================================================


def main():
    """Main script execution function with argument parsing."""

    parser = argparse.ArgumentParser(
        description="Parse biomedical research result files (.txt format) and perform gene analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Parse single file
  python parse_bioinfo_files.py spa_results.txt
  
  # Parse single file with output
  python parse_bioinfo_files.py hs_results.txt --output parsed_hs.csv
  
  # Parse multiple files
  python parse_bioinfo_files.py *.txt --output-dir parsed_results/
  
  # Parse with detailed logging
  python parse_bioinfo_files.py results.txt --verbose
  
  # Filter for Confirmed+Uncertain only
  python parse_bioinfo_files.py results.txt --filter-status Confirmed Uncertain
        """,
    )

    parser.add_argument(
        "files",
        nargs="+",
        type=str,
        help="Input file(s) to parse (.txt format). Accepts wildcards (e.g., *.txt)",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=None,
        help="Output CSV filename (for single file input)",
    )

    parser.add_argument(
        "-d",
        "--output-dir",
        type=str,
        default="parsed_results/",
        help="Output directory for results (default: parsed_results/)",
    )

    parser.add_argument(
        "-f",
        "--filter-status",
        nargs="+",
        choices=["Confirmed", "Uncertain", "Rejected", "Not assessed"],
        default=None,
        help="Filter entities by status (e.g., --filter-status Confirmed Uncertain)",
    )

    parser.add_argument(
        "-s",
        "--summary-only",
        action="store_true",
        help="Only generate summary statistics (no CSV output)",
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging"
    )

    args = parser.parse_args()

    # Configure logging
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")

    # Expand file patterns
    input_files = []
    for file_pattern in args.files:
        matched = list(Path(".").glob(file_pattern))
        if not matched:
            logger.warning(f"No files matched pattern: {file_pattern}")
        else:
            input_files.extend(matched)

    if not input_files:
        logger.error("No input files found")
        sys.exit(1)

    input_files = [f for f in input_files if f.is_file()]
    logger.info(f"Found {len(input_files)} file(s) to process")

    # Setup output directory
    output_dir = Path(args.output_dir)
    if not args.summary_only:
        output_dir.mkdir(exist_ok=True)
        logger.info(f"Output directory: {output_dir}")

    # Process files
    all_stats = []
    all_dfs = {}

    for input_file in sorted(input_files):
        try:
            # Parse file
            df = parse_bioinfo_result_file(input_file)

            # Normalize gene field
            df = normalize_gene_field(df)

            # Filter by status if requested
            if args.filter_status:
                original_count = len(df)
                df = df[df["Perplexity Status"].isin(args.filter_status)]
                logger.info(
                    f"Filtered to {len(df)} entities ({args.filter_status}) "
                    f"from {original_count} total"
                )

            # Generate statistics
            stats = analyze_file_statistics(df, str(input_file))
            all_stats.append(stats)
            all_dfs[input_file.stem] = df

            # Print summary
            print_file_summary(stats)

            # Save CSV if not summary-only
            if not args.summary_only:
                if args.output and len(input_files) == 1:
                    output_file = Path(args.output)
                else:
                    # Auto-generate output filename
                    status_suffix = (
                        "_" + "_".join(args.filter_status).lower()
                        if args.filter_status
                        else ""
                    )
                    output_file = (
                        output_dir / f"{input_file.stem}_parsed{status_suffix}.csv"
                    )

                df.to_csv(output_file, index=False)
                logger.info(f"Saved to: {output_file}")

        except Exception as e:
            logger.error(f"Error processing {input_file}: {e}", exc_info=args.verbose)
            continue

    # Generate master summary if multiple files
    if len(all_stats) > 1:
        print("\n" + "=" * 80)
        print("MASTER SUMMARY - ALL FILES")
        print("=" * 80)

        summary_df = pd.DataFrame(all_stats)

        # Display summary table
        summary_display = summary_df[
            [
                "filename",
                "total_entities",
                "genes",
                "variants",
                "confirmed",
                "uncertain",
                "rejected",
                "unique_genes",
            ]
        ].copy()

        summary_display["filename"] = summary_display["filename"].apply(
            lambda x: Path(x).name
        )

        print("\n" + summary_display.to_string(index=False))

        # Save master summary
        summary_file = output_dir / "master_summary.csv"
        summary_df.to_csv(summary_file, index=False)
        logger.info(f"Master summary saved to: {summary_file}")

        # Create merged dataset if requested
        if len(all_dfs) > 1 and not args.summary_only:
            merged_df = pd.concat(
                [
                    df.assign(disease_key=disease_key)
                    for disease_key, df in all_dfs.items()
                ],
                ignore_index=True,
            )

            merged_file = output_dir / "merged_all_files.csv"
            merged_df.to_csv(merged_file, index=False)

            logger.info(
                f"Merged dataset saved to: {merged_file} ({len(merged_df)} rows)"
            )
            print(f"\n✓ Merged all files: {merged_file}")

    logger.info(f"✓ Processing complete. Processed {len(all_stats)} file(s)")


if __name__ == "__main__":
    main()
