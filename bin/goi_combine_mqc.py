#!/usr/bin/env python3
# bin/goi_combine_mqc.py

"""
goi_combine_mqc.py

Merges per-GOI metric TSVs from goi_readcount_n50.py and goi_fulllength_count.py
into a single MultiQC-compatible table.

Each input TSV is wide-format with one row per sample. This script merges
them horizontally by sample and prepends the MultiQC custom content header,
producing a table section in the MultiQC report.

Collision detection ensures that duplicate column names across input files
(which would indicate non-unique gene names in the GOI BED) cause an early
error rather than silently overwriting values.

Usage:
    goi_combine_mqc.py --readcount-n50 sample_goi_readcount_n50.tsv
                       --fulllength sample_goi_fulllength.tsv
                       --sample SAMPLE_ID
                       --output sample_goi_summary_mqc.tsv
"""

import argparse
import csv
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine per-GOI metric TSVs into a single MultiQC table"
    )
    parser.add_argument("--readcount-n50", required=True, dest="readcount_n50",
                        help="Output of goi_readcount_n50.py")
    parser.add_argument("--fulllength",    required=True,
                        help="Output of goi_fulllength_count.py")
    parser.add_argument("--sample",        required=True,
                        help="Sample name")
    parser.add_argument("--output",        required=True,
                        help="Output MultiQC TSV path")
    return parser.parse_args()


# MultiQC custom content header template.
MQC_HEADER_TEMPLATE = """\
# id: 'goi_summary'
# section_name: 'Genes of Interest — Summary'
# description: 'Read count, percentage of total primary mapped reads, read N50, and transcript-boundary-spanning read count per gene of interest locus.'
# format: 'tsv'
# plot_type: 'table'
# pconfig:
#   id: 'goi_summary_table'
#   title: 'Genes of Interest Summary'
"""


def write_column_headers(out, cols):
    """
    Write MultiQC table header config in the same order as the TSV columns.

    MultiQC tables default to one decimal place, which hides small
    percentages. Defining all headers preserves the input column order.
    """
    out.write("# headers:\n")
    for col in cols:
        out.write(f"#   {col}:\n")
        out.write(f"#     title: '{col}'\n")
        if col.startswith("pct_total_primary_mapped_"):
            out.write("#     format: '{:,.4f}'\n")
            out.write("#     suffix: '%'\n")
        else:
            out.write("#     format: '{:,.0f}'\n")


def read_wide_tsv(path):
    """
    Read a wide-format TSV file and return the single data row as a dict.

    Each per-sample output TSV from the GOI metric scripts contains exactly
    one data row (one row per sample). Raises ValueError if this expectation
    is not met, catching truncated or malformed files early.

    Args:
        path: Path to the TSV file.

    Returns:
        Dict of column name -> value for the single data row.

    Raises:
        ValueError: If the file does not contain exactly one data row.
    """
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows   = list(reader)
    if len(rows) != 1:
        raise ValueError(
            f"Expected exactly 1 data row in {path}, got {len(rows)}."
        )
    return rows[0]


def validate_sample(row, path, expected_sample):
    """
    Warn if the sample name embedded in an input TSV does not match the
    expected sample name passed via --sample.

    A mismatch indicates that the wrong file has been passed for this sample,
    which could happen if Nextflow staging produces unexpected filename
    collisions.

    Args:
        row:             Dict of the data row from the TSV.
        path:            Path to the TSV (used in the warning message).
        expected_sample: Sample name from the --sample argument.
    """
    if row.get("Sample") != expected_sample:
        print(
            f"Warning: sample name in {path} ({row.get('Sample')}) "
            f"does not match expected sample ({expected_sample})",
            file=sys.stderr
        )


def check_collisions(named_dicts):
    """
    Check for duplicate column names across the input TSV dicts.

    Duplicate column names would arise if two GOIs share the same name in
    the GOI BED file, causing silent value overwriting during dict merge.
    Raises ValueError with a diagnostic message identifying which files
    share the conflicting column.

    Args:
        named_dicts: List of (name, dict) tuples where name is a label for
                     the input file and dict is the data row.

    Raises:
        ValueError: If any non-Sample column name appears in more than one input.
    """
    seen = {}
    for name, d in named_dicts:
        for col in d:
            if col == "Sample":
                continue
            if col in seen:
                raise ValueError(
                    f"Column '{col}' appears in both {seen[col]} and {name}. "
                    "Check that gene names are unique in your GOI BED file."
                )
            seen[col] = name


def main():
    args = parse_args()

    inputs = {
        "readcount_n50": read_wide_tsv(args.readcount_n50),
        "fulllength":    read_wide_tsv(args.fulllength),
    }

    # Validate sample names match across all input files
    for name, row in inputs.items():
        validate_sample(row, name, args.sample)

    # Check for column name collisions before merging
    check_collisions(list(inputs.items()))

    # Merge inputs — Sample column taken from readcount_n50, excluded from others
    merged = {**inputs["readcount_n50"]}
    for name in ["fulllength"]:
        merged.update({k: v for k, v in inputs[name].items() if k != "Sample"})

    cols = [c for c in merged if c != "Sample"]

    with open(args.output, "w") as out:
        out.write(MQC_HEADER_TEMPLATE)
        write_column_headers(out, cols)
        out.write("Sample\t"    + "\t".join(cols) + "\n")
        out.write(args.sample + "\t" + "\t".join(merged[c] for c in cols) + "\n")


if __name__ == "__main__":
    main()
