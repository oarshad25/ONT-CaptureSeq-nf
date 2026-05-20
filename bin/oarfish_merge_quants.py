#!/usr/bin/env python3

"""Merge per-sample Oarfish .quant files into a transcript x sample matrix."""

import argparse
import csv
import sys
from pathlib import Path


TRANSCRIPT_COLUMN_CANDIDATES = (
    "target_id",
    "target",
    "transcript_id",
    "transcript",
    "tx_id",
    "txp",
    "name",
    "Name",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge Oarfish .quant files into a transcript x sample count matrix."
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output transcript x sample TSV path.",
    )
    parser.add_argument(
        "--count-column",
        default="num_reads",
        help="Oarfish count column to merge. Default: num_reads.",
    )
    parser.add_argument(
        "--transcript-column",
        default=None,
        help="Transcript identifier column. If omitted, common names are tried before falling back to the first column.",
    )
    parser.add_argument(
        "quant_files",
        nargs="+",
        help="Per-sample Oarfish .quant files.",
    )
    return parser.parse_args()


def sample_id_from_path(path):
    name = Path(path).name
    return name[:-6] if name.endswith(".quant") else Path(path).stem


def find_transcript_column(fieldnames, requested):
    if requested:
        if requested not in fieldnames:
            raise ValueError(
                f"Transcript column '{requested}' not found. Available columns: {', '.join(fieldnames)}"
            )
        return requested

    for candidate in TRANSCRIPT_COLUMN_CANDIDATES:
        if candidate in fieldnames:
            return candidate

    if not fieldnames:
        raise ValueError("No header columns found")

    return fieldnames[0]


def read_quant(path, transcript_column, count_column):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path}: empty file or missing header")

        if count_column not in reader.fieldnames:
            raise ValueError(
                f"{path}: count column '{count_column}' not found. Available columns: {', '.join(reader.fieldnames)}"
            )

        tx_col = find_transcript_column(reader.fieldnames, transcript_column)
        counts = {}
        for row_number, row in enumerate(reader, start=2):
            transcript_id = row.get(tx_col, "")
            if transcript_id == "":
                raise ValueError(f"{path}: missing transcript id on line {row_number}")
            if transcript_id in counts:
                raise ValueError(f"{path}: duplicate transcript id '{transcript_id}'")
            counts[transcript_id] = row.get(count_column, "")

    return counts


def main():
    args = parse_args()

    sample_to_counts = {}
    for quant_file in args.quant_files:
        sample_id = sample_id_from_path(quant_file)
        if sample_id in sample_to_counts:
            raise ValueError(f"Duplicate sample id '{sample_id}' from {quant_file}")
        sample_to_counts[sample_id] = read_quant(
            quant_file,
            args.transcript_column,
            args.count_column,
        )

    samples = sorted(sample_to_counts)
    transcripts = sorted(
        {
            transcript_id
            for counts in sample_to_counts.values()
            for transcript_id in counts
        }
    )

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["transcript_id"] + samples)
        for transcript_id in transcripts:
            writer.writerow(
                [transcript_id]
                + [sample_to_counts[sample].get(transcript_id, "") for sample in samples]
            )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        sys.exit(f"ERROR: {error}")
