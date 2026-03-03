#!/usr/bin/env python3
"""Quantify reads flanked by R1-like sequence at both ends.

This script searches both ends of each read for local alignments to a primer
sequence and its reverse-complement, allowing mismatches and truncation.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import Dict, Generator, Tuple

import pandas as pd
from Bio import pairwise2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Quantify reads flanked by R1-like sequence on both ends."
    )
    parser.add_argument("--reads", required=True, help="Input FASTQ or FASTQ.GZ")
    parser.add_argument("--primer", required=True, help="Primer sequence (5'->3')")
    parser.add_argument(
        "--window",
        type=int,
        default=80,
        help="Number of bases to inspect at each read end (default: 80)",
    )
    parser.add_argument(
        "--min-aln-len",
        type=int,
        default=12,
        help="Minimum aligned primer bases for a side-hit (default: 12)",
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.80,
        help="Minimum identity for a side-hit (default: 0.80)",
    )
    parser.add_argument(
        "--out-prefix",
        required=True,
        help="Output prefix (writes .per_read.tsv, .summary.tsv, .threshold_grid.tsv)",
    )
    return parser.parse_args()


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def open_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("r")


def iter_fastq(path: Path) -> Generator[Tuple[str, str], None, None]:
    with open_maybe_gzip(path) as fh:
        while True:
            header = fh.readline().strip()
            if not header:
                break
            seq = fh.readline().strip()
            fh.readline()  # plus
            fh.readline()  # qual
            yield header[1:].split()[0], seq.upper()


def alignment_stats(query: str, motif: str, side_name: str, orientation: str) -> Dict[str, object]:
    # Local alignment on a small window; tuned to tolerate errors but keep
    # contiguous signal near primer-like regions.
    aln = pairwise2.align.localms(
        query, motif, 2.0, -1.0, -2.0, -0.5, one_alignment_only=True
    )
    if not aln:
        return {
            f"{side_name}_orientation": orientation,
            f"{side_name}_score": 0.0,
            f"{side_name}_aln_primer_bases": 0,
            f"{side_name}_matches": 0,
            f"{side_name}_mismatches": 0,
            f"{side_name}_identity": 0.0,
            f"{side_name}_primer_coverage": 0.0,
            f"{side_name}_window_start": -1,
            f"{side_name}_window_end": -1,
        }

    seq_a, seq_b, score, start, end = aln[0]
    aln_primer_bases = 0
    matches = 0
    mismatches = 0

    for a, b in zip(seq_a, seq_b):
        if b == "-":
            continue
        aln_primer_bases += 1
        if a == b:
            matches += 1
        elif a != "-":
            mismatches += 1
        else:
            # gap in query relative to primer counts as mismatch-like event
            mismatches += 1

    identity = (matches / aln_primer_bases) if aln_primer_bases else 0.0
    primer_coverage = aln_primer_bases / len(motif)

    return {
        f"{side_name}_orientation": orientation,
        f"{side_name}_score": float(score),
        f"{side_name}_aln_primer_bases": int(aln_primer_bases),
        f"{side_name}_matches": int(matches),
        f"{side_name}_mismatches": int(mismatches),
        f"{side_name}_identity": float(identity),
        f"{side_name}_primer_coverage": float(primer_coverage),
        f"{side_name}_window_start": int(start),
        f"{side_name}_window_end": int(end),
    }


def best_side_match(window_seq: str, primer: str, primer_rc: str, side_name: str) -> Dict[str, object]:
    best = None
    for orientation, motif in (("R1", primer), ("R1_rc", primer_rc)):
        stats = alignment_stats(window_seq, motif, side_name, orientation)
        if best is None:
            best = stats
            continue
        # Tie-break by score, then aligned primer bases, then identity.
        key_new = (
            stats[f"{side_name}_score"],
            stats[f"{side_name}_aln_primer_bases"],
            stats[f"{side_name}_identity"],
        )
        key_old = (
            best[f"{side_name}_score"],
            best[f"{side_name}_aln_primer_bases"],
            best[f"{side_name}_identity"],
        )
        if key_new > key_old:
            best = stats
    return best


def main() -> None:
    args = parse_args()
    reads_path = Path(args.reads)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    primer = args.primer.upper()
    primer_rc = revcomp(primer)

    rows = []
    for read_id, seq in iter_fastq(reads_path):
        left_window = seq[: args.window]
        right_window = seq[-args.window :] if len(seq) > args.window else seq

        left = best_side_match(left_window, primer, primer_rc, "left")
        right = best_side_match(right_window, primer, primer_rc, "right")

        row = {
            "read_id": read_id,
            "read_length": len(seq),
            **left,
            **right,
        }
        rows.append(row)

    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError("No reads found in FASTQ input.")

    df["left_hit"] = (df["left_aln_primer_bases"] >= args.min_aln_len) & (
        df["left_identity"] >= args.min_identity
    )
    df["right_hit"] = (df["right_aln_primer_bases"] >= args.min_aln_len) & (
        df["right_identity"] >= args.min_identity
    )
    df["flanked"] = df["left_hit"] & df["right_hit"]

    orient = []
    for _, row in df.iterrows():
        if row["flanked"]:
            orient.append(f"{row['left_orientation']}|{row['right_orientation']}")
        else:
            orient.append("NA")
    df["flanked_orientation"] = orient

    per_read_path = Path(f"{out_prefix}.per_read.tsv")
    df.to_csv(per_read_path, sep="\t", index=False)

    total = len(df)
    left_n = int(df["left_hit"].sum())
    right_n = int(df["right_hit"].sum())
    flanked_n = int(df["flanked"].sum())
    summary = pd.DataFrame(
        [
            {"metric": "total_reads", "value": total},
            {"metric": "left_hit_reads", "value": left_n},
            {"metric": "left_hit_pct", "value": left_n / total},
            {"metric": "right_hit_reads", "value": right_n},
            {"metric": "right_hit_pct", "value": right_n / total},
            {"metric": "flanked_reads", "value": flanked_n},
            {"metric": "flanked_pct", "value": flanked_n / total},
            {"metric": "window", "value": args.window},
            {"metric": "min_aln_len", "value": args.min_aln_len},
            {"metric": "min_identity", "value": args.min_identity},
        ]
    )
    summary_path = Path(f"{out_prefix}.summary.tsv")
    summary.to_csv(summary_path, sep="\t", index=False)

    # Sensitivity table across thresholds, computed from same alignment calls.
    grid_rows = []
    for min_len in (10, 12, 14, 16, 18, 20):
        for min_ident in (0.70, 0.75, 0.80, 0.85, 0.90, 0.95):
            left_hit = (df["left_aln_primer_bases"] >= min_len) & (
                df["left_identity"] >= min_ident
            )
            right_hit = (df["right_aln_primer_bases"] >= min_len) & (
                df["right_identity"] >= min_ident
            )
            flanked = left_hit & right_hit
            n = int(flanked.sum())
            grid_rows.append(
                {
                    "min_aln_len": min_len,
                    "min_identity": min_ident,
                    "flanked_reads": n,
                    "flanked_pct": n / total,
                }
            )
    grid = pd.DataFrame(grid_rows)
    grid_path = Path(f"{out_prefix}.threshold_grid.tsv")
    grid.to_csv(grid_path, sep="\t", index=False)

    orient_counts = (
        df.loc[df["flanked"], "flanked_orientation"].value_counts().rename_axis("orientation").reset_index(name="count")
    )
    orient_counts["pct_of_all_reads"] = orient_counts["count"] / total
    orient_counts["pct_of_flanked_reads"] = orient_counts["count"] / max(flanked_n, 1)
    orient_path = Path(f"{out_prefix}.orientation_counts.tsv")
    orient_counts.to_csv(orient_path, sep="\t", index=False)

    print(f"Wrote per-read calls: {per_read_path}")
    print(f"Wrote summary: {summary_path}")
    print(f"Wrote threshold grid: {grid_path}")
    print(f"Wrote orientation counts: {orient_path}")
    print(
        f"Flanked reads (default thresholds): {flanked_n}/{total} "
        f"({(flanked_n/total)*100:.2f}%)"
    )


if __name__ == "__main__":
    main()
