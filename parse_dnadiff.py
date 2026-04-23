#!/usr/bin/env python3
"""
parse_dnadiff.py

Parse a MUMmer4 show-diff output file and extract candidate insertion
sequences using samtools faidx.

Usage:
    python parse_dnadiff.py -d differences.txt -a assembly.fna [options]

Options:
    -d, --diff        Path to show-diff output file (required)
    -a, --assembly    Path to assembly FASTA file (required)
    -o, --output      Output FASTA file for extracted sequences [default: candidates.fna]
    -t, --types       Feature types to extract, comma-separated [default: GAP,DUP]
    -m, --min-len     Minimum feature length in bp to include [default: 1]
    -x, --max-len     Maximum feature length in bp (0 = no limit) [default: 0]
    --dry-run         Print samtools commands without running them
    --no-samtools     Write a BED file instead of running samtools
"""

import argparse
import subprocess
import sys
import os
from dataclasses import dataclass
from typing import Optional


# Feature types and what they mean
FEATURE_DESCRIPTIONS = {
    "GAP":  "Novel sequence in assembly (absent from reference) — top insertion candidate",
    "DUP":  "Duplicated region in assembly — may contain inserted gene copies",
    "BRK":  "Alignment breakpoint — marks discontinuity, usually not novel sequence",
    "JMP":  "Alignment jump — inversion or rearrangement, usually not novel sequence",
    "INV":  "Inversion",
    "SEQ":  "Sequence-level difference",
}

CANDIDATE_TYPES = {"GAP", "DUP"}  # Default types worth extracting


@dataclass
class DiffEntry:
    contig: str
    feat_type: str
    start: int
    end: int
    length: int

    @property
    def is_valid(self) -> bool:
        """Filter out negative-length entries (overlaps/inversions, not insertions)."""
        return self.length > 0 and self.start > 0 and self.end > 0

    @property
    def bed_start(self) -> int:
        """BED format is 0-based."""
        return self.start - 1

    def __str__(self):
        return f"{self.contig}:{self.start}-{self.end} ({self.feat_type}, {self.length} bp)"


def parse_diff_file(filepath: str) -> list[DiffEntry]:
    """Parse a show-diff output file into a list of DiffEntry objects."""
    entries = []
    header_skipped = False

    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            # Skip the two header lines (path lines and NUCMER marker)
            if line.startswith("/") or line == "NUCMER":
                continue
            # Skip the column header line
            if line.startswith("[SEQ]"):
                header_skipped = True
                continue

            parts = line.split("\t")
            if len(parts) < 5:
                continue

            contig    = parts[0]
            feat_type = parts[1]
            try:
                s1  = int(parts[2])
                e1  = int(parts[3])
                len1 = int(parts[4])
            except ValueError:
                continue

            entries.append(DiffEntry(
                contig=contig,
                feat_type=feat_type,
                start=min(s1, e1),   # normalise so start <= end
                end=max(s1, e1),
                length=abs(len1),
            ))

    return entries


def filter_entries(
    entries: list[DiffEntry],
    types: set[str],
    min_len: int,
    max_len: int,
) -> list[DiffEntry]:
    """Filter entries by type, validity, and length."""
    out = []
    for e in entries:
        if e.feat_type not in types:
            continue
        if not e.is_valid:
            continue
        if e.length < min_len:
            continue
        if max_len > 0 and e.length > max_len:
            continue
        out.append(e)
    return out


def print_summary(all_entries: list[DiffEntry], filtered: list[DiffEntry], types: set[str]):
    """Print a human-readable summary of what was found."""
    print("\n=== show-diff summary ===")
    type_counts: dict[str, int] = {}
    for e in all_entries:
        type_counts[e.feat_type] = type_counts.get(e.feat_type, 0) + 1

    for t, n in sorted(type_counts.items()):
        flag = " <-- extracting" if t in types else ""
        desc = FEATURE_DESCRIPTIONS.get(t, "")
        print(f"  {t:6s}  {n:3d}  {desc}{flag}")

    print(f"\n{len(filtered)} entries selected for extraction:\n")
    for e in sorted(filtered, key=lambda x: (-x.length, x.contig)):
        print(f"  {e}")
    print()


def write_bed(entries: list[DiffEntry], bed_path: str):
    """Write a BED file for use with bedtools getfasta."""
    with open(bed_path, "w") as fh:
        for e in entries:
            name = f"{e.contig}_{e.feat_type}_{e.start}_{e.end}"
            fh.write(f"{e.contig}\t{e.bed_start}\t{e.end}\t{name}\n")
    print(f"BED file written: {bed_path}")
    print(f"  bedtools getfasta -fi <assembly> -bed {bed_path} -fo candidates.fna -name")


def run_samtools(
    entries: list[DiffEntry],
    assembly: str,
    output: str,
    dry_run: bool,
):
    """Run samtools faidx to extract all candidate regions."""
    regions = [f"{e.contig}:{e.start}-{e.end}" for e in entries]

    cmd = ["samtools", "faidx", assembly] + regions

    if dry_run:
        print("=== samtools command (dry run) ===")
        print(" \\\n  ".join(cmd))
        print(f"  > {output}\n")
        return

    # Check samtools is available
    try:
        subprocess.run(["samtools", "--version"], capture_output=True, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        print("ERROR: samtools not found. Use --no-samtools to write a BED file instead.",
              file=sys.stderr)
        sys.exit(1)

    print(f"Running samtools faidx on {len(regions)} regions...")
    with open(output, "w") as out_fh:
        result = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(f"ERROR: samtools failed:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)

    print(f"Sequences written to: {output}")
    print(f"\nNext steps:")
    print(f"  prokka --outdir prokka_out --prefix candidates {output}")
    print(f"  blastn -query {output} -db nt -remote -out blast_nt.txt -outfmt 6")


def main():
    parser = argparse.ArgumentParser(
        description="Parse MUMmer4 show-diff output and extract candidate insertion sequences.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("-d", "--diff",     required=True,  help="show-diff output file")
    parser.add_argument("-a", "--assembly", required=False, help="Assembly FASTA file")
    parser.add_argument("-o", "--output",   default="candidates.fna", help="Output FASTA [default: candidates.fna]")
    parser.add_argument("-t", "--types",    default="GAP,DUP",
                        help="Feature types to extract, comma-separated [default: GAP,DUP]")
    parser.add_argument("-m", "--min-len",  type=int, default=1,
                        help="Minimum feature length in bp [default: 1]")
    parser.add_argument("-x", "--max-len",  type=int, default=0,
                        help="Maximum feature length in bp, 0=no limit [default: 0]")
    parser.add_argument("--dry-run",        action="store_true",
                        help="Print samtools command without running it")
    parser.add_argument("--no-samtools",    action="store_true",
                        help="Write a BED file instead of running samtools")
    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.diff):
        print(f"ERROR: diff file not found: {args.diff}", file=sys.stderr)
        sys.exit(1)

    if not args.no_samtools and not args.dry_run and not args.assembly:
        print("ERROR: --assembly is required unless --no-samtools or --dry-run is set.",
              file=sys.stderr)
        sys.exit(1)

    if args.assembly and not os.path.exists(args.assembly):
        print(f"ERROR: assembly file not found: {args.assembly}", file=sys.stderr)
        sys.exit(1)

    types = {t.strip().upper() for t in args.types.split(",")}

    # Parse and filter
    all_entries = parse_diff_file(args.diff)
    filtered    = filter_entries(all_entries, types, args.min_len, args.max_len)

    print_summary(all_entries, filtered, types)

    if not filtered:
        print("No entries passed filters. Exiting.")
        sys.exit(0)

    # Output
    if args.no_samtools:
        bed_path = args.output.replace(".fna", ".bed").replace(".fa", ".bed")
        if not bed_path.endswith(".bed"):
            bed_path += ".bed"
        write_bed(filtered, bed_path)
    else:
        run_samtools(filtered, args.assembly, args.output, args.dry_run)


if __name__ == "__main__":
    main()
