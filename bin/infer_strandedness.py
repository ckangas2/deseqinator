#!/usr/bin/env python3
"""
Call library strandedness from STAR's ReadsPerGene.out.tab, for every sample,
and refuse to continue unless the call is clear and the samples agree.

Wrong strandedness in featureCounts is silent. Count a reverse-stranded
library with -s 1 and roughly 90% of reads land on the "wrong" strand and
become Unassigned_NoFeatures; the matrix is still full of numbers, DESeq2
still runs, and the DE list is built on a few percent of the data, weighted
toward genes with antisense transcription. Count a stranded library with
-s 0 and reads from overlapping antisense genes become ambiguous.

STAR's GeneCounts columns give the answer for free, from the alignment that
has to happen anyway:
    col 2  unstranded
    col 3  read 1 on the transcript strand      (htseq -s yes  / featureCounts -s 1)
    col 4  read 1 on the opposite strand        (htseq -s reverse / featureCounts -s 2)

fraction_forward = col3 / (col3 + col4)
    >= threshold              forward
    <= 1 - threshold          reverse
    within 0.5 +/- band       unstranded
    anything else             ambiguous -> stop

Writes strandedness.tsv (one row per sample) and strandedness.txt (the call).
"""
import argparse
import sys
from pathlib import Path

FC_CODE = {"unstranded": 0, "forward": 1, "reverse": 2}


def sample_name(path):
    name = Path(path).name
    for suffix in (".ReadsPerGene.out.tab", "ReadsPerGene.out.tab"):
        if name.endswith(suffix):
            return name[: -len(suffix)].rstrip(".")
    return name


def tally(path):
    s3 = s4 = s2 = 0
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 or f[0].startswith("N_"):
                continue
            s2 += int(f[1]); s3 += int(f[2]); s4 += int(f[3])
    return s2, s3, s4


def call(frac, n, threshold, band, min_reads):
    if n < min_reads:
        return "ambiguous"
    if frac >= threshold:
        return "forward"
    if frac <= 1 - threshold:
        return "reverse"
    if abs(frac - 0.5) <= band:
        return "unstranded"
    return "ambiguous"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+", help="<sample>.ReadsPerGene.out.tab")
    ap.add_argument("--declared", default="auto",
                    choices=["auto", "unstranded", "forward", "reverse"])
    ap.add_argument("--threshold", type=float, default=0.8)
    ap.add_argument("--band", type=float, default=0.1)
    ap.add_argument("--min-reads", type=int, default=1000,
                    help="gene-assigned reads needed for a call")
    ap.add_argument("--out-tsv", default="strandedness.tsv")
    ap.add_argument("--out-call", default="strandedness.txt")
    args = ap.parse_args()

    rows = []
    for f in sorted(args.files):
        s2, s3, s4 = tally(f)
        n = s3 + s4
        frac = s3 / n if n else float("nan")
        c = call(frac, n, args.threshold, args.band, args.min_reads)
        rows.append((sample_name(f), s2, s3, s4, frac, c))

    with open(args.out_tsv, "w") as fh:
        fh.write("sample\tunstranded\tforward\treverse\tfraction_forward\tcall\n")
        for r in rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]:.4f}\t{r[5]}\n")
    for r in rows:
        print(f"  {r[0]:<30} fraction_forward={r[4]:.3f}  -> {r[5]}")

    problems = []
    amb = [r[0] for r in rows if r[5] == "ambiguous"]
    if amb:
        problems.append(f"no clear strandedness call for {amb}. Between thresholds, or too "
                        f"few gene-assigned reads (< {args.min_reads}). Check the reference "
                        f"annotation and the library, then set --strandedness explicitly.")
    calls = {r[5] for r in rows if r[5] != "ambiguous"}
    if len(calls) > 1:
        by = {c: [r[0] for r in rows if r[5] == c] for c in calls}
        problems.append(f"libraries disagree on strandedness: {by}. They were probably "
                        f"prepared with different kits; count them in separate runs.")
    inferred = calls.pop() if len(calls) == 1 and not amb else None
    if inferred and args.declared != "auto" and args.declared != inferred:
        problems.append(f"--strandedness {args.declared} was given but every library reads as "
                        f"{inferred}. Counting with the declared setting would discard most "
                        f"reads. Remove the flag or correct it.")

    if problems:
        for p in problems:
            print(f"FAIL: {p}", file=sys.stderr)
        return 1
    Path(args.out_call).write_text(inferred + "\n")
    print(f"strandedness: {inferred} (featureCounts -s {FC_CODE[inferred]})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
