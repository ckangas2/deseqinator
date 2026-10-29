#!/usr/bin/env python3
"""
Viral read fraction per sample, and a gate on it.

Each input row is <sample> <primary alignments> <primary alignments on viral
contigs>, from samtools on the coordinate-sorted BAM. Counting by contig
rather than by viral gene sidesteps the viral annotation entirely:
SARS-CoV-2 subgenomic mRNAs are nested and share their 3' ends, so a
gene-level count assigns many viral reads to more than one gene and
featureCounts discards them as ambiguous. Load per contig has no such loss.

The gate is the cheapest sample-swap and failed-infection check there is:
  - a mock sample carrying viral reads is contaminated or mislabelled
  - an infected sample with ~no viral reads did not get infected, or is
    mislabelled
Either way the DE result for that contrast would be wrong in a way DESeq2
cannot see, so the run stops.
"""
import argparse
import csv
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("loads", nargs="+")
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--mock-conditions", default="Mock",
                    help="comma-separated condition values that must carry no virus")
    ap.add_argument("--max-mock", type=float, default=0.001)
    ap.add_argument("--min-infected", type=float, default=0.001)
    ap.add_argument("--out", default="viral_load.tsv")
    args = ap.parse_args()

    mock = {c.strip() for c in args.mock_conditions.split(",") if c.strip()}
    with open(args.samplesheet, newline="") as fh:
        cond = {r["sample"].strip(): r["condition"].strip() for r in csv.DictReader(fh)}

    rows = []
    for path in args.loads:
        with open(path) as fh:
            for line in fh:
                f = line.split()
                if len(f) != 3:
                    continue
                total, viral = int(f[1]), int(f[2])
                rows.append((f[0], total, viral, viral / total if total else 0.0))

    problems = []
    missing = sorted(set(cond) - {r[0] for r in rows})
    if missing:
        problems.append(f"no viral load reported for {missing}")
    with open(args.out, "w") as fh:
        fh.write("sample\tcondition\tprimary_alignments\tviral_alignments\tviral_frac\tgate\n")
        for s, total, viral, frac in sorted(rows):
            c = cond.get(s)
            if c is None:
                problems.append(f"{s} is not in the samplesheet")
                continue
            if c in mock:
                verdict = "ok" if frac <= args.max_mock else "FAIL"
                if verdict == "FAIL":
                    problems.append(f"{s} is labelled {c} but {frac:.3%} of its alignments are "
                                    f"viral (max {args.max_mock:.3%}): contamination or a "
                                    f"sample swap")
            else:
                verdict = "ok" if frac >= args.min_infected else "FAIL"
                if verdict == "FAIL":
                    problems.append(f"{s} is labelled {c} but only {frac:.4%} of alignments are "
                                    f"viral (min {args.min_infected:.3%}): infection failed, "
                                    f"or the sample is mislabelled")
            fh.write(f"{s}\t{c}\t{total}\t{viral}\t{frac:.6f}\t{verdict}\n")
            print(f"  {s:<30} {c:<16} viral {frac:8.4%}  {verdict}")

    if problems:
        for p in problems:
            print(f"FAIL: {p}", file=sys.stderr)
        print("Thresholds: --max_mock_viral_frac / --min_infected_viral_frac. "
              "Conditions treated as uninfected: --mock_conditions.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
