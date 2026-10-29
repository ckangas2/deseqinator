#!/usr/bin/env python3
"""
Turn featureCounts output into a gene x sample matrix named by sample, and
stop if any library assigned too few reads to genes.

A low assignment rate is the common footprint of several silent failures:
wrong strandedness, a GTF from a different genome build, or contig-name
mismatch. featureCounts reports it in .summary, which nobody reads; the
matrix it writes looks normal. Here it is a hard threshold.

Column names in featureCounts output are the BAM paths as given. The BAMs
are named <sample>.bam by the alignment step; anything else is an error
rather than a guess.
"""
import argparse
import sys
from pathlib import Path


def sample_from_bam(col):
    name = Path(col).name
    if not name.endswith(".bam"):
        raise SystemExit(f"unexpected featureCounts column '{col}': expected <sample>.bam")
    return name[:-4]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--featurecounts", required=True)
    ap.add_argument("--summary", default="", help="defaults to <featurecounts>.summary")
    ap.add_argument("--min-assigned", type=float, default=0.3)
    ap.add_argument("--out-matrix", default="counts.tsv")
    ap.add_argument("--out-qc", default="assignment.tsv")
    args = ap.parse_args()

    summary = args.summary or args.featurecounts + ".summary"

    with open(args.featurecounts) as fh, open(args.out_matrix, "w") as out:
        header = None
        n = 0
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if header is None:
                if f[0] != "Geneid":
                    raise SystemExit(f"{args.featurecounts}: header does not start with Geneid")
                samples = [sample_from_bam(c) for c in f[6:]]
                if len(set(samples)) != len(samples):
                    raise SystemExit(f"duplicate sample columns: {samples}")
                header = samples
                out.write("gene_id\t" + "\t".join(samples) + "\n")
                continue
            out.write(f[0] + "\t" + "\t".join(f[6:]) + "\n")
            n += 1
    print(f"matrix: {n:,} genes x {len(header)} samples -> {args.out_matrix}")

    status = {}
    with open(summary) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        names = [sample_from_bam(c) for c in cols[1:]]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            status[f[0]] = [int(x) for x in f[1:]]

    totals = [sum(v[i] for v in status.values()) for i in range(len(names))]
    assigned = status.get("Assigned", [0] * len(names))
    reasons = [k for k in status if k != "Assigned"]

    problems = []
    with open(args.out_qc, "w") as fh:
        fh.write("sample\ttotal\tassigned\tassigned_frac\ttop_unassigned\n")
        for i, s in enumerate(names):
            frac = assigned[i] / totals[i] if totals[i] else 0.0
            top = max(reasons, key=lambda k: status[k][i]) if reasons else ""
            fh.write(f"{s}\t{totals[i]}\t{assigned[i]}\t{frac:.4f}\t{top}\n")
            print(f"  {s:<30} assigned {frac:6.1%}   largest loss: {top}")
            if frac < args.min_assigned:
                problems.append(f"{s}: {frac:.1%} of alignments assigned to genes "
                                f"(< {args.min_assigned:.0%}); largest category {top}")

    if problems:
        for p in problems:
            print(f"FAIL: {p}", file=sys.stderr)
        print("Usual causes: wrong strandedness, annotation from a different genome build, "
              "or GTF/FASTA contig-name mismatch.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
