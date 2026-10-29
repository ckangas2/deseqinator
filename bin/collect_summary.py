#!/usr/bin/env python3
"""
Two run-level tables:

  summary.tsv     one row per sample: reads in/out of fastp, STAR unique and
                  multi-mapping rates, featureCounts assignment, viral load,
                  strandedness call
  contrasts.tsv   one row per contrast: DE counts, label-swap baseline,
                  top gene sets, rank-shuffle control
"""
import argparse
import csv
import json
from pathlib import Path


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def star_log(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if "|" in line:
                k, v = line.split("|", 1)
                out[k.strip()] = v.strip()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastp", nargs="*", default=[])
    ap.add_argument("--star", nargs="*", default=[])
    ap.add_argument("--assignment", required=True)
    ap.add_argument("--viral", required=True)
    ap.add_argument("--strandedness", required=True)
    ap.add_argument("--deseq", nargs="*", default=[])
    ap.add_argument("--fgsea", nargs="*", default=[])
    ap.add_argument("--out-samples", default="summary.tsv")
    ap.add_argument("--out-contrasts", default="contrasts.tsv")
    args = ap.parse_args()

    rows = {}
    def row(s):
        return rows.setdefault(s, {"sample": s})

    for f in args.fastp:
        s = Path(f).name.removesuffix(".fastp.json")
        j = json.loads(Path(f).read_text())["summary"]
        row(s)["reads_raw"] = j["before_filtering"]["total_reads"]
        row(s)["reads_trimmed"] = j["after_filtering"]["total_reads"]
    for f in args.star:
        s = Path(f).name.removesuffix(".Log.final.out")
        log = star_log(f)
        row(s)["star_unique_pct"] = log.get("Uniquely mapped reads %", "")
        row(s)["star_multi_pct"] = log.get("% of reads mapped to multiple loci", "")
    for r in read_tsv(args.assignment):
        row(r["sample"])["assigned_frac"] = r["assigned_frac"]
    for r in read_tsv(args.viral):
        row(r["sample"])["condition"] = r["condition"]
        row(r["sample"])["viral_frac"] = r["viral_frac"]
    for r in read_tsv(args.strandedness):
        row(r["sample"])["fraction_forward"] = r["fraction_forward"]
        row(r["sample"])["strandedness"] = r["call"]

    cols = ["sample", "condition", "reads_raw", "reads_trimmed", "star_unique_pct",
            "star_multi_pct", "assigned_frac", "viral_frac", "fraction_forward", "strandedness"]
    with open(args.out_samples, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for s in sorted(rows):
            w.writerow({c: rows[s].get(c, "") for c in cols})

    contrasts = {}
    for f in args.deseq:
        for r in read_tsv(f):
            contrasts.setdefault(r["id"], {}).update(r)
    for f in args.fgsea:
        for r in read_tsv(f):
            contrasts.setdefault(r["id"], {}).update(r)
    ccols = ["id", "cell_line", "reference", "target", "design", "n_reference", "n_target",
             "genes_tested", "n_up", "n_down", "swap_max_de", "sets_sig_up", "sets_sig_down",
             "rank_shuffle_hits", "median_set_coverage", "top_up", "top_down"]
    with open(args.out_contrasts, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ccols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for k in sorted(contrasts):
            w.writerow({c: contrasts[k].get(c, "") for c in ccols})
    print(Path(args.out_samples).read_text())
    print(Path(args.out_contrasts).read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
