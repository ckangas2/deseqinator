#!/usr/bin/env python3
"""
Compare this run's counts against the authors' deposited count matrix and
catch sample swaps.

For each sample with a `deposited_id` in the samplesheet, correlate its
counts (Spearman, genes expressed in both) against EVERY deposited column.
The declared column must be the best match. If another deposited column
fits better, the samplesheet maps runs to samples wrongly, which is the
single most damaging error in a reanalysis: DESeq2 cannot detect it, and it
flips or erases the effect being measured.

Deposited matrices are keyed by symbol or by Ensembl ID; the key type is
detected, and our counts are aggregated to the same key via the GTF table.
Column names are compared after lowercasing and stripping punctuation, so
"Series7_Calu3_Mock_1" matches "series7-calu3-mock-1".
"""
import argparse
import csv
import re
import sys

import numpy as np
import pandas as pd


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def spearman(a, b):
    ra = pd.Series(a).rank().to_numpy()
    rb = pd.Series(b).rank().to_numpy()
    return float(np.corrcoef(ra, rb)[0, 1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ours", required=True, help="gene_id x sample counts")
    ap.add_argument("--annotation", required=True, help="gene_annotation.tsv")
    ap.add_argument("--deposited", required=True)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--min-mean", type=float, default=10.0)
    ap.add_argument("--min-rho", type=float, default=0.9)
    ap.add_argument("--out", default="reconciliation.tsv")
    ap.add_argument("--out-matrix", default="reconciliation_matrix.tsv")
    args = ap.parse_args()

    sep = "," if args.deposited.removesuffix(".gz").endswith(".csv") else "\t"
    dep = pd.read_csv(args.deposited, sep=sep, index_col=0)
    dep = dep.select_dtypes(include=[np.number])
    dep.index = dep.index.astype(str)

    ours = pd.read_csv(args.ours, sep="\t", index_col=0)
    ann = pd.read_csv(args.annotation, sep="\t", dtype=str).fillna("")

    ens_like = np.mean([i.startswith("ENS") for i in dep.index[:1000]]) > 0.5
    if ens_like:
        dep.index = [i.split(".")[0] for i in dep.index]
        key = ann.set_index("gene_id")["gene_id_unversioned"]
        keytype = "Ensembl gene ID"
    else:
        key = ann.set_index("gene_id")["gene_name"]
        keytype = "gene symbol"
    ours = ours.groupby(key.reindex(ours.index).fillna("").to_numpy()).sum()
    ours = ours.loc[ours.index != ""]
    dep = dep.groupby(level=0).sum()

    shared = ours.index.intersection(dep.index)
    print(f"deposited matrix keyed by {keytype}; {len(shared):,} shared genes")
    if len(shared) < 1000 and len(shared) < 0.5 * min(len(ours), len(dep)):
        print("FAIL: too few shared genes to compare; wrong deposited file or key type?",
              file=sys.stderr)
        return 1
    o, d = ours.loc[shared], dep.loc[shared]
    keep = (o.mean(axis=1) >= args.min_mean) & (d.mean(axis=1) >= args.min_mean)
    o, d = o.loc[keep], d.loc[keep]

    with open(args.samplesheet, newline="") as fh:
        ss = [{k.strip(): (v or "").strip() for k, v in r.items()} for r in csv.DictReader(fh)]
    declared = {r["sample"]: r.get("deposited_id", "") for r in ss}
    dep_norm = {norm(c): c for c in d.columns}

    rho = pd.DataFrame({s: [spearman(o[s], d[c]) for c in d.columns] for s in o.columns},
                       index=d.columns).T
    rho.to_csv(args.out_matrix, sep="\t", float_format="%.4f")

    problems = []
    with open(args.out, "w") as fh:
        fh.write("sample\tdeposited_id\trho_declared\tbest_match\trho_best\tverdict\n")
        for s in o.columns:
            want = declared.get(s, "")
            if not want:
                fh.write(f"{s}\t\t\t\t\tskipped\n")
                continue
            col = dep_norm.get(norm(want))
            best = rho.loc[s].idxmax()
            if col is None:
                problems.append(f"{s}: deposited_id '{want}' is not a column of the deposited matrix")
                fh.write(f"{s}\t{want}\t\t{best}\t{rho.loc[s, best]:.4f}\tFAIL\n")
                continue
            r_dec, r_best = rho.loc[s, col], rho.loc[s, best]
            verdict = "ok"
            if best != col:
                verdict = "FAIL"
                problems.append(f"{s}: declared as '{col}' (rho {r_dec:.3f}) but matches "
                                f"'{best}' better (rho {r_best:.3f}). Sample swap?")
            elif r_dec < args.min_rho:
                verdict = "low"
                print(f"warn: {s} vs {col}: rho {r_dec:.3f} < {args.min_rho}")
            fh.write(f"{s}\t{col}\t{r_dec:.4f}\t{best}\t{r_best:.4f}\t{verdict}\n")
            print(f"  {s:<30} -> {col:<30} rho {r_dec:.4f}  best {best}  {verdict}")

    if problems:
        for p in problems:
            print(f"FAIL: {p}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
