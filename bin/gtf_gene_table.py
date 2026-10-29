#!/usr/bin/env python3
"""
gene_id -> gene_name / gene_type table, taken from the same GTF featureCounts
counted against.

Symbols come from the annotation that produced the counts, not from an
org.Hs.eg.db or biomaRt lookup at analysis time: those resolve against
whatever release is installed, so the mapping drifts silently between
machines and years. One GTF, one mapping.
"""
import argparse
import gzip
import re
import sys

ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out", default="gene_annotation.tsv")
    args = ap.parse_args()

    op = gzip.open if args.gtf.endswith(".gz") else open
    genes = {}
    with op(args.gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if len(f) < 9 or f[2] not in ("gene", "exon"):
                continue
            a = dict(ATTR_RE.findall(f[8]))
            gid = a.get("gene_id")
            if not gid or (gid in genes and f[2] == "exon"):
                continue
            genes[gid] = (a.get("gene_name", ""),
                          a.get("gene_type", a.get("gene_biotype", "")), f[0])

    with open(args.out, "w") as fh:
        fh.write("gene_id\tgene_id_unversioned\tgene_name\tgene_type\tchrom\n")
        for gid, (name, gtype, chrom) in genes.items():
            fh.write(f"{gid}\t{gid.split('.')[0]}\t{name}\t{gtype}\t{chrom}\n")
    named = sum(1 for v in genes.values() if v[0])
    print(f"{len(genes):,} genes, {named:,} with gene_name -> {args.out}")
    return 0 if genes else 1


if __name__ == "__main__":
    sys.exit(main())
