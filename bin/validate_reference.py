#!/usr/bin/env python3
"""
Pre-flight checks on the host + virus reference, run before STAR
genomeGenerate. Each failure below gives a successful index and wrong counts.

  1. Host GTF and host FASTA disagree on contig names.
     Ensembl GTF ("1", "MT") against a UCSC/GENCODE FASTA ("chr1", "chrM"):
     STAR drops every annotation row whose contig is absent, without an
     error. The index builds, reads align, and ReadsPerGene / featureCounts
     assign almost nothing.

  2. Viral contig name collides with a host contig name.
     Two sequences called the same thing in one FASTA: the index is built,
     and which sequence a read was counted against is undefined.

  3. Viral FASTA empty or trivially short (wrong file, failed download).
     Infected samples would then show zero viral load and fail the viral
     gate later for the wrong reason; better to fail here.

  4. Host GTF lacks gene_name on most genes.
     DE still runs, but GSEA keys on symbols, so most of the ranking would
     be dropped. Warned here, enforced as a coverage check in fgsea.

Contig lists are passed in as files of FASTA header names (one per line),
extracted with `grep '^>'` by the calling process: reading a 3 GB FASTA in
Python for its headers is needlessly slow.
"""
import argparse
import gzip
import re
import sys
from collections import Counter

GENE_ID_RE = re.compile(r'gene_id\s+"([^"]+)"')
GENE_NAME_RE = re.compile(r'gene_name\s+"([^"]+)"')


def open_maybe_gz(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p, "rt")


def read_headers(path):
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                out.append(line.lstrip(">").split()[0])
    return out


def fasta_length(path):
    n = 0
    with open_maybe_gz(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                n += len(line.strip())
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host-contigs", required=True, help="file: host FASTA header names")
    ap.add_argument("--host-gtf", required=True)
    ap.add_argument("--viral-fasta", required=True)
    ap.add_argument("--min-exon-on-fasta", type=float, default=0.99)
    ap.add_argument("--min-viral-length", type=int, default=1000)
    ap.add_argument("--min-gene-name-frac", type=float, default=0.5)
    args = ap.parse_args()

    fail, warn = [], []
    def ok(m):  print(f"  [ ok ] {m}")
    def bad(m): print(f"  [FAIL] {m}"); fail.append(m)
    def wrn(m): print(f"  [warn] {m}"); warn.append(m)

    host = read_headers(args.host_contigs)
    host_set = set(host)
    print(f"[1] host GTF contigs vs host FASTA ({len(host)} FASTA contigs)")
    if not host:
        bad("host FASTA contains no contigs")
    exon_contigs = Counter()
    genes, named = set(), set()
    with open_maybe_gz(args.host_gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if len(f) < 9 or f[2] != "exon":
                continue
            exon_contigs[f[0]] += 1
            m = GENE_ID_RE.search(f[8])
            if m:
                genes.add(m.group(1))
                if GENE_NAME_RE.search(f[8]):
                    named.add(m.group(1))
    n_exon = sum(exon_contigs.values())
    if n_exon == 0:
        bad("host GTF has no exon rows: featureCounts -t exon would count nothing")
    else:
        on = sum(n for c, n in exon_contigs.items() if c in host_set)
        frac = on / n_exon
        absent = sorted((c for c in exon_contigs if c not in host_set),
                        key=lambda c: -exon_contigs[c])
        if frac < args.min_exon_on_fasta:
            hint = ""
            if absent and any(("chr" + c) in host_set for c in absent[:5]):
                hint = " (GTF uses '1', FASTA uses 'chr1': Ensembl annotation on a UCSC/GENCODE genome)"
            elif absent and any(c.startswith("chr") and c[3:] in host_set for c in absent[:5]):
                hint = " (GTF uses 'chr1', FASTA uses '1')"
            bad(f"only {frac:.1%} of exon rows are on contigs present in the FASTA; "
                f"STAR silently drops the rest{hint}. First absent: {absent[:5]}")
        else:
            ok(f"{frac:.2%} of {n_exon:,} exon rows on FASTA contigs")
            if absent:
                wrn(f"exon rows on {len(absent)} contig(s) absent from the FASTA "
                    f"(dropped by STAR): {absent[:5]}")

    print("[2] host GTF gene_name coverage")
    if genes:
        gfrac = len(named) / len(genes)
        msg = f"{len(named):,}/{len(genes):,} genes carry gene_name ({gfrac:.1%})"
        (ok if gfrac >= args.min_gene_name_frac else wrn)(msg)

    print("[3] viral FASTA")
    viral = []
    with open_maybe_gz(args.viral_fasta) as fh:
        viral = [l[1:].split()[0] for l in fh if l.startswith(">")]
    if not viral:
        bad("viral FASTA contains no contigs")
    else:
        vlen = fasta_length(args.viral_fasta)
        print(f"      contigs: {viral}  total length: {vlen:,} bp")
        if vlen < args.min_viral_length:
            bad(f"viral genome is {vlen} bp; expected a full genome (wrong or truncated file?)")
        else:
            ok("viral genome length plausible")
        dup = sorted(set(viral) & host_set)
        if dup:
            bad(f"viral contig name(s) also used by the host FASTA: {dup}")
        else:
            ok("no host/virus contig-name collision")
        if len(set(viral)) != len(viral):
            bad("duplicate contig names inside the viral FASTA")

    print("=" * 60)
    print(f"RESULT: {len(fail)} failure(s), {len(warn)} warning(s)")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
