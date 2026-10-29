#!/usr/bin/env bash
# Host and viral reference for -profile gse147507, plus the deposited counts.
#
#   GRCh38 primary assembly + GENCODE v44 primary-assembly annotation:
#     same provider for FASTA and GTF, so contig names agree ("chr1" in
#     both); GENCODE gene_name is the HGNC symbol MSigDB keys on.
#   SARS-CoV-2 USA-WA1/2020 (MN985325.1): the isolate Blanco-Melo et al.
#     infected with. NC_045512.2 (Wuhan-Hu-1) also works; set VIRUS_ACC.
#   GSE147507 raw counts: the authors' matrix, for sample-swap detection.
#
# Each download is checked by content, not exit status: a redirect page or
# a truncated file downloads "successfully".
set -euo pipefail
cd "$(dirname "$0")/.."
REF=reference
GENCODE=44
VIRUS_ACC="${VIRUS_ACC:-MN985325.1}"
mkdir -p "$REF"

get() { [ -s "$2" ] && { echo "  have $(basename "$2")"; return; }
        echo "  -> $(basename "$2")"
        if command -v aria2c >/dev/null 2>&1; then
          aria2c -x16 -s16 -c -d "$(dirname "$2")" -o "$(basename "$2")" "$1"
        else curl -fL --retry 3 -o "$2" "$1"; fi; }

base="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE}"
get "$base/GRCh38.primary_assembly.genome.fa.gz"                      "$REF/GRCh38.primary_assembly.genome.fa.gz"
get "$base/gencode.v${GENCODE}.primary_assembly.annotation.gtf.gz"    "$REF/gencode.v${GENCODE}.primary_assembly.annotation.gtf.gz"
n=$(zcat "$REF/GRCh38.primary_assembly.genome.fa.gz" | grep -c '^>')
[ "$n" -ge 25 ] || { echo "FATAL: GRCh38 FASTA has $n contigs" >&2; exit 1; }
echo "  GRCh38 contigs: $n"

vf="$REF/SARS-CoV-2_USA-WA1-2020.fa"
if [ ! -s "$vf" ]; then
  curl -fsS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${VIRUS_ACC}&rettype=fasta&retmode=text" > "$vf"
fi
len=$(grep -v '^>' "$vf" | tr -d '\n' | wc -c)
[ "$len" -gt 29000 ] || { echo "FATAL: $vf is $len bp; expected ~29.9 kb" >&2; exit 1; }
echo "  $(head -1 "$vf" | cut -c1-70)  ($len bp)"

# The supplementary file is located by listing the GEO directory, not by
# assuming its exact name.
geo="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147507/suppl/"
dep=$(curl -fsS "$geo" | grep -o 'GSE147507_RawReadCounts_Human[^"<]*' | head -1 || true)
if [ -n "$dep" ]; then
  get "$geo$dep" "$REF/GSE147507_RawReadCounts_Human.tsv.gz"
else
  echo "  WARN: no RawReadCounts_Human file found at $geo; run without --deposited_counts" >&2
fi
