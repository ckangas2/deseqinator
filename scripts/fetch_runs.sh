#!/usr/bin/env bash
# Select and download the GSE147507 runs used by -profile gse147507.
#
#   bash scripts/fetch_runs.sh                 # Calu-3 (Series7) + A549-ACE2 (Series6)
#   bash scripts/fetch_runs.sh --list          # show how every sample title parses
#   SERIES=Series7 bash scripts/fetch_runs.sh  # Calu-3 only
#
# Run accessions and URLs come from ENA's run report at run time, not from
# a hand-copied table; select_runs.py refuses the selection unless each
# series is exactly 3 Mock + 3 SARS-CoV-2 of one cell line.
# HTTPS with aria2c (16 connections) where available: ENA's anonymous FTP
# is capped per stream. Every file is checked against ENA's MD5.
set -euo pipefail
cd "$(dirname "$0")/.."

PROJECT="${PROJECT:-PRJNA615032}"   # GSE147507
SERIES="${SERIES:-Series6,Series7}"
FASTQ_DIR="${FASTQ_DIR:-fastq}"
CONN="${CONN:-16}"
mkdir -p config "$FASTQ_DIR"

report=config/ena_runs.tsv
curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${PROJECT}&result=read_run&fields=run_accession,sample_alias,sample_title,library_layout,read_count,fastq_ftp,fastq_md5,fastq_bytes&format=tsv&limit=0" > "$report"
echo "ENA: $(($(wc -l < "$report") - 1)) runs in $PROJECT"

if [ "${1:-}" = "--list" ]; then
  python3 scripts/select_runs.py "$report" --series "$SERIES" --list | column -t -s $'\t'
  exit 0
fi
python3 scripts/select_runs.py "$report" --series "$SERIES" --fastq-dir "$FASTQ_DIR"

while IFS=$'\t' read -r u md5 f; do
  if [ -f "$f" ] && [ "$(md5sum "$f" | cut -d' ' -f1)" = "$md5" ]; then
    echo "  ok (verified): $(basename "$f")"; continue; fi
  echo "  -> $(basename "$f")"
  if command -v aria2c >/dev/null 2>&1; then
    aria2c -x"$CONN" -s"$CONN" -c --summary-interval=30 -d "$(dirname "$f")" -o "$(basename "$f")" "$u"
  else
    curl -fL --retry 3 -C - -o "$f" "$u"
  fi
  [ "$(md5sum "$f" | cut -d' ' -f1)" = "$md5" ] || { echo "MD5 MISMATCH: $f" >&2; exit 1; }
done < config/downloads.tsv
echo "samplesheet.csv written; FASTQs in $FASTQ_DIR/"
