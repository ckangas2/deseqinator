# deseqinator runtime: fastp, STAR, samtools, featureCounts, MultiQC,
# DESeq2 / apeglm / fgsea / msigdbr, Python.
#
#   docker build -t ghcr.io/ckangas2/deseqinator:0.1.0 .
#
# Versions are pinned in environment.yml, shared with `-profile conda`.
# The pipeline's own scripts are NOT copied in: Nextflow puts bin/ on PATH
# at run time, and a baked-in copy would sit earlier on PATH and silently
# run stale code after every edit.
FROM mambaorg/micromamba:1.5.10-bookworm-slim

LABEL org.opencontainers.image.source="https://github.com/ckangas2/deseqinator" \
      org.opencontainers.image.description="Bulk RNA-seq DE + GSEA for virus infection"

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba clean --all --yes \
 && rm /tmp/environment.yml

ENV PATH="/opt/conda/bin:${PATH}"

# Fail the build, not the first pipeline run, if anything is missing.
RUN STAR --version | grep -qx 2.7.11b \
 && fastp --version 2>&1 | grep -q 0.23.4 \
 && featureCounts -v 2>&1 | grep -q 2.0.6 \
 && samtools --version | head -1 | grep -q 1.21 \
 && multiqc --version \
 && python -c "import numpy, pandas" \
 && Rscript -e 'suppressPackageStartupMessages({library(DESeq2); library(apeglm); library(fgsea); library(msigdbr); library(data.table); library(ggplot2)}); stopifnot(nrow(msigdbr(species="Homo sapiens", category="H")) > 7000)' \
 && command -v ps zcat >/dev/null
