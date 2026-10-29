# deseqinator

[![CI](https://github.com/ckangas2/deseqinator/actions/workflows/ci.yml/badge.svg)](https://github.com/ckangas2/deseqinator/actions/workflows/ci.yml)

A Nextflow pipeline for bulk RNA-seq of virus infection: fastp → STAR → featureCounts → DESeq2 → fgsea against MSigDB Hallmarks. Reads align to a concatenated host + virus genome, so each library's viral load is measured from the same alignment that produces host counts.

The target dataset is GSE147507 (Blanco-Melo et al., *Cell* 2020): SARS-CoV-2 vs mock in Calu-3 and A549-ACE2, biological triplicates. It has been reanalysed many times, so there are published DE counts to land near, and the authors' raw count matrix to reconcile against.

Each stage is the standard tool with standard settings. What this repo adds is the checks between the stages. Every failure in the table below produces a complete results directory with plausible numbers in it; here, each one stops the run.

## Try it

You need [Nextflow](https://www.nextflow.io/docs/latest/install.html) ≥ 24.04 and Docker.

```
nextflow run ckangas2/deseqinator -profile test,docker
```

This runs nine synthetic single-end libraries (200 genes, 40 planted as 4-fold DE, infected libraries ~20% viral) through the whole pipeline. `tests/check_test_results.py` then asserts the planted truth comes back: DE genes recovered with the right sign, the up/down gene sets enriched in the planted direction, the unchanged set not, decoy samples kept out of the contrast. CI also runs seven negative cases, one per gate, and requires each to fail.

No Docker? `-profile test,singularity`, `test,apptainer` or `test,conda`.

## Run it on GSE147507

```
git clone https://github.com/ckangas2/deseqinator && cd deseqinator
bash scripts/fetch_reference.sh      # GRCh38 + GENCODE v44, SARS-CoV-2 USA-WA1/2020, deposited counts
bash scripts/fetch_runs.sh --list    # how every GEO sample title parses
bash scripts/fetch_runs.sh           # Calu-3 + A549-ACE2, 12 runs, writes samplesheet.csv
nextflow run . -profile gse147507,docker
```

Run accessions are not hand-copied. `fetch_runs.sh` pulls ENA's run report for the project and parses each sample title. It refuses the selection unless each requested series is exactly 3 Mock + 3 SARS-CoV-2 of one cell line.

## Run it on your data

**Samplesheet** (CSV). `fastq_2` is empty for single-end. `batch` and `deposited_id` are optional.

```
sample,fastq_1,fastq_2,cell_line,condition,replicate,batch,deposited_id
Calu3_mock_1,fastq/SRR1.fastq.gz,,Calu3,Mock,1,Series7,Series7_Calu3_Mock_1
```

**Contrasts** (CSV). One row per comparison, scoped to one cell line. `*` pools all lines and adds `cell_line` to the design.

```
id,cell_line,reference,target
calu3_sars2_vs_mock,Calu3,Mock,SARS-CoV-2
```

**Run**

```
nextflow run ckangas2/deseqinator -profile docker \
    --input samplesheet.csv --contrasts contrasts.csv \
    --host_fasta genome.fa.gz --host_gtf genes.gtf.gz --viral_fasta virus.fa
```

`nextflow run ckangas2/deseqinator --help` lists every parameter.

## What gets caught

| failure | what you see | what it actually is | caught by |
|---|---|---|---|
| cell lines pooled in `~ condition` | a DE list | cell-line differences inflate residual variance, or leak into the condition effect if unbalanced | contrasts are per cell line; `*` adds `cell_line` to the design (`CHECK_DESIGN`) |
| a third arm (e.g. Ruxolitinib) in the contrast | a DE list | treated samples counted as "infected" | a contrast selects exactly two conditions; exclusions are printed |
| batch confounded with condition | an error deep in DESeq2, or with partial confounding the effect quietly assigned to batch | the design cannot separate them | model-matrix rank check, before alignment |
| n = 1 or 2 per arm | confident p-values | dispersion from the trend alone | `--min_replicates`, before alignment |
| wrong strandedness in featureCounts | a full matrix | most reads discarded as `Unassigned_NoFeatures`; DE built on a few percent of the data | strand called from STAR GeneCounts; contradicting `--strandedness` fails (`INFER_STRANDEDNESS`) |
| GTF from another build, or `1` vs `chr1` | a full matrix of small numbers | STAR drops annotation on absent contigs without warning | contig check before indexing (`VALIDATE_REFERENCE`); assignment-rate floor (`FEATURECOUNTS`) |
| mock contaminated, or samples swapped | a DE list | the comparison is not the one you think | viral load per library: mock must be ~0, infected must not be (`VIRAL_GATE`) |
| runs mapped to the wrong GEO samples | a DE list | a swap DESeq2 cannot detect | correlation against the authors' deposited counts; declared column must be the best match (`RECONCILE`) |
| sets keyed differently from the ranks (Ensembl vs symbol, mouse vs human) | "no significant gene sets" | every set fell below `minSize` | median set coverage floor (`FGSEA`) |
| enrichment not driven by the ranking | significant gene sets | e.g. a broken ranking that still sorts | rank-shuffle control must give ~0 hits (`FGSEA`) |
| DE indistinguishable from noise | a DE list | n = 3 is small | every balanced relabelling refitted; real DE count reported against the label-swap maximum (`DESEQ2`) |

## Design decisions worth stating

**Viral load by contig, not by viral gene.** SARS-CoV-2 subgenomic mRNAs are nested and share 3' ends. A gene-level count makes most viral reads multi-gene, and featureCounts discards them as ambiguous. Counting primary alignments on the viral contig has no such loss and needs no viral GTF.

**Host genes only in the DE matrix.** Infected libraries can be a third viral. Median-of-ratios on host genes is unaffected by that; any normalisation that sees total reads is not. The viral contig carries no annotation, so viral reads never enter `counts.tsv`.

**Condition is recoded before fitting.** `SARS-CoV-2` becomes `target` inside the model, so the coefficient is always `condition_target_vs_reference`. Raw labels pass through `make.names()` in `resultsNames()`, and a coefficient looked up by a mangled name is an easy silent miss.

**GSEA ranks on the Wald statistic.** Raw log2FC puts noisy low-count genes at both ends of the list. Shrunken LFC (apeglm) is reported for plots and effect sizes.

**Symbols come from the GTF that was counted against.** Not from `org.Hs.eg.db` at analysis time, which resolves against whatever release is installed. When several Ensembl IDs share a symbol, the one with the largest |stat| is kept, and the number dropped is reported.

**Gene sets are pinned.** `msigdbr` 7.5.1 ships MSigDB 7.5.1 inside the package, so `msigdb:H` needs no network and gives the same 50 sets every time. Pass `--gene_sets file.gmt` for anything else.

**The design is decided once.** `check_design.py` writes `contrasts.validated.tsv` (samples, arms, formula) before any alignment. DESeq2 reads it and re-derives nothing, so the samples analysed are exactly the ones the design report listed.

## Prior art

[nf-core/rnaseq](https://nf-co.re/rnaseq) does fastp/STAR/featureCounts with far more options, and [nf-core/differentialabundance](https://nf-co.re/differentialabundance) does DESeq2 + GSEA. For production work, use those. This repo is narrower: one infection design, end to end, with the checks above wired in as hard stops, and a validation dataset with known answers.

## Outputs

```
results/
  summary.tsv                  per sample: reads, mapping, assignment, viral load, strand
  contrasts.tsv                per contrast: DE counts, label-swap baseline, top gene sets
  design/                      contrasts.validated.tsv, design_report.txt
  qc/                          strandedness.tsv, viral_load.tsv, reconciliation.tsv
  counts/                      counts.tsv (gene_id x sample), assignment.tsv, featureCounts raw
  deseq2/<id>/                 results, normalised counts, size factors, label swaps,
                               PCA / MA / dispersion / sample-distance plots
  fgsea/<id>/                  results, rank-shuffle control, set coverage, plots
  reference/                   validation log, gene_annotation.tsv, index check
  multiqc/  fastp/  star/  pipeline_info/
```

## Layout

```
main.nf                    design -> reference -> per sample -> gates -> per contrast
nextflow.config            parameters, resources, profiles
conf/                      test.config, gse147507.config
modules/local/             one process per file, each with a stub block
bin/
  check_design.py          samplesheet + contrasts -> validated design
  validate_reference.py    contig agreement, viral/host name clash
  infer_strandedness.py    strand call from STAR GeneCounts
  tidy_counts.py           featureCounts -> matrix; assignment floor
  viral_gate.py            viral load per library vs its label
  gtf_gene_table.py        gene_id -> symbol from the counted GTF
  run_deseq2.R             one contrast, QC plots, label swaps
  run_fgsea.R              preranked GSEA, coverage check, rank-shuffle control
  reconcile_counts.py      ours vs deposited counts; swap detection
  collect_summary.py       run-level tables
scripts/                   reference and run download, run selection
tests/
  make_test_data.py        regenerates tests/data/ byte-for-byte (CI checks this)
  test_checks.py           unit tests: one planted failure per check
  check_test_results.py    asserts the test profile against the planted truth
```

## Development

```
bash tests/run_tests.sh                          # unit tests, seconds, no STAR/R
nextflow run . -profile test -stub-run           # wiring only, no tools needed
nextflow run . -profile test,docker              # end to end
python3 tests/check_test_results.py results
```

## Status

- [x] Unit tests (28) pass; `nextflow lint` clean; stub run completes, and design failures stop it before alignment
- [x] Container build and conda solve
- [x] `-profile test,docker` end to end locally (Nextflow 26.10.4): 17/17 result checks pass
- [ ] The same in CI, plus the seven negative cases
- [ ] GSE147507 run; DE counts vs published reanalyses; Hallmark results; reconciliation against the deposited matrix → RESULTS.md
