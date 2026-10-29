# Results: GSE147507

Not run yet. This file is the template the GSE147507 run fills in, with the comparison points chosen in advance so the result is judged against them, not explained after.

## What to compare against

**DE counts.** One published reanalysis (DESeq2, p < 0.05, |log2FC| ≥ 1) reports:

| contrast | up | down |
|---|---|---|
| Calu-3 SARS-CoV-2 vs mock | 1,264 | 660 |
| A549-ACE2 SARS-CoV-2 vs mock | 1,083 | 352 |

That reanalysis thresholded on p, not padj. Report both, and expect the padj counts to be lower.

**Biology.** Blanco-Melo et al. describe a strong chemokine/inflammatory response with comparatively muted type I/III interferon. In Calu-3 the Hallmark expectations are `INTERFERON_ALPHA_RESPONSE`, `INTERFERON_GAMMA_RESPONSE` and `TNFA_SIGNALING_VIA_NFKB` enriched up. If none of them is, check upstream before interpreting.

**Reconciliation.** Every sample should match its own column of `GSE147507_RawReadCounts_Human` best. Spearman ρ will be below 1, since the authors used a different aligner and annotation.

## Sections to fill

1. Run QC: `summary.tsv` (mapping, assignment, strand call, viral load per library)
2. Design report: samples per contrast, what was excluded
3. DE: counts at padj < 0.05 and at p < 0.05, label-swap baseline, PCA
4. Hallmarks: top sets up/down, rank-shuffle control, set coverage
5. Reconciliation table
6. Anything that disagreed with expectations, and what it turned out to be
