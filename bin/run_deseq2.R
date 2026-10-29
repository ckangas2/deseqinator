#!/usr/bin/env Rscript
# One DESeq2 contrast, defined entirely by contrasts.validated.tsv.
#
# The design, the sample list and the two arms come from check_design.py.
# Nothing is re-derived here, so the samples analysed are exactly the ones
# the design report printed before alignment started.
#
# Choices worth stating:
#   - Condition is recoded to reference/target before fitting. Raw labels
#     such as "SARS-CoV-2" are mangled by make.names() inside resultsNames(),
#     and a coefficient looked up by a mangled name is an easy silent miss.
#     The coefficient is always condition_target_vs_reference.
#   - Counts are host genes only (viral reads sit on the viral contig, which
#     the host GTF does not annotate). Infected libraries can be >30% viral;
#     median-of-ratios on host genes is unaffected by that, total-count
#     scaling would not be.
#   - DEG counts use the unshrunken LFC and padj; the apeglm-shrunken LFC is
#     reported for ranking and plotting. GSEA ranks on the Wald statistic.
#   - Label swaps: the same model refitted on every balanced relabelling of
#     the samples. With n = 3 per arm there are 9 distinct ones, each mixing
#     the true arms. A real effect should give far more DE genes than any of
#     them; if not, the result is not distinguishable from relabelling.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
opt <- list()
i <- 1
while (i <= length(args)) {
  k <- sub("^--", "", args[i]); opt[[k]] <- args[i + 1]; i <- i + 2
}
need <- c("counts", "samplesheet", "contrasts", "id", "annotation",
          "fdr", "lfc", "min_count", "label_swaps", "seed")
miss <- setdiff(need, names(opt))
if (length(miss)) stop("missing arguments: ", paste(miss, collapse = ", "))
fdr       <- as.numeric(opt$fdr)
lfc_thr   <- as.numeric(opt$lfc)
min_count <- as.integer(opt$min_count)
do_swaps  <- tolower(opt$label_swaps) %in% c("true", "1", "yes")
set.seed(as.integer(opt$seed))
id <- opt$id

ct  <- read.delim(opt$contrasts, stringsAsFactors = FALSE, colClasses = "character")
row <- ct[ct$id == id, ]
if (nrow(row) != 1) stop("contrast '", id, "' not found exactly once in ", opt$contrasts)

counts <- read.delim(opt$counts, row.names = 1, check.names = FALSE)
ss  <- read.csv(opt$samplesheet, stringsAsFactors = FALSE, check.names = FALSE,
                strip.white = TRUE, colClasses = "character")
ann <- read.delim(opt$annotation, stringsAsFactors = FALSE, colClasses = "character")

samples <- strsplit(row$samples, ",")[[1]]
absent  <- setdiff(samples, colnames(counts))
if (length(absent)) stop("samples in the design but not in the count matrix: ",
                         paste(absent, collapse = ", "))

cd <- ss[match(samples, ss$sample), , drop = FALSE]
rownames(cd) <- cd$sample
cd$condition <- factor(ifelse(cd$condition == row$reference, "reference", "target"),
                       levels = c("reference", "target"))
covars <- if (nzchar(row$covariates)) strsplit(row$covariates, ",")[[1]] else character(0)
for (v in covars) cd[[v]] <- factor(cd[[v]])
design <- as.formula(paste("~", paste(c(covars, "condition"), collapse = " + ")))
stopifnot(identical(gsub(" ", "", deparse(design)), gsub(" ", "", row$design)))

m <- as.matrix(counts[, samples, drop = FALSE])
storage.mode(m) <- "integer"

# Keep genes with >= min_count reads in at least as many samples as the
# smaller arm: a gene expressed in only one arm is kept, which is the point.
smallest <- min(table(cd$condition))
keep <- rowSums(m >= min_count) >= smallest
message(sprintf("[%s] %s ; %d samples ; %d / %d genes pass the expression filter",
                id, deparse(design), ncol(m), sum(keep), nrow(m)))

fit <- function(mat, coldata) {
  dds <- DESeqDataSetFromMatrix(mat, coldata, design)
  DESeq(dds, quiet = TRUE)
}
n_de <- function(res) {
  sig <- !is.na(res$padj) & res$padj < fdr
  c(up = sum(sig & res$log2FoldChange >= lfc_thr),
    down = sum(sig & res$log2FoldChange <= -lfc_thr))
}

dds  <- fit(m[keep, , drop = FALSE], cd)
coef <- "condition_target_vs_reference"
if (!coef %in% resultsNames(dds)) stop("coefficient ", coef, " not in model: ",
                                       paste(resultsNames(dds), collapse = ", "))
res <- results(dds, name = coef, alpha = fdr)
shr <- lfcShrink(dds, coef = coef, type = "apeglm", quiet = TRUE)

out <- data.frame(
  gene_id   = rownames(res),
  gene_name = ann$gene_name[match(rownames(res), ann$gene_id)],
  gene_type = ann$gene_type[match(rownames(res), ann$gene_id)],
  baseMean  = res$baseMean,
  log2FoldChange = res$log2FoldChange,
  lfcSE     = res$lfcSE,
  stat      = res$stat,
  pvalue    = res$pvalue,
  padj      = res$padj,
  log2FoldChange_shrunk = shr$log2FoldChange,
  stringsAsFactors = FALSE)
out <- out[order(out$pvalue, na.last = TRUE), ]
write.table(out, sprintf("%s.deseq2.tsv", id), sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")

nc <- counts(dds, normalized = TRUE)
write.table(data.frame(gene_id = rownames(nc), nc, check.names = FALSE),
            sprintf("%s.normalized_counts.tsv", id), sep = "\t", quote = FALSE,
            row.names = FALSE)
sf <- sizeFactors(dds)
write.table(data.frame(sample = names(sf), condition = cd[names(sf), "condition"],
                       size_factor = sf),
            sprintf("%s.size_factors.tsv", id), sep = "\t", quote = FALSE,
            row.names = FALSE)

# ---------------------------------------------------------------- QC plots
vsd <- if (nrow(dds) >= 1000) vst(dds, blind = TRUE) else
  varianceStabilizingTransformation(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
g <- ggplot(pca, aes(PC1, PC2, colour = condition, label = name)) +
  geom_point(size = 3) + geom_text(vjust = -1, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = c(reference = "#4C72B0", target = "#C44E52"),
                      labels = c(reference = row$reference, target = row$target)) +
  xlab(sprintf("PC1: %d%% variance", pv[1])) + ylab(sprintf("PC2: %d%% variance", pv[2])) +
  ggtitle(sprintf("%s (vst, blind)", id)) + theme_bw()
ggsave(sprintf("%s.pca.png", id), g, width = 6, height = 5, dpi = 120)

png(sprintf("%s.dispersion.png", id), width = 700, height = 600)
plotDispEsts(dds, main = sprintf("%s dispersion estimates", id))
invisible(dev.off())

png(sprintf("%s.ma.png", id), width = 700, height = 600)
plotMA(shr, alpha = fdr, main = sprintf("%s (apeglm-shrunken LFC)", id))
invisible(dev.off())

d <- as.matrix(dist(t(assay(vsd))))
png(sprintf("%s.sample_distances.png", id), width = 700, height = 650)
heatmap(d, symm = TRUE, margins = c(10, 10), main = sprintf("%s sample distances", id))
invisible(dev.off())

# ---------------------------------------------------------------- label swaps
real <- n_de(res)
swaps <- data.frame()
if (do_swaps) {
  n  <- ncol(m)
  nt <- sum(cd$condition == "target")
  truth <- which(cd$condition == "target")
  combos <- combn(n, nt, simplify = FALSE)
  seen <- character(0)
  for (cmb in combos) {
    comp <- setdiff(seq_len(n), cmb)
    k1 <- paste(cmb, collapse = ","); k2 <- paste(comp, collapse = ",")
    if (identical(sort(cmb), truth) || (length(comp) == nt && identical(sort(comp), truth))) next
    if (k1 %in% seen || k2 %in% seen) next
    seen <- c(seen, k1)
    fake <- cd
    fake$condition <- factor(ifelse(seq_len(n) %in% cmb, "target", "reference"),
                             levels = c("reference", "target"))
    r <- tryCatch(results(fit(m[keep, , drop = FALSE], fake), name = coef, alpha = fdr),
                  error = function(e) NULL)
    if (is.null(r)) next   # relabelling confounded with a covariate
    x <- n_de(r)
    swaps <- rbind(swaps, data.frame(
      fake_target = paste(rownames(cd)[cmb], collapse = ","),
      true_targets_in_fake_target = sum(cmb %in% truth),
      n_up = x[["up"]], n_down = x[["down"]], n_total = sum(x)))
    if (nrow(swaps) >= 30) break
  }
  write.table(swaps, sprintf("%s.label_swaps.tsv", id), sep = "\t", quote = FALSE,
              row.names = FALSE)
}
swap_max <- if (nrow(swaps)) max(swaps$n_total) else NA
swap_med <- if (nrow(swaps)) median(swaps$n_total) else NA
if (!is.na(swap_max) && sum(real) <= swap_max)
  message(sprintf("WARNING [%s]: %d DE genes, but a label swap gave %d. The result is not ",
                  id, sum(real), swap_max), "clearly separable from relabelling noise.")

summ <- data.frame(
  id = id, cell_line = row$cell_line, reference = row$reference, target = row$target,
  design = row$design, n_reference = sum(cd$condition == "reference"),
  n_target = sum(cd$condition == "target"), genes_tested = nrow(dds),
  fdr = fdr, lfc_threshold = lfc_thr,
  n_up = real[["up"]], n_down = real[["down"]],
  n_sig_any_lfc = sum(!is.na(res$padj) & res$padj < fdr),
  label_swaps = nrow(swaps), swap_max_de = swap_max, swap_median_de = swap_med,
  size_factor_min = round(min(sf), 3), size_factor_max = round(max(sf), 3))
write.table(summ, sprintf("%s.deseq2_summary.tsv", id), sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")
message(sprintf("[%s] up %d, down %d (padj < %g, |LFC| >= %g); label swaps max %s",
                id, real[["up"]], real[["down"]], fdr, lfc_thr, swap_max))

writeLines(capture.output(sessionInfo()), sprintf("%s.sessionInfo.txt", id))
