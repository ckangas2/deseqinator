#!/usr/bin/env Rscript
# Preranked GSEA (fgsea) on one DESeq2 contrast.
#
# Ranking: the Wald statistic, not log2FoldChange. Raw LFC puts noisy
# low-count genes at the extremes of the list; the statistic already weighs
# effect size by its standard error.
#
# Three checks that each catch a run that would otherwise "find nothing" or
# "find something" for the wrong reason:
#   1. Symbol mapping and duplicates are counted and reported. Several
#      Ensembl IDs share a symbol (PAR genes, readthroughs); fgsea needs
#      unique names, and the rule here is: keep the entry with the largest
#      |stat|.
#   2. Gene-set coverage. If the ranked list is keyed by Ensembl IDs while
#      the sets are symbols (or human sets meet mouse genes), every set
#      drops below minSize and fgsea returns an empty table, which reads as
#      "no enrichment". Median coverage below --min_coverage stops the run.
#   3. Rank-shuffle control. The same sets against the same statistics
#      randomly reassigned to genes. Anything significant there means the
#      significance in the real run is not coming from the ranking.

suppressPackageStartupMessages({
  library(fgsea)
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
opt <- list()
i <- 1
while (i <= length(args)) {
  k <- sub("^--", "", args[i]); opt[[k]] <- args[i + 1]; i <- i + 2
}
need <- c("deseq", "id", "gene_sets", "gmt", "min_size", "max_size",
          "min_coverage", "max_null_hits", "fdr", "seed")
miss <- setdiff(need, names(opt))
if (length(miss)) stop("missing arguments: ", paste(miss, collapse = ", "))
id        <- opt$id
min_size  <- as.integer(opt$min_size)
max_size  <- as.integer(opt$max_size)
min_cov   <- as.numeric(opt$min_coverage)
max_null  <- as.integer(opt$max_null_hits)
fdr       <- as.numeric(opt$fdr)
seed      <- as.integer(opt$seed)
log_lines <- character(0)
say <- function(...) { s <- sprintf(...); message(s); log_lines <<- c(log_lines, s) }

# ---------------------------------------------------------------- ranks
res <- fread(opt$deseq, na.strings = "NA")
n_all <- nrow(res)
res <- res[!is.na(stat)]
n_stat <- nrow(res)
res <- res[!is.na(gene_name) & gene_name != ""]
n_named <- nrow(res)
res[, abs_stat := abs(stat)]
setorder(res, -abs_stat)
dups <- res[duplicated(gene_name)]
res <- res[!duplicated(gene_name)]
stats <- setNames(res$stat, res$gene_name)
stats <- sort(stats, decreasing = TRUE)
say("[%s] ranked list: %d genes tested, %d with a statistic, %d with a symbol, %d after dropping %d duplicate symbols (kept largest |stat|)",
    id, n_all, n_stat, n_named, length(stats), nrow(dups))

# ---------------------------------------------------------------- sets
spec <- opt$gene_sets
if (startsWith(spec, "msigdb:")) {
  collection <- sub("^msigdb:", "", spec)
  suppressPackageStartupMessages(library(msigdbr))
  df <- if ("collection" %in% names(formals(msigdbr))) {
    msigdbr(species = "Homo sapiens", collection = collection)
  } else {
    msigdbr(species = "Homo sapiens", category = collection)
  }
  pathways <- lapply(split(df$gene_symbol, df$gs_name), unique)
  source_desc <- sprintf("msigdbr %s, collection %s", packageVersion("msigdbr"), collection)
} else {
  pathways <- gmtPathways(opt$gmt)
  source_desc <- sprintf("GMT %s", basename(opt$gmt))
}
say("[%s] %d gene sets from %s", id, length(pathways), source_desc)

coverage <- vapply(pathways, function(g) mean(unique(g) %in% names(stats)), numeric(1))
med_cov <- median(coverage)
say("[%s] gene-set coverage: median %.1f%% of set genes present in the ranked list",
    id, 100 * med_cov)
fwrite(data.table(pathway = names(coverage), size = lengths(pathways), coverage = coverage),
       sprintf("%s.set_coverage.tsv", id), sep = "\t")
if (med_cov < min_cov) {
  writeLines(log_lines, sprintf("%s.fgsea.log", id))
  stop(sprintf(paste0(
    "median gene-set coverage %.1f%% is below %.0f%%. The ranked list and the gene sets ",
    "are probably keyed differently (Ensembl IDs vs symbols, or the wrong species). ",
    "Ranked names look like: %s. Set genes look like: %s."),
    100 * med_cov, 100 * min_cov,
    paste(head(names(stats), 3), collapse = ", "),
    paste(head(pathways[[1]], 3), collapse = ", ")))
}

# ---------------------------------------------------------------- fgsea
run <- function(s, sd) {
  set.seed(sd)
  fgsea(pathways = pathways, stats = s, minSize = min_size, maxSize = max_size,
        eps = 0, nproc = 1)
}
fg <- run(stats, seed)
setorder(fg, padj, -NES)
fg_out <- copy(fg)
fg_out[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
fwrite(fg_out, sprintf("%s.fgsea.tsv", id), sep = "\t", na = "NA")

shuf <- setNames(sample(unname(stats)), names(stats))
shuf <- sort(shuf, decreasing = TRUE)
fg_null <- run(shuf, seed + 1L)
null_hits <- sum(!is.na(fg_null$padj) & fg_null$padj < fdr)
fg_null[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
fwrite(fg_null, sprintf("%s.fgsea_rank_shuffle.tsv", id), sep = "\t", na = "NA")

sig <- fg[!is.na(padj) & padj < fdr]
say("[%s] %d sets tested, %d with padj < %g (%d up, %d down); rank-shuffle control: %d",
    id, nrow(fg), nrow(sig), fdr, sum(sig$NES > 0), sum(sig$NES < 0), null_hits)

# ---------------------------------------------------------------- plots
top <- head(fg[!is.na(padj)][order(padj)], 20)
if (nrow(top)) {
  top[, label := sub("^HALLMARK_", "", pathway)]
  top[, direction := ifelse(NES > 0, "up", "down")]
  top[, label := factor(label, levels = label[order(NES)])]
  g <- ggplot(top, aes(NES, label, fill = direction, alpha = padj < fdr)) +
    geom_col() +
    scale_fill_manual(values = c(up = "#C44E52", down = "#4C72B0")) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), name = sprintf("padj < %g", fdr)) +
    labs(x = "normalised enrichment score", y = NULL,
         title = sprintf("%s: top gene sets by padj", id)) +
    theme_bw()
  ggsave(sprintf("%s.fgsea_top.png", id), g, width = 8, height = 0.3 * nrow(top) + 1.5, dpi = 120)
}
for (dir in c("up", "down")) {
  pick <- if (dir == "up") sig[NES > 0][order(padj)] else sig[NES < 0][order(padj)]
  if (nrow(pick)) {
    p <- pick$pathway[1]
    g <- plotEnrichment(pathways[[p]], stats) +
      labs(title = sprintf("%s: %s", id, p),
           subtitle = sprintf("NES %.2f, padj %.2g", pick$NES[1], pick$padj[1]))
    ggsave(sprintf("%s.enrichment_%s.png", id, dir), g, width = 7, height = 4, dpi = 120)
  }
}

top_up   <- paste(head(sig[NES > 0][order(padj)]$pathway, 3), collapse = ";")
top_down <- paste(head(sig[NES < 0][order(padj)]$pathway, 3), collapse = ";")
fwrite(data.table(id = id, gene_sets = source_desc, ranked_genes = length(stats),
                  duplicate_symbols_dropped = nrow(dups), unnamed_dropped = n_stat - n_named,
                  median_set_coverage = round(med_cov, 4), sets_tested = nrow(fg),
                  sets_sig = nrow(sig), sets_sig_up = sum(sig$NES > 0),
                  sets_sig_down = sum(sig$NES < 0), rank_shuffle_hits = null_hits,
                  top_up = top_up, top_down = top_down),
       sprintf("%s.fgsea_summary.tsv", id), sep = "\t")
writeLines(log_lines, sprintf("%s.fgsea.log", id))
writeLines(capture.output(sessionInfo()), sprintf("%s.fgsea.sessionInfo.txt", id))

if (null_hits > max_null) {
  stop(sprintf(paste0(
    "rank-shuffle control gave %d significant sets (max %d). Significance is not coming ",
    "from the ranking; inspect %s.fgsea_rank_shuffle.tsv."), null_hits, max_null, id))
}
