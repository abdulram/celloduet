#!/usr/bin/env Rscript
# 11a - TGF-beta target panel: per-gene DESeq2 effects for both contrasts under every model.
#
# Panel (user-supplied, with confidence weight): Serpine1, Smad7, Tgfbi, Pmepa1 (High); Skil, Ccn2 (=Ctgf), Fn1
# (Medium); Smurf2 (Low). Context genes (pathway input, not readouts): Tgfb1-3, Tgfbr1-2, Ltbp1.
# Models (all 18 animals, full transcriptome so size factors / dispersions are genome-wide):
#   baseline ~ sex + group ; comp1 ~ sex + Composition_PC1 + group ; comp13 ~ sex + Composition_PC1 + Composition_PC3 + group
#   tech ~ sex + pct_mito + group
# Contrasts: 5xFAD_NTC vs WT_NTC (genotype) and 5xFAD_SPP1 vs 5xFAD_NTC (treatment, same library profile).

suppressPackageStartupMessages({library(DESeq2)})
set.seed(7)
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
CA <- normalizePath(file.path(here, "..")); OUT <- file.path(CA, "11_tgfb_panel"); TAB <- file.path(CA, "tables")
source(file.path(here, "de_helpers.R"))

obj <- readRDS(file.path(CA, "01_QC", "qc_objects.rds"))
meta <- read.delim(file.path(TAB, "sample_metadata_with_composition.tsv"), row.names = 1, check.names = FALSE)
meta$group <- factor(meta$group, levels = c("WT_NTC", "5xFAD_NTC", "5xFAD_SPP1")); meta$sex <- factor(meta$sex)
for (v in c("Composition_PC1", "Composition_PC3", "pct_mito")) meta[[v]] <- as.numeric(scale(meta[[v]]))
genes <- obj$genes
PANEL <- c("Serpine1", "Smad7", "Tgfbi", "Pmepa1", "Skil", "Smurf2", "Ccn2", "Fn1", "Tgfb1", "Tgfb2", "Tgfb3", "Tgfbr1", "Tgfbr2", "Ltbp1")
ids <- rownames(genes)[genes$gene_name %in% PANEL & rownames(genes) %in% rownames(obj$counts)]

MODELS <- list(baseline = ~ sex + group, comp1 = ~ sex + Composition_PC1 + group,
               comp13 = ~ sex + Composition_PC1 + Composition_PC3 + group, tech = ~ sex + pct_mito + group)
CONTRASTS <- list(genotype_5xFAD_vs_WT = c("group", "5xFAD_NTC", "WT_NTC"), SPP1_vs_NTC = c("group", "5xFAD_SPP1", "5xFAD_NTC"))
rows <- list()
for (m in names(MODELS)) {
  fit <- run_deseq(obj$counts[, rownames(meta)], meta, MODELS[[m]])
  for (k in names(CONTRASTS)) {
    r <- results_tbl(fit, CONTRASTS[[k]], genes)
    r <- r[r$gene_id %in% ids, ]
    rows[[paste(m, k)]] <- data.frame(model = m, contrast = k, r[, c("gene_id", "gene_name", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")])
  }
}
res <- do.call(rbind, rows)
write.table(res, file.path(TAB, "tgfb_panel_DE_all_models.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(res[order(res$contrast, res$gene_name, res$model), c("contrast", "gene_name", "model", "log2FoldChange", "pvalue", "padj")], digits = 3)
