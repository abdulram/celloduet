#!/usr/bin/env Rscript
# 04 - Baseline (UNADJUSTED) genotype association.
#
# Model: ~ sex + group  (group = WT_NTC / 5xFAD_NTC / 5xFAD_SPP1), all 18 animals (SPP1 animals help estimate
# dispersion). Genotype contrast = 5xFAD_NTC vs WT_NTC (same treatment). Age and batch are unavailable and
# therefore cannot be modelled. Pathway test: limma::camera on blind-VST with the same design (competitive test
# that accounts for inter-gene correlation; sample-level design rather than gene-permutation GSEA).
# This is the UNADJUSTED genotype association, NOT a cell-intrinsic disease effect.

suppressPackageStartupMessages({library(DESeq2); library(limma); library(ggplot2); library(pheatmap)})
set.seed(7)
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
CA <- normalizePath(file.path(here, "..")); OUT <- file.path(CA, "04_baseline_DE"); TAB <- file.path(CA, "tables")
source(file.path(here, "common.R")); source(file.path(here, "de_helpers.R"))

obj <- readRDS(file.path(CA, "01_QC", "qc_objects.rds"))
meta <- read.delim(file.path(TAB, "sample_metadata_with_composition.tsv"), row.names = 1, check.names = FALSE)
meta$group <- factor(meta$group, levels = c("WT_NTC", "5xFAD_NTC", "5xFAD_SPP1")); meta$sex <- factor(meta$sex)
cf <- obj$counts[, rownames(meta)]; genes <- obj$genes

fit <- run_deseq(cf, meta, ~ sex + group)
res <- results_tbl(fit, c("group", "5xFAD_NTC", "WT_NTC"), genes)
write.table(res, file.path(TAB, "DE_baseline_5xFAD_vs_WT.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("baseline padj<0.05:", sum(res$padj < 0.05, na.rm = TRUE), " up:", sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE),
    " down:", sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE), "\n")
saveRDS(fit, file.path(OUT, "dds_baseline.rds"))

ma_volcano(res, OUT, "baseline", "Baseline ~ sex + group: 5xFAD_NTC vs WT_NTC")
top_heatmap(obj$vst, res, meta, OUT, "baseline", genes)
cam <- camera_sets(obj$vst, meta, ~ sex + group, "group5xFAD_NTC", file.path(TAB, "gene_sets.tsv"))
write.table(cam, file.path(TAB, "pathways_camera_baseline.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(head(cam, 15))
writeLines(capture.output(sessionInfo()), file.path(CA, "logs", "sessionInfo_04.txt"))
