#!/usr/bin/env Rscript
# 13a - Genome-wide SPP1 vs NTC within 5xFAD, adjusted for fibroblast abundance.
# Model (12 5xFAD animals): ~ sex + fibroblast_stromal + treatment, where fibroblast_stromal is the validated
# fibroblast marker-signature score (step 02), standardised. Rationale: SPP1 samples carry less fibroblast RNA
# (step 12b, perm p = 0.025), so this model asks whether SPP1 changes expression beyond fibroblast content.
# NOTE: genes that are themselves fibroblast-signature markers (Pdgfra, Dcn, Fbln1, Mfap4, Igfbp6, Slc47a1,
# Col14a1, Smoc2) are partly adjusted for themselves (circular) - flagged in the output.

suppressPackageStartupMessages({library(DESeq2)})
set.seed(7)
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
CA <- normalizePath(file.path(here, "..")); TAB <- file.path(CA, "tables")
source(file.path(here, "common.R")); source(file.path(here, "de_helpers.R"))
obj <- readRDS(file.path(CA, "01_QC", "qc_objects.rds"))
meta <- read.delim(file.path(TAB, "sample_metadata_with_composition.tsv"), row.names = 1, check.names = FALSE)
sig <- read.delim(file.path(TAB, "composition_signature_scores_all.tsv"), row.names = 1, check.names = FALSE)
f5 <- meta[meta$genotype == "5xFAD", ]
f5$treatment <- factor(f5$treatment, levels = c("NTC", "SPP1")); f5$sex <- factor(f5$sex)
f5$fibro <- as.numeric(scale(sig[rownames(f5), "fibroblast_stromal"]))
fit <- run_deseq(obj$counts[, rownames(f5)], f5, ~ sex + fibro + treatment)
r <- results_tbl(fit, NULL, obj$genes, name = "treatment_SPP1_vs_NTC")
markers <- c("Pdgfra", "Dcn", "Fbln1", "Mfap4", "Igfbp6", "Slc47a1", "Col14a1", "Smoc2")
r$fibroblast_signature_marker <- r$gene_name %in% markers
write.table(r, file.path(TAB, "DE_SPP1_vs_NTC_fibroblast_adjusted.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("fibroblast-adjusted padj<0.05:", sum(r$padj < 0.05, na.rm = TRUE), "\n")
print(head(r[, c("gene_name", "log2FoldChange", "pvalue", "padj")], 15), digits = 3)
