#!/usr/bin/env Rscript
# 12a - PRIMARY treatment contrast: 5xFAD_SPP1 vs 5xFAD_NTC, genome-wide, under every model.
# SPP1 and NTC libraries share the same technical profile, so this contrast is not affected by the
# genotype/processing confound. Models (all 18 animals for dispersion; contrast SPP1 vs NTC):
#   baseline ~ sex + group ; comp1 ~ sex + Composition_PC1 + group ; comp13 ~ sex + Composition_PC1 + Composition_PC3 + group
#   tech ~ sex + pct_mito + group
# plus a 5xFAD-only model (n = 12) ~ sex + treatment, and limma camera pathway tests on the 5xFAD-only design.

suppressPackageStartupMessages({library(DESeq2); library(limma)})
set.seed(7)
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
CA <- normalizePath(file.path(here, "..")); OUT <- file.path(CA, "12_spp1_vs_ntc"); TAB <- file.path(CA, "tables")
dir.create(OUT, showWarnings = FALSE)
source(file.path(here, "common.R")); source(file.path(here, "de_helpers.R"))

obj <- readRDS(file.path(CA, "01_QC", "qc_objects.rds"))
meta <- read.delim(file.path(TAB, "sample_metadata_with_composition.tsv"), row.names = 1, check.names = FALSE)
meta$group <- factor(meta$group, levels = c("WT_NTC", "5xFAD_NTC", "5xFAD_SPP1")); meta$sex <- factor(meta$sex)
for (v in c("Composition_PC1", "Composition_PC3", "pct_mito")) meta[[v]] <- as.numeric(scale(meta[[v]]))
genes <- obj$genes; cf <- obj$counts[, rownames(meta)]

MODELS <- list(baseline = ~ sex + group, comp1 = ~ sex + Composition_PC1 + group,
               comp13 = ~ sex + Composition_PC1 + Composition_PC3 + group, tech = ~ sex + pct_mito + group)
summ <- list(); all <- NULL
for (m in names(MODELS)) {
  fit <- run_deseq(cf, meta, MODELS[[m]])
  r <- results_tbl(fit, c("group", "5xFAD_SPP1", "5xFAD_NTC"), genes)
  write.table(r, file.path(TAB, paste0("DE_SPP1_vs_NTC_", m, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  summ[[m]] <- data.frame(model = m, formula = deparse(MODELS[[m]]), padj05 = sum(r$padj < 0.05, na.rm = TRUE),
                          padj10 = sum(r$padj < 0.1, na.rm = TRUE), up05 = sum(r$padj < 0.05 & r$log2FoldChange > 0, na.rm = TRUE),
                          down05 = sum(r$padj < 0.05 & r$log2FoldChange < 0, na.rm = TRUE))
  x <- r[, c("gene_id", "gene_name", "baseMean", "log2FoldChange", "pvalue", "padj")]
  names(x)[4:6] <- paste0(c("LFC_", "p_", "padj_"), m)
  all <- if (is.null(all)) x else merge(all, x[, c(1, 4:6)], by = "gene_id")
}
# 5xFAD-only model
f5 <- meta[meta$genotype == "5xFAD", ]; f5$treatment <- factor(f5$treatment, levels = c("NTC", "SPP1"))
fit5 <- run_deseq(cf[, rownames(f5)], f5, ~ sex + treatment)
r5 <- results_tbl(fit5, NULL, genes, name = "treatment_SPP1_vs_NTC")
write.table(r5, file.path(TAB, "DE_SPP1_vs_NTC_5xFADonly.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
summ[["5xFADonly"]] <- data.frame(model = "5xFADonly", formula = "~ sex + treatment (12 animals)", padj05 = sum(r5$padj < 0.05, na.rm = TRUE),
                                  padj10 = sum(r5$padj < 0.1, na.rm = TRUE), up05 = sum(r5$padj < 0.05 & r5$log2FoldChange > 0, na.rm = TRUE),
                                  down05 = sum(r5$padj < 0.05 & r5$log2FoldChange < 0, na.rm = TRUE))
x <- r5[, c("gene_id", "log2FoldChange", "pvalue", "padj")]; names(x)[2:4] <- c("LFC_5xFADonly", "p_5xFADonly", "padj_5xFADonly")
all <- merge(all, x, by = "gene_id")
all$min_padj_any_model <- apply(all[, grep("^padj_", names(all))], 1, min, na.rm = TRUE)
all$n_models_padj05 <- rowSums(all[, grep("^padj_", names(all))] < 0.05, na.rm = TRUE)
all$sign_consistent <- apply(sign(all[, grep("^LFC_", names(all))]), 1, function(s) length(unique(s)) == 1)
all <- all[order(all$p_baseline), ]
write.table(all, gzfile(file.path(TAB, "SPP1_vs_NTC_all_models.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
s <- do.call(rbind, summ); write.table(s, file.path(TAB, "SPP1_vs_NTC_DE_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(s)
print(head(all[all$n_models_padj05 > 0, c("gene_name", "baseMean", "LFC_baseline", "padj_baseline", "padj_comp1", "padj_comp13", "padj_tech", "padj_5xFADonly", "n_models_padj05")], 40), digits = 3)

# camera on the 5xFAD-only design
v5 <- obj$vst[, rownames(f5)]
cam <- camera_sets(v5, f5, ~ sex + treatment, "treatmentSPP1", file.path(TAB, "gene_sets.tsv"))
write.table(cam, file.path(TAB, "pathways_camera_SPP1_vs_NTC.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(head(cam, 12))
ma_volcano(r5, OUT, "SPP1_vs_NTC_5xFADonly", "5xFAD_SPP1 vs 5xFAD_NTC (~ sex + treatment)")
writeLines(capture.output(sessionInfo()), file.path(CA, "logs", "sessionInfo_12a.txt"))
