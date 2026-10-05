#!/usr/bin/env Rscript
# 05 - Composition-adjusted genotype models (parsimonious) + collinearity diagnostics + comparison to baseline.
#
# Models (all 18 animals, genotype contrast = 5xFAD_NTC vs WT_NTC):
#   baseline : ~ sex + group                                   (step 04)
#   comp1    : ~ sex + Composition_PC1 + group
#   comp12   : ~ sex + Composition_PC1 + Composition_PC2 + group
#   comp13   : ~ sex + Composition_PC1 + Composition_PC3 + group  (PC3 = bone axis, also genotype-associated)
#   tech     : ~ sex + pct_mito + group                         (known technical covariate; sensitivity)
#   comp1tech: ~ sex + Composition_PC1 + pct_mito + group
# Covariates are standardised. Rank and VIF are checked for every model.
# Gene classification (baseline padj<0.05 genes; retention = adjusted LFC / baseline LFC):
#   robust                 retention >= 0.5 (same direction, at least half the effect kept)
#   composition-sensitive  retention < 0.5 (incl. sign reversal)
#   potentially masked     baseline padj >= 0.05, adjusted padj < 0.05 and |adjusted LFC| > |baseline LFC|
# NOTE: retention in a model whose covariate is nearly collinear with genotype is uninformative (VIF reported).

suppressPackageStartupMessages({library(DESeq2); library(limma); library(ggplot2); library(car)})
set.seed(7)
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
CA <- normalizePath(file.path(here, "..")); OUT <- file.path(CA, "05_composition_adjusted_DE"); TAB <- file.path(CA, "tables")
source(file.path(here, "common.R")); source(file.path(here, "de_helpers.R"))

obj <- readRDS(file.path(CA, "01_QC", "qc_objects.rds"))
meta <- read.delim(file.path(TAB, "sample_metadata_with_composition.tsv"), row.names = 1, check.names = FALSE)
meta$group <- factor(meta$group, levels = c("WT_NTC", "5xFAD_NTC", "5xFAD_SPP1")); meta$sex <- factor(meta$sex)
for (v in c("Composition_PC1", "Composition_PC2", "Composition_PC3", "pct_mito")) meta[[v]] <- as.numeric(scale(meta[[v]]))
cf <- obj$counts[, rownames(meta)]; genes <- obj$genes
base <- read.delim(file.path(TAB, "DE_baseline_5xFAD_vs_WT.tsv"))

MODELS <- list(comp1 = ~ sex + Composition_PC1 + group, comp12 = ~ sex + Composition_PC1 + Composition_PC2 + group,
               comp13 = ~ sex + Composition_PC1 + Composition_PC3 + group, tech = ~ sex + pct_mito + group,
               comp1tech = ~ sex + Composition_PC1 + pct_mito + group)

# ---- collinearity diagnostics
ntc <- meta$group %in% c("WT_NTC", "5xFAD_NTC")
g5 <- as.numeric(meta$genotype == "5xFAD")
diag <- do.call(rbind, lapply(names(MODELS), function(m) {
  f <- MODELS[[m]]; X <- model.matrix(f, meta)
  dummy <- lm(update(f, rnorm(nrow(meta)) ~ .), data = meta)
  v <- car::vif(dummy); v <- if (is.matrix(v)) v[, "GVIF^(1/(2*Df))"]^2 else v
  covs <- setdiff(all.vars(f), c("sex", "group"))
  data.frame(model = m, formula = deparse(f), rank = qr(X)$rank, n_coef = ncol(X),
             max_VIF = round(max(v), 2), VIF_group = round(v[["group"]], 2),
             covariates = paste(covs, collapse = "+"),
             cor_with_genotype_NTC = paste(sapply(covs, function(cv) round(cor(meta[ntc, cv], g5[ntc]), 2)), collapse = ";"))
}))
write.table(diag, file.path(TAB, "adjusted_models_collinearity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(diag)

# ---- fit models
all_res <- list(baseline = base)
for (m in names(MODELS)) {
  fit <- run_deseq(cf, meta, MODELS[[m]])
  r <- results_tbl(fit, c("group", "5xFAD_NTC", "WT_NTC"), genes)
  write.table(r, file.path(TAB, paste0("DE_", m, "_5xFAD_vs_WT.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  all_res[[m]] <- r
  cat(m, "padj<0.05:", sum(r$padj < 0.05, na.rm = TRUE), "\n")
  cam <- camera_sets(obj$vst, meta, MODELS[[m]], "group5xFAD_NTC", file.path(TAB, "gene_sets.tsv"))
  write.table(cam, file.path(TAB, paste0("pathways_camera_", m, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}

# ---- per-gene comparison table + classification
cmp <- base[, c("gene_id", "gene_name", "baseMean", "log2FoldChange", "padj")]
names(cmp)[4:5] <- c("LFC_baseline", "padj_baseline")
for (m in names(MODELS)) {
  r <- all_res[[m]][, c("gene_id", "log2FoldChange", "padj")]; names(r)[2:3] <- paste0(c("LFC_", "padj_"), m)
  cmp <- merge(cmp, r, by = "gene_id", all.x = TRUE)
  ret <- cmp[[paste0("LFC_", m)]] / cmp$LFC_baseline
  cls <- rep("not genotype-associated at baseline", nrow(cmp))
  bs <- !is.na(cmp$padj_baseline) & cmp$padj_baseline < 0.05
  cls[bs & ret >= 0.5] <- "robust (>=50% of effect retained)"
  cls[bs & (ret < 0.5 | is.na(ret))] <- "composition/covariate-sensitive"
  mk <- !bs & !is.na(cmp[[paste0("padj_", m)]]) & cmp[[paste0("padj_", m)]] < 0.05 & abs(cmp[[paste0("LFC_", m)]]) > abs(cmp$LFC_baseline)
  cls[mk] <- "potentially masked"
  cmp[[paste0("retention_", m)]] <- ret; cmp[[paste0("class_", m)]] <- cls
}
cmp <- cmp[order(cmp$padj_baseline), ]
write.table(cmp, gzfile(file.path(TAB, "genotype_effects_baseline_vs_adjusted.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
summ <- do.call(rbind, lapply(names(MODELS), function(m) {
  t <- table(factor(cmp[[paste0("class_", m)]], levels = c("robust (>=50% of effect retained)", "composition/covariate-sensitive", "potentially masked")))
  bs <- !is.na(cmp$padj_baseline) & cmp$padj_baseline < 0.05
  data.frame(model = m, baseline_sig = sum(bs), robust = t[[1]], sensitive = t[[2]], masked = t[[3]],
             adjusted_sig = sum(cmp[[paste0("padj_", m)]] < 0.05, na.rm = TRUE),
             retained_and_sig = sum(bs & cmp[[paste0("padj_", m)]] < 0.05 & cmp[[paste0("retention_", m)]] >= 0.5, na.rm = TRUE),
             median_retention = round(median(cmp[[paste0("retention_", m)]][bs], na.rm = TRUE), 2),
             LFC_correlation = round(cor(cmp$LFC_baseline[bs], cmp[[paste0("LFC_", m)]][bs], use = "complete"), 2))
}))
write.table(summ, file.path(TAB, "adjusted_models_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ)

# ---- Figure 7: baseline vs adjusted LFC
long <- do.call(rbind, lapply(names(MODELS), function(m) data.frame(model = m, base = cmp$LFC_baseline, adj = cmp[[paste0("LFC_", m)]],
                                                                    cls = cmp[[paste0("class_", m)]])))
long <- long[long$cls != "not genotype-associated at baseline", ]
long$model <- factor(long$model, levels = names(MODELS),
                     labels = sapply(names(MODELS), function(m) paste0(m, "  (VIF group=", diag$VIF_group[diag$model == m], ")")))
g <- ggplot(long, aes(base, adj, colour = cls)) + geom_hline(yintercept = 0, colour = INK2) + geom_vline(xintercept = 0, colour = INK2) +
  geom_abline(slope = 1, linetype = 3, colour = INK2) + geom_point(size = 0.5, alpha = 0.6) + facet_wrap(~model, ncol = 3) +
  scale_colour_manual(values = c("robust (>=50% of effect retained)" = "#1baf7a", "composition/covariate-sensitive" = "#eb6834", "potentially masked" = "#4a3aa7")) +
  coord_cartesian(xlim = c(-4, 4), ylim = c(-4, 4)) + labs(x = "baseline genotype log2FC", y = "adjusted genotype log2FC", colour = NULL,
  title = "Genotype effects before vs after adjustment (baseline-significant genes)") + theme_ca() + theme(legend.position = "bottom")
ggsave(file.path(OUT, "Fig7_baseline_vs_adjusted_LFC.png"), g, width = 12, height = 8.5, dpi = 150)
sl <- reshape(summ[, c("model", "robust", "sensitive", "masked")], direction = "long", varying = c("robust", "sensitive", "masked"),
              v.names = "n", timevar = "class", times = c("robust", "sensitive", "masked"))
g2 <- ggplot(sl, aes(model, n, fill = class)) + geom_col(position = "dodge", colour = SURFACE) +
  scale_fill_manual(values = c(robust = "#1baf7a", sensitive = "#eb6834", masked = "#4a3aa7")) +
  labs(title = "Gene categories per adjustment model", y = "genes", x = NULL) + theme_ca()
ggsave(file.path(OUT, "significance_categories.png"), g2, width = 8, height = 4.5, dpi = 150)
writeLines(capture.output(sessionInfo()), file.path(CA, "logs", "sessionInfo_05.txt"))
