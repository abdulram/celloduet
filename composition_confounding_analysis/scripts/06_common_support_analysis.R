#!/usr/bin/env Rscript
# 06 - Composition-matched / common-support sensitivity analyses (pre-specified algorithms; no hand-picking).
#
# Pool: WT_NTC (n=6) vs 5xFAD_NTC (n=6).
# A. RESTRICTION: keep animals whose Composition_PC1 AND Composition_PC3 (the two genotype-associated composition
#    axes, step 03) both lie inside the common-support range [max of group minima, min of group maxima].
# B. MATCHING: MatchIt nearest-neighbour 1:1 without replacement, Mahalanobis distance on Composition_PC1 +
#    Composition_PC3, caliper 1.0 SD (applied per covariate). Unmatched animals are dropped.
# C. TECHNICAL COMMON SUPPORT: only animals with the high-multimapping library profile (tech_profile B):
#    WT 616 + 623 vs all 6 5xFAD_NTC. Composition is NOT matched here; this isolates the technical axis.
# Each cohort: DESeq2 ~ sex + genotype (sex dropped if not estimable). Compared with baseline by LFC correlation
# and direction concordance among baseline-significant genes, plus camera pathway tests.

suppressPackageStartupMessages({library(DESeq2); library(MatchIt); library(ggplot2); library(limma)})
set.seed(7)
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
CA <- normalizePath(file.path(here, "..")); OUT <- file.path(CA, "06_common_support_matching"); TAB <- file.path(CA, "tables")
source(file.path(here, "common.R")); source(file.path(here, "de_helpers.R"))

obj <- readRDS(file.path(CA, "01_QC", "qc_objects.rds"))
meta <- read.delim(file.path(TAB, "sample_metadata_with_composition.tsv"), row.names = 1, check.names = FALSE)
pool <- meta[meta$group %in% c("WT_NTC", "5xFAD_NTC"), ]
pool$genotype <- factor(pool$genotype, levels = c("WT", "5xFAD")); pool$sex <- factor(pool$sex)
pool$treat <- as.integer(pool$genotype == "5xFAD")
base <- read.delim(file.path(TAB, "DE_baseline_5xFAD_vs_WT.tsv"))
bsig <- base$gene_id[!is.na(base$padj) & base$padj < 0.05]

# ---- A. restriction
cs <- function(v) c(max(tapply(pool[[v]], pool$genotype, min)), min(tapply(pool[[v]], pool$genotype, max)))
r1 <- cs("Composition_PC1"); r3 <- cs("Composition_PC3")
inA <- pool$Composition_PC1 >= r1[1] & pool$Composition_PC1 <= r1[2] & pool$Composition_PC3 >= r3[1] & pool$Composition_PC3 <= r3[2]
coh <- list(restricted = rownames(pool)[inA])

# ---- B. matching
m <- matchit(treat ~ Composition_PC1 + Composition_PC3, data = pool, method = "nearest", distance = "mahalanobis",
             caliper = c(Composition_PC1 = 1, Composition_PC3 = 1), std.caliper = TRUE, replace = FALSE, ratio = 1)
md <- match.data(m)
coh$matched <- rownames(md)
capture.output(summary(m), file = file.path(OUT, "matchit_summary.txt"))

# ---- C. technical common support
coh$tech_profile_B <- rownames(pool)[pool$tech_profile == "B_high_multimap"]

sel <- do.call(rbind, lapply(names(coh), function(k) data.frame(cohort = k, animal = pool[coh[[k]], "animal"], group = pool[coh[[k]], "group"],
                                                                 sex = pool[coh[[k]], "sex"], Composition_PC1 = round(pool[coh[[k]], "Composition_PC1"], 2),
                                                                 Composition_PC3 = round(pool[coh[[k]], "Composition_PC3"], 2), tech_profile = pool[coh[[k]], "tech_profile"])))
write.table(sel, file.path(TAB, "common_support_cohorts_animals.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Common-support ranges: PC1", round(r1, 2), " PC3", round(r3, 2), "\n")
print(table(sel$cohort, sel$group))

# ---- DE per cohort
out <- list(); camres <- list(); rows <- list()
for (k in names(coh)) {
  s <- coh[[k]]; d <- droplevels(pool[s, ])
  nwt <- sum(d$genotype == "WT"); n5 <- sum(d$genotype == "5xFAD")
  if (nwt < 2 || n5 < 2) { rows[[k]] <- data.frame(cohort = k, n_WT = nwt, n_5xFAD = n5, note = "too few animals for DE"); next }
  des <- if (all(table(d$sex, d$genotype) > 0) && nrow(d) > 4) ~ sex + genotype else ~ genotype
  fit <- run_deseq(obj$counts[, s], d, des)
  r <- results_tbl(fit, NULL, obj$genes, name = "genotype_5xFAD_vs_WT")
  write.table(r, file.path(TAB, paste0("DE_", k, "_5xFAD_vs_WT.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  out[[k]] <- r
  j <- merge(base[, c("gene_id", "log2FoldChange")], r[, c("gene_id", "log2FoldChange", "padj")], by = "gene_id", suffixes = c("_base", "_coh"))
  jb <- j[j$gene_id %in% bsig, ]
  rows[[k]] <- data.frame(cohort = k, n_WT = nwt, n_5xFAD = n5, design = deparse(des), padj05 = sum(r$padj < 0.05, na.rm = TRUE),
                          LFC_cor_baseline_sig = round(cor(jb$log2FoldChange_base, jb$log2FoldChange_coh), 2),
                          direction_concordance_baseline_sig = round(mean(sign(jb$log2FoldChange_base) == sign(jb$log2FoldChange_coh)), 3),
                          median_retention = round(median(jb$log2FoldChange_coh / jb$log2FoldChange_base), 2),
                          baseline_sig_also_sig = sum(jb$padj < 0.05, na.rm = TRUE), note = "")
  camres[[k]] <- camera_sets(obj$vst[, s], d, des, "genotype5xFAD", file.path(TAB, "gene_sets.tsv"))
  write.table(camres[[k]], file.path(TAB, paste0("pathways_camera_", k, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}
summ <- do.call(rbind, lapply(rows, function(x) { for (c in c("design", "padj05", "LFC_cor_baseline_sig", "direction_concordance_baseline_sig",
                                                               "median_retention", "baseline_sig_also_sig")) if (!c %in% names(x)) x[[c]] <- NA; x }))
write.table(summ, file.path(TAB, "common_support_DE_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ)

# ---- Figure 10
p <- meta; p$in_matched <- rownames(p) %in% coh$matched; p$in_restricted <- rownames(p) %in% coh$restricted
p <- p[p$group %in% c("WT_NTC", "5xFAD_NTC"), ]
p$status <- ifelse(p$in_matched, "matched", ifelse(p$in_restricted, "restricted only", "excluded"))
g1 <- ggplot(p, aes(Composition_PC1, Composition_PC3)) +
  annotate("rect", xmin = r1[1], xmax = r1[2], ymin = r3[1], ymax = r3[2], fill = "#cde2fb", alpha = 0.6) +
  geom_point(aes(colour = group, shape = status), size = 3.5, stroke = 1) + geom_text(aes(label = animal), size = 2.7, vjust = -1.2, colour = INK2) +
  scale_colour_manual(values = PALS$group) + scale_shape_manual(values = c(matched = 16, `restricted only` = 1, excluded = 4)) +
  labs(title = "Composition common support (shaded) and matched animals", subtitle = "WT_NTC vs 5xFAD_NTC; Mahalanobis 1:1, caliper 1 SD") + theme_ca()
ggsave(file.path(OUT, "Fig10a_common_support_selection.png"), g1, width = 7.5, height = 6, dpi = 150)
if (length(out)) {
  lf <- do.call(rbind, lapply(names(out), function(k) { j <- merge(base[base$gene_id %in% bsig, c("gene_id", "log2FoldChange")], out[[k]][, c("gene_id", "log2FoldChange")], by = "gene_id")
    data.frame(cohort = k, base = j$log2FoldChange.x, coh = j$log2FoldChange.y) }))
  g2 <- ggplot(lf, aes(base, coh)) + geom_hline(yintercept = 0, colour = INK2) + geom_vline(xintercept = 0, colour = INK2) +
    geom_point(size = 0.4, alpha = 0.4, colour = "#2a78d6") + geom_abline(slope = 1, linetype = 3) + facet_wrap(~cohort) +
    coord_cartesian(xlim = c(-4, 4), ylim = c(-4, 4)) + labs(x = "baseline (full cohort) log2FC", y = "cohort log2FC",
    title = "Do baseline genotype effects reproduce in common-support cohorts?") + theme_ca()
  ggsave(file.path(OUT, "Fig10b_common_support_LFC.png"), g2, width = 11, height = 4.5, dpi = 150)
}
writeLines(capture.output(sessionInfo()), file.path(CA, "logs", "sessionInfo_06.txt"))
