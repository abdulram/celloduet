#!/usr/bin/env Rscript
# 13b - SPP1 vs NTC pathway testing for the pathway-link map (5xFAD animals, n = 12).
# Gene sets: all 50 MSigDB Hallmark v7.0 (human symbols -> mouse by case-insensitive match), the curated modules of
# step 03b, and mechanistic modules added for the SPP1 hypotheses (literature-based, small; see MECH below).
# Test: limma camera (sample-level, competitive) on blind VST, reported twice: inter-gene correlation fixed at the
#   limma default 0.01 (liberal when set genes co-vary) and estimated per set (inter.gene.cor = NA, conservative;
#   used for claims), under
#   unadjusted  ~ sex + treatment
#   fibro_adj   ~ sex + fibroblast_stromal + treatment
# Leading genes per set and model: set members with DESeq2 p < 0.05 in the set's direction (used for link map).

suppressPackageStartupMessages({library(limma)})
set.seed(7)
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
CA <- normalizePath(file.path(here, "..")); ROOT <- dirname(CA); TAB <- file.path(CA, "tables")
OUT <- file.path(CA, "13_pathway_links"); dir.create(OUT, showWarnings = FALSE)
obj <- readRDS(file.path(CA, "01_QC", "qc_objects.rds"))
meta <- read.delim(file.path(TAB, "sample_metadata_with_composition.tsv"), row.names = 1, check.names = FALSE)
sig <- read.delim(file.path(TAB, "composition_signature_scores_all.tsv"), row.names = 1, check.names = FALSE)
genes <- obj$genes
f5 <- meta[meta$genotype == "5xFAD", ]; f5$treatment <- factor(f5$treatment, levels = c("NTC", "SPP1")); f5$sex <- factor(f5$sex)
f5$fibro <- as.numeric(scale(sig[rownames(f5), "fibroblast_stromal"]))
V <- obj$vst[, rownames(f5)]
sym <- genes[rownames(V), "gene_name"]

# ---- gene sets
sets <- list()
for (l in readLines(file.path(ROOT, "data", "genesets", "h.all.v7.0.symbols.gmt"))) {
  f <- strsplit(l, "\t")[[1]]; sets[[paste0("H_", sub("HALLMARK_", "", f[1]))]] <- f[-(1:2)]
}
sets <- lapply(sets, function(g) unique(sym[toupper(sym) %in% toupper(g)]))
gs <- read.delim(file.path(TAB, "gene_sets.tsv"))
for (s in unique(gs$set[grepl("^CURATED_|^TECHNICAL_", gs$set)])) sets[[sub("CURATED_", "C_", s)]] <- gs$gene_name[gs$set == s]
MECH <- list(
  M_SPP1_RECEPTORS_INTEGRIN_FAK = c("Spp1", "Cd44", "Itgav", "Itgb1", "Itgb3", "Itgb5", "Itga4", "Itga5", "Itga8", "Itga9", "Ptk2", "Src", "Ilk", "Pxn", "Tln1", "Vcl", "Bcar1", "Crk"),
  M_YAP_TAZ_CORE_TARGETS = c("Ccn1", "Ccn2", "Ankrd1", "Amotl2", "Axl", "Ajuba", "Lats2", "Crim1", "Cdc42ep3", "Arhgef17", "Tgfb2", "Fstl1", "Serpine1", "Nuak2", "Dock5", "Gadd45b", "Foxf2", "Wwc1", "Thbs1", "Tead1"),
  M_YAP_TAZ_PATHWAY = c("Yap1", "Wwtr1", "Tead1", "Tead2", "Tead3", "Tead4", "Lats1", "Lats2", "Stk3", "Stk4", "Sav1", "Mob1a", "Nf2", "Amot", "Amotl1", "Amotl2", "Vgll4"),
  M_MECHANOTRANSDUCTION_ACTOMYOSIN = c("Rhoa", "Rock1", "Rock2", "Myl9", "Myh9", "Myh10", "Acta2", "Tagln", "Cnn1", "Cnn2", "Zyx", "Lpp", "Flna", "Actn1", "Actn4", "Piezo1", "Trpv4"),
  M_TGFB_SMAD_TARGETS = c("Serpine1", "Smad7", "Tgfbi", "Pmepa1", "Skil", "Ccn2", "Fn1", "Col1a1", "Col3a1", "Ltbp2", "Postn", "Timp1", "Lox", "Junb", "Id1", "Bhlhe40"),
  M_TYPE_I_IFN_EXTENDED = c("Rsad2", "Isg20", "Ifit1", "Ifit2", "Ifit3", "Ifit3b", "Isg15", "Irf7", "Irf9", "Stat1", "Stat2", "Oasl1", "Oasl2", "Oas1a", "Mx1", "Mx2", "Bst2", "Usp18", "Rtp4", "Xaf1", "Ifi44", "Ifi27l2a", "Zbp1", "Ddx58", "Ifih1"),
  M_CHOLESTEROL_SREBP = c("Insig1", "Hmgcr", "Hmgcs1", "Ldlr", "Sqle", "Fdft1", "Idi1", "Cyp51", "Dhcr7", "Dhcr24", "Msmo1", "Sc5d", "Lss", "Mvd", "Fdps", "Srebf2", "Pcsk9", "Stard4"),
  M_EOSINOPHIL_GRANULE = c("Prg2", "Prg3", "Epx", "Ear1", "Ear2", "Ear6", "Ccr3", "Siglecf", "Il5ra", "Alox15"),
  M_ARACHNOID_BARRIER = c("Cldn11", "Cdh1", "Tjp1", "Ocln", "Slc22a6", "Abcg2", "Dpp4", "Crabp2", "Aldh1a2", "Slc47a1")
)
sets <- c(sets, MECH)
idx <- lapply(sets, function(g) which(sym %in% g)); idx <- idx[sapply(idx, length) >= 5]
long <- do.call(rbind, lapply(names(idx), function(s) data.frame(set = s, gene = sym[idx[[s]]], gene_id = rownames(V)[idx[[s]]])))
write.table(long, file.path(TAB, "spp1_pathway_sets_long.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

de <- list(unadjusted = read.delim(file.path(TAB, "DE_SPP1_vs_NTC_5xFADonly.tsv")),
           fibro_adj = read.delim(file.path(TAB, "DE_SPP1_vs_NTC_fibroblast_adjusted.tsv")))
designs <- list(unadjusted = model.matrix(~ sex + treatment, f5), fibro_adj = model.matrix(~ sex + fibro + treatment, f5))
res <- NULL
for (m in names(designs)) {
  X <- designs[[m]]
  cam <- camera(V, idx, design = X, contrast = which(colnames(X) == "treatmentSPP1"))
  camc <- camera(V, idx, design = X, contrast = which(colnames(X) == "treatmentSPP1"), inter.gene.cor = NA)[rownames(cam), ]
  fit <- eBayes(lmFit(V, X)); t <- fit$t[, "treatmentSPP1"]
  d <- de[[m]]; rownames(d) <- d$gene_id
  lead <- sapply(rownames(cam), function(s) {
    ii <- idx[[s]]; dir <- ifelse(cam[s, "Direction"] == "Up", 1, -1)
    g <- rownames(V)[ii]; ok <- !is.na(d[g, "pvalue"]) & d[g, "pvalue"] < 0.05 & sign(d[g, "log2FoldChange"]) == dir
    paste(sort(sym[ii][ok]), collapse = ";")
  })
  res <- rbind(res, data.frame(model = m, set = rownames(cam), n_genes = cam$NGenes, direction = cam$Direction,
                               mean_t = sapply(rownames(cam), function(s) mean(t[idx[[s]]])), p = cam$PValue, FDR = cam$FDR,
                               p_conservative = camc$PValue, FDR_conservative = camc$FDR,
                               leading_genes = lead, n_leading = sapply(strsplit(lead, ";"), length)))
}
write.table(res, file.path(TAB, "spp1_pathway_camera_both_models.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (m in names(designs)) { cat("\n==", m, "\n"); r <- res[res$model == m, ]; print(head(r[order(r$p_conservative), c("set", "direction", "mean_t", "p", "FDR", "p_conservative", "FDR_conservative", "n_leading")], 25), digits = 3) }
