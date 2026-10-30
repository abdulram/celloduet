#!/usr/bin/env Rscript
# 14b - CD11b+ microglia, SPP1 vs NTC (n = 3 vs 3, all male 5xFAD): DE and pathway tests mirroring dura step 13.
# DESeq2 models
#   unadjusted   ~ treatment
#   carry_adj    ~ carry_astro + treatment  (log astrocyte-marker share; near-collinear with treatment, sensitivity only)
# Gene-set tests: limma camera on CD11b VST, default inter-gene correlation (0.01, liberal) and estimated per set
# (inter.gene.cor = NA, conservative; used for claims), in two gene universes
#   all_genes     every tested gene
#   mg_intrinsic  only microglia-intrinsic genes (step 14a): carry-over cannot shift them > ~5%
# Gene sets: dura step-13 sets (Hallmark, dura curated, mechanistic) + CD11b modules (data/modules.gmt).
# Leading genes = members with unadjusted DESeq2 p < 0.05 in the set's direction.

suppressPackageStartupMessages({library(DESeq2); library(limma)})
set.seed(7)
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
root <- normalizePath(file.path(here, "..")); TAB <- file.path(root, "tables")
dura <- file.path(dirname(root), "composition_confounding_analysis", "tables")
OUT <- file.path(root, "14_pathway_links"); dir.create(OUT, showWarnings = FALSE)

m <- read.delim(gzfile(file.path(root, "data", "55ZS67-expression-matrix.tsv.gz")), check.names = FALSE)
m <- m[!grepl("^ERCC", m$gene_id), ]
cnt <- round(as.matrix(m[, grepl("_count$", colnames(m))])); rownames(cnt) <- m$gene_id
colnames(cnt) <- sub("_count$", "", colnames(cnt)); sym_all <- setNames(m$gene_name, m$gene_id)
qc <- read.delim(file.path(root, "01_QC", "qc_metrics.tsv"), row.names = 1)
pur <- read.delim(file.path(root, "02_purity_celltyping", "purity_carryover_per_sample.tsv"), row.names = 1)
s <- rownames(qc)[qc$fraction == "CD11b"]
cd <- data.frame(treatment = factor(qc[s, "treatment"], levels = c("NTC", "Spp1")),
                 carry_astro = as.numeric(scale(log(pur[s, "astro_marker_pct_vs_ACSA2"]))), row.names = s)
x <- cnt[, s]; x <- x[rowSums(x >= 10) >= 3, ]
cat("carry_astro vs treatment r =", round(cor(cd$carry_astro, as.numeric(cd$treatment)), 2), "\n")

de <- list()
for (mname in c("unadjusted", "carry_adj")) {
  f <- if (mname == "unadjusted") ~ treatment else ~ carry_astro + treatment
  dds <- DESeq(DESeqDataSetFromMatrix(x, cd, f), quiet = TRUE)
  r <- results(dds, name = "treatment_Spp1_vs_NTC")
  de[[mname]] <- data.frame(gene_id = rownames(r), gene_name = sym_all[rownames(r)], baseMean = r$baseMean,
                            log2FoldChange = r$log2FoldChange, lfcSE = r$lfcSE, stat = r$stat, pvalue = r$pvalue, padj = r$padj)
  if (mname == "unadjusted") V <- assay(vst(dds, blind = TRUE))
}
intr <- read.delim(file.path(TAB, "CD11b_microglia_intrinsic_genes.tsv"), row.names = 1)
intr$microglia_intrinsic <- as.logical(intr$microglia_intrinsic)
for (mname in names(de)) {
  d <- de[[mname]]; d$microglia_intrinsic <- intr[d$gene_id, "microglia_intrinsic"]
  d$reason_excluded <- intr[d$gene_id, "reason_excluded"]; d <- d[order(d$pvalue), ]
  write.table(d, file.path(TAB, paste0("DE_CD11b_SPP1_vs_NTC_", mname, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(mname, ": padj<0.05 =", sum(d$padj < 0.05, na.rm = TRUE), "; intrinsic padj<0.05 =",
      sum(d$padj < 0.05 & d$microglia_intrinsic, na.rm = TRUE), "\n")
  de[[mname]] <- d
}

# ---- gene sets (same definitions as dura 13b, plus CD11b modules)
sym <- sym_all[rownames(V)]
sets <- list()
for (l in readLines(file.path(root, "data", "h.all.v7.0.symbols.gmt"))) {
  f <- strsplit(l, "\t")[[1]]; sets[[paste0("H_", sub("HALLMARK_", "", f[1]))]] <- f[-(1:2)]
}
sets <- lapply(sets, function(g) unique(sym[toupper(sym) %in% toupper(g)]))
gs <- read.delim(file.path(dura, "gene_sets.tsv"))
for (st in unique(gs$set[grepl("^CURATED_", gs$set)])) sets[[sub("CURATED_", "C_", st)]] <- gs$gene_name[gs$set == st]
dl <- read.delim(file.path(dura, "spp1_pathway_sets_long.tsv"))
for (st in unique(dl$set[grepl("^M_", dl$set)])) sets[[st]] <- unique(dl$gene[dl$set == st])
for (l in readLines(file.path(root, "data", "modules.gmt"))) {
  f <- strsplit(l, "\t")[[1]]; if (grepl("^L_", f[1])) next  # lineage panels are purity markers, not pathways
  sets[[paste0("G_", f[1])]] <- f[-(1:2)]
}
intr_ids <- rownames(intr)[intr$microglia_intrinsic]
universes <- list(all_genes = rownames(V), mg_intrinsic = intersect(rownames(V), intr_ids))
X <- model.matrix(~ treatment, cd); ctr <- which(colnames(X) == "treatmentSpp1")
du <- de$unadjusted; rownames(du) <- du$gene_id
res <- NULL; long <- NULL
for (u in names(universes)) {
  Vu <- V[universes[[u]], ]; su <- sym[rownames(Vu)]
  idx <- lapply(sets, function(g) which(su %in% g)); idx <- idx[sapply(idx, length) >= 5]
  long <- rbind(long, do.call(rbind, lapply(names(idx), function(st) data.frame(universe = u, set = st, gene = su[idx[[st]]], gene_id = rownames(Vu)[idx[[st]]]))))
  cam <- camera(Vu, idx, design = X, contrast = ctr)
  camc <- camera(Vu, idx, design = X, contrast = ctr, inter.gene.cor = NA)[rownames(cam), ]
  t <- eBayes(lmFit(Vu, X))$t[, ctr]
  lead <- sapply(rownames(cam), function(st) {
    g <- rownames(Vu)[idx[[st]]]; dir <- ifelse(cam[st, "Direction"] == "Up", 1, -1)
    ok <- !is.na(du[g, "pvalue"]) & du[g, "pvalue"] < 0.05 & sign(du[g, "log2FoldChange"]) == dir
    paste(sort(su[idx[[st]]][ok]), collapse = ";")
  })
  res <- rbind(res, data.frame(model = u, set = rownames(cam), n_genes = cam$NGenes, direction = cam$Direction,
                               mean_t = sapply(rownames(cam), function(st) mean(t[idx[[st]]])), p = cam$PValue, FDR = cam$FDR,
                               p_conservative = camc$PValue, FDR_conservative = camc$FDR,
                               leading_genes = lead, n_leading = sapply(strsplit(lead, ";"), length)))
}
write.table(long, file.path(TAB, "CD11b_pathway_sets_long.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res, file.path(TAB, "CD11b_pathway_camera_both_universes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (u in names(universes)) { cat("\n==", u, "\n"); r <- res[res$model == u, ]
  print(head(r[order(r$p_conservative), c("set", "n_genes", "direction", "p", "FDR", "p_conservative", "FDR_conservative", "n_leading")], 22), digits = 3) }
