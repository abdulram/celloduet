# 04: DESeq2 per fraction, Spp1 vs NTC (~ treatment), plus sensitivity models,
# leave-one-out robustness, VST/PCA input, and limma-voom camera on curated modules + Hallmark.
suppressPackageStartupMessages({library(DESeq2); library(limma); library(edgeR)})
set.seed(1)
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE))), ".."))
m <- read.delim(gzfile(file.path(root, "data", "55ZS67-expression-matrix.tsv.gz")), check.names=FALSE)
m <- m[!grepl("^ERCC", m$gene_id), ]
cnt <- as.matrix(m[, grepl("_count$", colnames(m))]); rownames(cnt) <- m$gene_id
colnames(cnt) <- sub("_count$", "", colnames(cnt))
cnt <- round(cnt)  # facility splits multimappers fractionally; DESeq2 needs integers
ann <- m[, 1:3]; rownames(ann) <- ann$gene_id
qc <- read.delim(file.path(root, "01_QC", "qc_metrics.tsv"), row.names=1)
pur <- read.delim(file.path(root, "02_purity_celltyping", "purity_carryover_per_sample.tsv"), row.names=1)
gmt <- strsplit(readLines(file.path(root, "data", "modules.gmt")), "\t")
mods <- setNames(lapply(gmt, function(x) x[-(1:2)]), sapply(gmt, `[`, 1))
hg <- strsplit(readLines(file.path(root, "data", "h.all.v7.0.symbols.gmt")), "\t")
hall <- setNames(lapply(hg, function(x) x[-(1:2)]), sapply(hg, `[`, 1))
up1 <- function(x) paste0(substr(x,1,1), tolower(substring(x,2)))
hall <- lapply(hall, up1)  # human->mouse by case-insensitive symbol match

run_frac <- function(frac, outdir) {
  s <- rownames(qc)[qc$fraction == frac]
  cd <- data.frame(treatment=factor(qc[s, "treatment"], levels=c("NTC","Spp1")),
                   assigned=scale(qc[s, "pct_assigned_of_unique"])[,1],
                   unmapped=scale(qc[s, "pct_unmapped"])[,1],
                   carry_mg=scale(log(pur[s, "microglial_RNA_pct_vs_CD11b"]))[,1],
                   carry_astro=scale(log(pur[s, "astro_marker_pct_vs_ACSA2"]))[,1],
                   oligo=scale(log(pur[s, "Plp1_over_Atp1b2"]))[,1],
                   ieg=scale(pur[s, "dissociation_IEG_log2cpm_mean"])[,1], row.names=s)
  x <- cnt[, s]; x <- x[rowSums(x >= 10) >= 3, ]
  dds <- DESeqDataSetFromMatrix(x, cd, ~ treatment); dds <- DESeq(dds, quiet=TRUE)
  r <- results(dds, name="treatment_Spp1_vs_NTC")
  sh <- lfcShrink(dds, coef="treatment_Spp1_vs_NTC", type="normal", quiet=TRUE)
  out <- data.frame(gene_id=rownames(r), gene=ann[rownames(r), "gene_name"], baseMean=r$baseMean,
                    log2FC=r$log2FoldChange, lfcSE=r$lfcSE, log2FC_shrunk=sh$log2FoldChange,
                    stat=r$stat, pvalue=r$pvalue, padj=r$padj)
  # sensitivity: one covariate at a time (5 residual df -> 3 df); covariates are near-collinear with treatment in places
  cov <- if (frac == "CD11b") c("assigned", "carry_astro", "ieg") else c("unmapped", "carry_mg", "oligo", "ieg")
  for (cv in cov) {
    cat(frac, cv, "cor with treatment:", round(cor(cd[[cv]], as.numeric(cd$treatment)), 2), "\n")
    d2 <- DESeqDataSetFromMatrix(x, cd, as.formula(paste("~", cv, "+ treatment")))
    d2 <- DESeq(d2, quiet=TRUE); r2 <- results(d2, name="treatment_Spp1_vs_NTC")
    out[[paste0("log2FC_adj_", cv)]] <- r2$log2FoldChange; out[[paste0("padj_adj_", cv)]] <- r2$padj
  }
  # leave-one-out: max padj and min |LFC| across the 6 drop-one fits
  loo_p <- loo_l <- matrix(NA, nrow(x), length(s), dimnames=list(rownames(x), s))
  for (d in s) {
    k <- setdiff(s, d); d3 <- DESeqDataSetFromMatrix(x[, k], cd[k, ], ~ treatment)
    d3 <- DESeq(d3, quiet=TRUE); r3 <- results(d3, name="treatment_Spp1_vs_NTC")
    loo_p[, d] <- r3$pvalue; loo_l[, d] <- r3$log2FoldChange
  }
  out$loo_max_pvalue <- apply(loo_p, 1, max, na.rm=TRUE)
  out$loo_same_sign <- apply(sign(loo_l) == sign(out$log2FC), 1, all)
  out <- out[order(out$pvalue), ]
  write.table(out, file.path(root, outdir, paste0("DE_Spp1_vs_NTC_", frac, ".tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  vs <- assay(vst(dds, blind=TRUE))
  write.table(vs, gzfile(file.path(root, "04_exploration", paste0("vst_", frac, ".tsv.gz"))), sep="\t", quote=FALSE)
  # camera (limma-voom), inter-gene correlation estimated
  y <- DGEList(x); y <- calcNormFactors(y); des <- model.matrix(~ treatment, cd)
  v <- voom(y, des); sym <- ann[rownames(v), "gene_name"]
  cam <- function(sets, tag) {
    idx <- ids2indices(sets, sym); idx <- idx[sapply(idx, length) >= 3]
    cr <- camera(v, idx, des, contrast=2); cr$set <- rownames(cr); cr$collection <- tag
    cr$n_tested <- sapply(idx[rownames(cr)], length); cr }
  cr <- rbind(cam(mods, "curated"), cam(hall, "Hallmark"))
  cr <- cr[order(cr$PValue), c("collection","set","n_tested","Direction","PValue","FDR")]
  write.table(cr, file.path(root, "07_modules", paste0("camera_Spp1_vs_NTC_", frac, ".tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  cat("\n==", frac, ": genes tested", nrow(out), "| padj<0.05:", sum(out$padj < 0.05, na.rm=TRUE),
      "| padj<0.1:", sum(out$padj < 0.1, na.rm=TRUE), "\n")
  print(head(out[, c("gene","baseMean","log2FC","padj","loo_max_pvalue","loo_same_sign")], 30), row.names=FALSE)
  print(head(cr, 15), row.names=FALSE)
}
run_frac("CD11b", "05_DE_CD11b")
run_frac("ACSA2", "06_DE_ACSA2")
writeLines(capture.output(sessionInfo()), file.path(root, "logs", "sessionInfo_R.txt"))
