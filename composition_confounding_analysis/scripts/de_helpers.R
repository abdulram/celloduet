# Shared DE helpers (sourced by steps 04-06).
suppressPackageStartupMessages({library(DESeq2); library(limma); library(ggplot2); library(pheatmap)})

run_deseq <- function(counts, meta, design) {
  dds <- DESeqDataSetFromMatrix(counts, meta, design = design)
  stopifnot(qr(model.matrix(design, meta))$rank == ncol(model.matrix(design, meta)))   # full-rank check
  DESeq(dds, quiet = TRUE)
}

results_tbl <- function(dds, contrast, genes, name = NULL) {
  r <- if (is.null(name)) results(dds, contrast = contrast, alpha = 0.05) else results(dds, name = name, alpha = 0.05)
  d <- as.data.frame(r); d$gene_id <- rownames(d); d$gene_name <- genes[d$gene_id, "gene_name"]
  d[order(d$pvalue), c("gene_id", "gene_name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
}

ma_volcano <- function(res, out, tag, title) {
  d <- res[!is.na(res$padj), ]
  d$sig <- ifelse(d$padj < 0.05, ifelse(d$log2FoldChange > 0, "up", "down"), "ns")
  cols <- c(up = "#e34948", down = "#2a78d6", ns = "#c9c8c3")
  lab <- head(d[d$sig != "ns", ], 15)
  g1 <- ggplot(d, aes(log10(baseMean + 1), log2FoldChange, colour = sig)) + geom_point(size = 0.6, alpha = 0.7) +
    scale_colour_manual(values = cols) + geom_hline(yintercept = 0, colour = INK2) + labs(title = paste("MA:", title), colour = NULL) + theme_ca()
  g2 <- ggplot(d, aes(log2FoldChange, -log10(pvalue), colour = sig)) + geom_point(size = 0.7, alpha = 0.7) +
    scale_colour_manual(values = cols, labels = c(down = sprintf("down (%d)", sum(d$sig == "down")), ns = "ns", up = sprintf("up (%d)", sum(d$sig == "up")))) +
    geom_text(data = lab, aes(label = gene_name), size = 2.6, colour = INK, vjust = -0.6, check_overlap = TRUE) +
    labs(title = paste("Volcano:", title), colour = "padj<0.05") + theme_ca()
  ggsave(file.path(out, paste0("MA_", tag, ".png")), g1, width = 7, height = 5, dpi = 150)
  ggsave(file.path(out, paste0("volcano_", tag, ".png")), g2, width = 7, height = 5.5, dpi = 150)
}

top_heatmap <- function(vst, res, meta, out, tag, genes, n = 40) {
  ids <- head(res$gene_id[!is.na(res$padj) & res$padj < 0.05], n)
  if (length(ids) < 2) return(invisible(NULL))
  ord <- rownames(meta)[order(meta$group, meta$sex)]
  m <- vst[ids, ord]; m <- t(scale(t(m))); rownames(m) <- genes[ids, "gene_name"]
  colnames(m) <- paste(meta[ord, "animal"], meta[ord, "group"], meta[ord, "sex"])
  ann <- data.frame(group = meta[ord, "group"], tech_profile = meta[ord, "tech_profile"],
                    Composition_PC1 = meta[ord, "Composition_PC1"], row.names = colnames(m))
  pheatmap(m, cluster_cols = FALSE, annotation_col = ann, annotation_colors = list(group = PALS$group, tech_profile = PALS$tech_profile),
           color = colorRampPalette(c("#184f95", "#86b6ef", "#f0efec", "#f0a3a2", "#a8322f"))(100), breaks = seq(-2.5, 2.5, length.out = 101),
           fontsize_row = 7, fontsize_col = 7, main = paste("Top genes:", tag), filename = file.path(out, paste0("heatmap_top_", tag, ".png")), width = 9, height = 9)
}

camera_sets <- function(vst, meta, design, coef, gs_file, keep = colnames(vst)) {
  gs <- read.delim(gs_file)
  idx <- lapply(split(gs$gene_id, gs$set), function(g) which(rownames(vst) %in% g))
  idx <- idx[sapply(idx, length) >= 5]
  X <- model.matrix(design, meta[keep, , drop = FALSE])
  cam <- camera(vst[, keep], idx, design = X, contrast = which(colnames(X) == coef), inter.gene.cor = NA)
  cam$set <- rownames(cam); cam$n_genes <- cam$NGenes
  # effect size: mean moderated-t of the set genes (direction + magnitude)
  fit <- eBayes(lmFit(vst[, keep], X))
  t <- fit$t[, coef]
  cam$mean_t <- sapply(cam$set, function(s) mean(t[idx[[s]]]))
  cam[, c("set", "n_genes", "Direction", "mean_t", "PValue", "FDR")]
}
