#!/usr/bin/env Rscript
# 01 - Load counts + metadata, QC, VST, transcriptome PCA.
#
# Inputs (read-only):
#   data/VT3YXB-expression-matrix.tsv     gene x sample counts + CPM (facility output)
#   data/VT3YXB-mapping-stats.csv         uniquely / multi / unmapped reads per sample
#   data/VT3YXB-gene-biotype-summary.csv  genes with >=5 reads per biotype
# Metadata: parsed from facility sample names "<animal>_<genotype>_<sex>_<treatment>".
#   Age, batch, dissection date/person and pathology are NOT available anywhere in the project.
#   Two proxies are recorded and treated as unverified:
#     id_series   animal-ID numbering (leading "0" = 04xx series vs 6xx series) - possible cohort
#     tech_profile library mapping profile (see below) - an observed technical axis
# Outputs: composition_confounding_analysis/01_QC, tables/, intermediate RDS objects.

suppressPackageStartupMessages({library(DESeq2); library(ggplot2); library(pheatmap)})
set.seed(7)
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
ROOT <- normalizePath(file.path(here, "..", ".."))
CA <- file.path(ROOT, "composition_confounding_analysis")
OUT <- file.path(CA, "01_QC"); TAB <- file.path(CA, "tables")
source(file.path(here, "common.R"))

# ------------------------------------------------------------------ load
raw <- read.delim(file.path(ROOT, "data", "VT3YXB-expression-matrix.tsv"), check.names = FALSE)
raw <- raw[!grepl("^ERCC", raw$gene_id), ]          # ERCC spike-ins: zero reads in every sample
lib_idx <- 1:18
libs <- paste0("VT3YXB_", lib_idx)
sample_names <- c("0453_5xFAD_M_SPP1", "612_5xFAD_M_SPP1", "613_5xFAD_M_SPP1", "621_5xFAD_F_SPP1", "0463_5xFAD_F_SPP1",
                  "0455_5xFAD_F_SPP1", "616_WT_M_NTC", "611_WT_M_NTC", "626_WT_M_NTC", "624_WT_F_NTC", "623_WT_F_NTC",
                  "0466_WT_F_NTC", "617_5xFAD_M_NTC", "0452_5xFAD_M_NTC", "0451_5xFAD_M_NTC", "629_5xFAD_F_NTC",
                  "0469_5xFAD_F_NTC", "620_5xFAD_F_NTC")
counts <- as.matrix(raw[, paste0(libs, "_count")]); storage.mode(counts) <- "integer"
rownames(counts) <- raw$gene_id; colnames(counts) <- libs
genes <- data.frame(gene_id = raw$gene_id, gene_name = ifelse(is.na(raw$gene_name) | raw$gene_name == "", raw$gene_id, raw$gene_name),
                    biotype = raw$gene_biotype, row.names = raw$gene_id)

p <- do.call(rbind, strsplit(sample_names, "_"))
meta <- data.frame(library = libs, sample = sample_names, animal = p[, 1], genotype = p[, 2], sex = p[, 3],
                   treatment = p[, 4], row.names = libs)
meta$group <- paste(meta$genotype, meta$treatment, sep = "_")
meta$id_series <- ifelse(substr(meta$animal, 1, 1) == "0", "04xx", "6xx")
meta$age <- NA; meta$batch <- NA; meta$dissection <- NA       # not provided

ms <- read.csv(file.path(ROOT, "data", "VT3YXB-mapping-stats.csv"), row.names = 1, check.names = FALSE)[sample_names, ]
meta$total_reads <- rowSums(ms)
meta$pct_unique <- 100 * ms[["Uniquely Mapped"]] / meta$total_reads
meta$pct_multi <- 100 * ms[["Multi-mapped"]] / meta$total_reads
meta$assigned_counts <- colSums(counts)
meta$genes_detected <- colSums(counts >= 5)
mt <- grepl("^mt-", genes$gene_name)
meta$pct_mito <- 100 * colSums(counts[mt, ]) / colSums(counts)
hb <- genes$gene_name %in% c("Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2")
meta$pct_hemoglobin <- 100 * colSums(counts[hb, ]) / colSums(counts)
xist <- counts[genes$gene_name == "Xist", ]; ychr <- colSums(counts[genes$gene_name %in% c("Ddx3y", "Uty", "Kdm5d", "Eif2s3y"), ])
meta$sex_check <- ifelse(xist > ychr, "F", "M"); meta$sex_ok <- meta$sex_check == meta$sex
# technical profile: two-cluster split on % multi-mapped (bimodal: ~10% vs ~16-20%); recorded, not assumed biological
meta$tech_profile <- ifelse(meta$pct_multi < 13, "A_low_multimap", "B_high_multimap")

# ------------------------------------------------------------------ filter + VST
small <- genes$biotype %in% c("rRNA", "Mt_rRNA", "Mt_tRNA", "misc_RNA", "snRNA", "snoRNA", "scaRNA", "miRNA", "ribozyme")
keep <- rowSums(counts >= 10) >= 3 & !small
cf <- counts[keep, ]
cat("genes kept:", nrow(cf), "of", nrow(counts), "\n")
meta$group <- factor(meta$group, levels = c("WT_NTC", "5xFAD_NTC", "5xFAD_SPP1"))
meta$sex <- factor(meta$sex)
dds <- DESeqDataSetFromMatrix(cf, meta, design = ~ sex + group)
dds <- estimateSizeFactors(dds)
meta$size_factor <- sizeFactors(dds)
vsd <- vst(dds, blind = TRUE)
V <- assay(vsd)

# ------------------------------------------------------------------ PCA (top 2000 variable genes)
top <- order(apply(V, 1, var), decreasing = TRUE)[1:2000]
pc <- prcomp(t(V[top, ]), center = TRUE, scale. = FALSE)
ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
for (i in 1:5) meta[[paste0("Tx_PC", i)]] <- pc$x[, i]
load <- data.frame(gene = genes[rownames(V)[top], "gene_name"], pc$rotation[, 1:5])
write.table(load, file.path(OUT, "transcriptome_PCA_loadings.tsv"), sep = "\t", quote = FALSE)
write.table(data.frame(PC = paste0("Tx_PC", 1:10), pct_variance = round(ve[1:10], 2)),
            file.path(TAB, "transcriptome_PC_variance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
top_load <- sapply(1:5, function(i) paste(head(load$gene[order(-abs(load[[paste0("PC", i)]]))], 8), collapse = ","))
write.table(data.frame(PC = paste0("Tx_PC", 1:5), pct_variance = round(ve[1:5], 1), top_loading_genes = top_load),
            file.path(OUT, "transcriptome_PC_top_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------------ plots
lab <- function(i) sprintf("Tx_PC%d (%.1f%%)", i, ve[i])
pca_plot <- function(col, title, discrete = TRUE, a = 1, b = 2) {
  d <- meta; d$x <- d[[paste0("Tx_PC", a)]]; d$y <- d[[paste0("Tx_PC", b)]]
  g <- ggplot(d, aes(x, y)) + geom_point(aes(colour = .data[[col]], shape = sex), size = 3.2, stroke = 0.8) +
    geom_text(aes(label = animal), size = 2.7, vjust = -1.1, colour = INK2) +
    labs(x = lab(a), y = lab(b), title = title, colour = col) + theme_ca()
  if (discrete && col %in% names(PALS)) g <- g + scale_colour_manual(values = PALS[[col]])
  if (!discrete) g <- g + scale_colour_gradient(low = "#cde2fb", high = "#0d366b")
  g
}
ggsave(file.path(OUT, "PCA_by_genotype.png"), pca_plot("genotype", "Transcriptome PCA - genotype"), width = 7, height = 5.5, dpi = 150)
ggsave(file.path(OUT, "PCA_by_group.png"), pca_plot("group", "Transcriptome PCA - genotype x treatment"), width = 7, height = 5.5, dpi = 150)
ggsave(file.path(OUT, "PCA_by_sex.png"), pca_plot("sex", "Transcriptome PCA - sex"), width = 7, height = 5.5, dpi = 150)
ggsave(file.path(OUT, "PCA_by_id_series.png"), pca_plot("id_series", "Transcriptome PCA - animal-ID series (possible cohort; unverified)"), width = 7, height = 5.5, dpi = 150)
ggsave(file.path(OUT, "PCA_by_tech_profile.png"), pca_plot("tech_profile", "Transcriptome PCA - library mapping profile"), width = 7, height = 5.5, dpi = 150)
ggsave(file.path(OUT, "PCA_by_pct_mito.png"), pca_plot("pct_mito", "Transcriptome PCA - % mitochondrial reads", FALSE), width = 7, height = 5.5, dpi = 150)
ggsave(file.path(OUT, "PCA_PC3_PC4_by_group.png"), pca_plot("group", "Transcriptome PCA (PC3 vs PC4) - group", TRUE, 3, 4), width = 7, height = 5.5, dpi = 150)
png(file.path(OUT, "PCA_age_batch_NOT_AVAILABLE.png"), width = 700, height = 200); par(mar = c(0, 0, 0, 0)); plot.new()
text(0.5, 0.5, "Age, batch and dissection metadata are not available for this dataset;\nPCA by age/batch cannot be produced. See PCA_by_id_series / PCA_by_tech_profile.", cex = 1.1); dev.off()

ord <- order(meta$group, meta$sex)
qc_long <- rbind(
  data.frame(metric = "total reads (M)", value = meta$total_reads / 1e6, animal = meta$animal, group = meta$group),
  data.frame(metric = "assigned counts (M)", value = meta$assigned_counts / 1e6, animal = meta$animal, group = meta$group),
  data.frame(metric = "% uniquely mapped", value = meta$pct_unique, animal = meta$animal, group = meta$group),
  data.frame(metric = "% multi-mapped", value = meta$pct_multi, animal = meta$animal, group = meta$group),
  data.frame(metric = "% mitochondrial", value = meta$pct_mito, animal = meta$animal, group = meta$group),
  data.frame(metric = "genes with >=5 reads", value = meta$genes_detected, animal = meta$animal, group = meta$group))
qc_long$animal <- factor(qc_long$animal, levels = meta$animal[ord])
g <- ggplot(qc_long, aes(animal, value, fill = group)) + geom_col(width = 0.8, colour = SURFACE, linewidth = 0.6) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) + scale_fill_manual(values = PALS$group) +
  labs(x = NULL, y = NULL, title = "Library QC per animal") + theme_ca() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7))
ggsave(file.path(OUT, "library_QC.png"), g, width = 12, height = 6.5, dpi = 150)

cm <- cor(V[top, ]); colnames(cm) <- rownames(cm) <- paste(meta$animal, meta$group, meta$sex)
ann <- data.frame(group = meta$group, sex = meta$sex, tech_profile = meta$tech_profile, row.names = rownames(cm))
pheatmap(cm, annotation_col = ann, annotation_row = ann, annotation_colors = list(group = PALS$group, sex = PALS$sex, tech_profile = PALS$tech_profile),
         color = colorRampPalette(c("#f0efec", "#86b6ef", "#1c5cab", "#0d366b"))(100), fontsize = 7,
         main = "Sample-sample Pearson correlation (VST, top 2000 genes)", filename = file.path(OUT, "sample_correlation_heatmap.png"), width = 10, height = 9)

# ------------------------------------------------------------------ save
write.table(meta[, setdiff(names(meta), "library")], file.path(TAB, "sample_metadata_qc_pcs.tsv"), sep = "\t", quote = FALSE, col.names = NA)
saveRDS(list(counts = cf, genes = genes, meta = meta, vst = V, pca = pc, pct_var = ve), file.path(CA, "01_QC", "qc_objects.rds"))
write.table(V, gzfile(file.path(CA, "01_QC", "vst_blind.tsv.gz")), sep = "\t", quote = FALSE, col.names = NA)
write.table(cf, gzfile(file.path(CA, "01_QC", "filtered_counts.tsv.gz")), sep = "\t", quote = FALSE, col.names = NA)
cat("Tx PC variance %:", round(ve[1:5], 1), "\n")
print(table(meta$group, meta$tech_profile))
writeLines(capture.output(sessionInfo()), file.path(CA, "logs", "sessionInfo_01.txt"))
