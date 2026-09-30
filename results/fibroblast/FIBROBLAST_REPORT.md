# Fibroblast-focused re-analysis (minimizing contamination)

**Question:** can we restrict the analysis to fibroblast-specific genes so that bone, olfactory, choroid plexus, blood and immune carry-over stops driving the results?

**Short answer:** yes, partly. The analysis uses 443 genes whose bulk signal comes mostly from fibroblasts. It confirms two points:
- The **WT vs 5xFAD difference is still dominated by the technical/batch axis**, even within fibroblast genes.
- In SPP1 animals there is a **modest, sex-consistent reduction of a fibroblast matrix (ECM) program**, with Cldn11, Fbln1, Thbd and Loxl3 down. Part of this effect overlaps with bone (osteoblast) content, so it needs validation.

Code: `analysis/reference/build_reference.py` (single-cell reference) and `analysis/fibroblast_analysis.py`.

---

## Method

1. **Single-cell reference.** We built a reference from the CELLxGENE Census (release 2025-01-30, read directly from its public S3 bucket). It has 33 mouse cell classes, about 1,200 cells each, balanced across studies:
   - **Fibroblasts:** dura (1 study), leptomeningeal/arachnoid-pia (7 studies) and brain perivascular (4 studies).
   - **Everything else that could be in this tissue:** osteoblast, osteoclast, chondrocyte, choroid plexus epithelium, olfactory neurons/epithelium/ensheathing cells, erythrocytes, neutrophils, monocytes, macrophages, microglia, B/T/NK/plasma/mast/dendritic cells, endothelium, pericytes, smooth muscle, skeletal muscle, neurons, glia, Schwann cells, keratinocytes and adipocytes.
   - Osteoclast, choroid plexus and olfactory profiles exist in the Census only in an embryonic atlas, so those three come from embryonic cells.
2. **Deconvolution.** Each bulk sample was deconvolved against the reference using NNLS on 25 marker genes per class (`F2_deconvolution.png`).
3. **Fibroblast signal share per gene.** Using each sample's estimated composition, we calculated what fraction of each gene's predicted bulk signal comes from fibroblasts. This weights every contaminating cell type by how abundant it actually is in these samples.
   - **Strict:** median share ≥ 80% and ≥ 60% in every sample → **151 genes**.
   - **Relaxed:** median share ≥ 60% and ≥ 40% in every sample → **292 more genes**.
   - **Bulk safety net:** we dropped 785 more genes whose expression still tracks a contaminant module (|r| ≥ 0.5; mostly the muscle and olfactory modules). This catches cell types the reference lacks, such as calvarial periosteum.
   - We first tried the simpler rule "10× higher in fibroblasts than in *every* atlas cell type". Only 31 genes passed. Fibroblasts share most matrix genes (Col1a1, Dcn, Lum, Fmod) with chondrocytes, osteoblasts and pericytes, and that rule ignores how rare those cells are here.
4. **Fibroblast-normalized DESeq2.** DESeq2 was run on the 443 genes only, so the size factors come from fibroblast genes. Effects are therefore per unit of fibroblast content, and the non-fibroblast part of the transcriptome can't influence normalization. Two models: `~ sex + group` and `~ sex + pct_mito + group`.
5. **Enrichment.** Preranked GSEA within the fibroblast gene universe (Hallmark plus custom sets).

**The gene set does capture fibroblasts:** the 443 genes are dominated by the Hallmark EMT/ECM program (33 genes, FDR 2e-16). The most highly expressed strict genes include Gpx3, Gsn, Ogn, Serpinf1, Mfap4, Col3a1, Igfbp4, Lum, Fmod, Fbln1, Cxcl12, Smoc2 and Cldn11 (`F1_reference_specificity_heatmap.png`). By best-matching class, 296 are dura-type, 59 leptomeningeal and 88 perivascular.

## Findings

### Composition (deconvolution)
The tissue is **dura-dominated**: dura-type fibroblasts account for about 50% of marker signal.

* The remainder is mostly osteoblasts (~9%), macrophages (~9%) and endothelium (~5–10%).
* Flagged samples:
  * **0453:** 47% plasma cells and 22% olfactory epithelium.
  * **613:** 24% plasma cells and 10% olfactory epithelium. The plasma cells plus olfactory epithelium point to nasal mucosa carry-over.
  * **620:** 19% osteoclasts and 22% osteoblasts, which is why its Spp1 is so extreme.
* Fibroblast fraction does not differ between groups (5xFAD vs WT p = 0.82; SPP1 vs NTC p = 0.48). So a change in fibroblast abundance can't explain group differences.

### 5xFAD NTC vs WT NTC: still confounded
| Model (fibroblast genes only) | padj < 0.05 |
|---|---|
| `~ sex + group` | 198 / 443 |
| `~ sex + pct_mito + group` | **0** |

Restricting to fibroblast genes doesn't remove the batch-like axis described in the main report. It acts on every cell type's RNA, not just on composition. **The genotype comparison remains uninterpretable without batch-balanced samples.**

### 5xFAD SPP1 vs 5xFAD NTC: a modest fibroblast matrix signal
| Gene | log2FC (fibroblast-normalized) | padj `~sex+group` | padj `+pct_mito` |
|---|---|---|---|
| **Cldn11** (arachnoid barrier) | −0.37 | **0.011** | **0.009** |
| **Fbln1** | −0.35 | 0.20 | **0.021** |
| **Thbd** | −0.39 | 0.12 | **0.049** |
| **Loxl3** (collagen crosslinking) | −0.52 | 0.20 | **0.045** |
| **Ifit3b** (interferon) | +0.53 | 0.21 | **0.005** |

* Cldn11 was already the most consistent hit in the whole-tissue analysis. It holds in every model and in both sexes.
* **GSEA: the Hallmark EMT/ECM set is lower in SPP1** within fibroblast genes (NES −1.99, FDR 0.006; with mito adjustment NES −1.82, FDR 0.044). The leading edge is Fbln1, Col3a1, Mmp2, Lum, Ecm1, Col16a1, Timp3, Col4a1, Serpinh1 and others.

  So the matrix-down signal from the whole-tissue GSEA is **not simply bone or osteoclast contamination**: it persists inside fibroblast-specific genes. The interferon-up signal mostly doesn't carry over; interferon genes are largely non-fibroblast.
* **The effect is consistent across sexes** (`F7_spp1_fibroblast_ecm.png`). A score built from these genes is lower in SPP1 females and males. Within 5xFAD, it tracks bulk Spp1 (Spearman ρ = 0.83, p = 0.001; ρ = 0.77 without outlier 620).

**Robustness check** (score ~ SPP1 + sex + covariate, n = 12; `ecm_score_spp1_robustness.csv`):

| Adjustment | SPP1 effect | p |
|---|---|---|
| none | −0.23 | 0.024 |
| + olfactory contamination | −0.20 | 0.037 |
| + osteoclast fraction | −0.17 | 0.008 |
| excluding 620 | −0.17 | 0.032 |
| + osteoblast fraction | −0.10 | 0.11 |

The effect is robust to olfactory and osteoclast carry-over and to the outlier. It **roughly halves when adjusting for osteoblast content**, which has two readings:
- Residual confounding: skull/periosteal matrix cells share these genes with dura fibroblasts.
- A real biological link: Spp1 (osteopontin) is a matrix protein that regulates bone and mesenchymal cells, and the SPP1 group has slightly less osteoblast signal (median 8.7% vs 10.7%).

Bulk data can't tell these apart.

## Bottom line and next steps
* The **fibroblast-restricted view supports a real but modest SPP1 effect on meningeal fibroblast matrix and barrier genes** (Cldn11, Fbln1, Loxl3, collagens). It's hypothesis-level, with n = 6 per group and a partial overlap with bone content.
* It **does not rescue the genotype comparison**. That needs new batch-balanced WT samples.
* **Best validation routes:**
  1. Peel the dura cleanly off the calvaria, or FACS-sort PDGFRα+ or CD45−CD31−CD140a+ fibroblasts, then repeat on sorted cells.
  2. Stain dura whole-mounts for CLDN11 and FBLN1/collagen IV.
  3. Check Spp1 knockdown directly in the cell type your SPP1 construct targets.

## Caveats
* The single-cell reference mixes platforms (10x and sci-RNA-seq) and developmental stages. The osteoclast, choroid plexus and olfactory profiles are embryonic, and dura fibroblasts come from a single study. Signal shares are estimates.
* Normalizing on fibroblast genes shifts all fold changes by a per-contrast constant relative to whole-tissue normalization. The ranking of genes within a contrast is essentially unchanged; what changes is which genes are allowed into the test.
* Preranked GSEA permutes genes, not samples, so its q-values are optimistic at n = 6 vs 6.

## Files
| File | Contents |
|---|---|
| `fibroblast_specific_genes.csv` | The 443 genes: tier, fibroblast signal share, reference expression, dominant contaminant |
| `fibroblast_specificity_all_genes.csv.gz` | Same metrics for every gene |
| `deconvolution_proportions.csv`, `fibroblast_scores_per_sample.csv` | Per-sample composition |
| `DE_<contrast>_<fibnorm\|fibnorm_mitoadj>.csv`, `DE_summary.csv` | Fibroblast-normalized DESeq2 |
| `GSEA_*.csv` | GSEA within fibroblast genes |
| `ecm_score_*.csv` | SPP1 ECM score, correlation with Spp1, robustness models |
| `compare_fibnorm_vs_wholetissue_*.csv` | Effect sizes vs the whole-tissue analysis |
| `fibroblast_normalised_vst.csv` | Normalized expression for plotting |
| `F1`–`F7` figures | Reference specificity, deconvolution, technical axis, volcanoes, GSEA, heatmaps, ECM score |
| `../../data/reference/` | Census pseudobulk reference and the list of cells used |
