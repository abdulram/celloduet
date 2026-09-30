# Fibroblast-focused re-analysis (minimizing contamination)

**Question:** can we restrict the analysis to fibroblast-specific genes so that bone, olfactory, blood and immune carry-over stops driving the results?

**Short answer.** Yes. There are 151 genes whose bulk signal is predicted to come mostly from fibroblasts (dominated by the matrix/EMT program, FDR 6e-10). Restricting the analysis to them shows:
- The **WT vs 5xFAD difference is still the technical/batch axis**.
- In SPP1 animals, a **small set of fibroblast genes shifts: Cldn11, Thbd and Fbln1 down; Matn4 and Ifit3b up**. This holds across models and across two versions of the gene set.
- A **broader "matrix program down" claim is not robust**. It was significant with an earlier deconvolution, which had a flaw (see *Revision note*), and falls to p ≈ 0.1 after the fix.

Code: `analysis/reference/build_reference.py` (single-cell reference) and `analysis/fibroblast_analysis.py`.

---

## Deconvolution: how the cell representation was estimated

1. **Reference.** We built pseudobulk expression profiles for 33 mouse cell classes from the CELLxGENE Census (release 2025-01-30). Each class uses about 1,200 cells sampled across studies; each cell is normalized to counts per 10k, then averaged per study and then across studies.
   - **Fibroblasts:** dura, leptomeningeal and brain perivascular.
   - **Other classes:** everything plausibly in the tissue, including bone, osteoclast, cartilage, choroid plexus, olfactory, blood, immune, vascular, muscle and neural/glial cells.
   - Osteoclast, choroid plexus and olfactory profiles exist only in an embryonic atlas.
2. **Marker genes.** For each class, the marker genes are those expressed at ≥ 1 CP10K and ≥ 3–5× higher than in every other class (top 25–50). We excluded mitochondrial and ribosomal genes, because they follow the technical axis, and antibody/TCR variable genes, because they are clone-specific.
   - The dura fibroblast class absorbs all fibroblast signal. The leptomeningeal and perivascular fibroblast classes have too few genes that separate them from dura to be fitted separately.
3. **Fit.** Each bulk sample is fitted as a non-negative mixture of the class profiles (NNLS on the marker genes, weighted 1/√expression).
4. **Trimming.** Markers that the fit misses by more than 3–4× in a sample are dropped and the fit repeated (4–5 rounds). This removes genes from cell types the reference lacks. For example, nasal gland genes (Scgb1c1 at 5,400 CPM, Agr2) had been pulled into "plasma cells".
5. **Ensemble.** The final estimate averages 4 settings, because single settings proved sensitive. It agreed best with independent canonical-marker scores:

   | Class | Spearman ρ with marker score |
   |---|---|
   | Osteoclast | 0.88 |
   | Erythrocyte | 0.93 |
   | Neutrophil | 0.75 |
   | Muscle | 0.75 |
   | Olfactory | 0.55 |
   | Plasma | 0.57 |
   | Macrophage | 0.39 |
   | Osteoblast | 0.31 |
   | Choroid plexus | ≈ 0 |

   Marker-level fit quality is ρ ≈ 0.88–0.92 in every sample.

**How to read it.** Values are shares of *marker-gene signal*, not cell counts, because cell types differ in RNA content. Class-level estimates are approximate. The "Neural / CP / ependymal" slice (~15%) is the least trustworthy: canonical brain genes are low in these samples, so it probably absorbs unmodeled cell types.

### Results (`F2b_cell_representation_per_sample.png`, `F2_deconvolution.png`)
| Median % of signal | WT NTC | 5xFAD NTC | 5xFAD SPP1 |
|---|---|---|---|
| Fibroblast | 32 | 26 | 24 |
| Osteoblast / chondrocyte | 14 | 15 | 11 |
| Osteoclast | 4 | 7 | 5 |
| Erythrocyte | 4 | 7 | 13 |
| Immune | 15 | 12 | 12 |
| Vascular | 8 | 7 | 5 |

Per-sample standouts, all confirmed by raw reads:
* **620:** 46% osteoclast. The next-highest are **0455** (26%) and **0452** (23%). These are exactly the three samples with the highest Spp1, so bulk Spp1 mostly reports osteoclast content.
* **623, 616, 0455, 0452:** 20–32% erythrocyte. These are the samples with the most hemoglobin reads (7–12%).
* **0453:** 28% olfactory/nasal. **613:** about 7% olfactory/nasal.
* **Fibroblast fraction:** somewhat lower in 5xFAD than WT (p = 0.09), again aligned with the technical axis. It does not differ between SPP1 and NTC (p = 0.70).

## Selecting fibroblast-specific genes
Using each sample's estimated composition, we calculated for every gene the fraction of its predicted bulk signal that comes from fibroblasts. This weights each contaminating cell type by how abundant it actually is here. A simpler rule, "10× higher than every atlas cell type", keeps only about 30 genes, because fibroblasts share matrix genes with chondrocytes, osteoblasts and pericytes.

* **Strict:** median share ≥ 80% and ≥ 60% in every sample → **66 genes**.
* **Relaxed:** median share ≥ 60% and ≥ 40% in every sample → **85 genes**.
* **Removed:** genes whose bulk expression still tracks a contaminant module (|r| ≥ 0.5), and canonical bone, osteoclast, choroid plexus, olfactory, blood, muscle and brain markers (blacklisted).
* **What remains:** Ogn, Serpinf1, Mfap4, Igfbp4, Lum, Fmod, Fbln1, Cpxm1/2, Aebp1, Col6a2, Col14a1, Smoc2, Pi16, Ccl19, Cldn11, Thbd and others. 119 are dura-type, 14 leptomeningeal and 18 perivascular.

## Differential expression on fibroblast genes (fibroblast-normalized)
DESeq2 was run on the 151 genes only, so the size factors come from fibroblast genes.

| Contrast | `~sex + group` | `~sex + pct_mito + group` |
|---|---|---|
| 5xFAD NTC vs WT NTC | 57 / 151 | **0** |
| 5xFAD SPP1 vs 5xFAD NTC | 3 | 5 |

**5xFAD vs WT:** restricting to fibroblast genes does not remove the batch-like axis. It acts on every cell type's RNA. The genotype comparison needs batch-balanced samples.

**SPP1 vs NTC:**

| Gene | log2FC | padj (base) | padj (+mito) |
|---|---|---|---|
| Cldn11 (arachnoid barrier) | −0.33 | 0.048 | 0.039 |
| Thbd | −0.36 | 0.048 | 0.039 |
| Matn4 | +0.55 | 0.048 | 0.039 |
| Fbln1 | −0.36 | 0.14 | 0.039 |
| Ifit3b (interferon) | +0.66 | 0.094 | 0.001 |

These genes were also the top hits with the earlier gene set and in the whole-tissue analysis. Cldn11, Fbln1 and Thbd are lower in SPP1 in both sexes.

**Program-level matrix (ECM) change.** Within this gene set, GSEA finds no gene set at FDR < 0.25. A score of the matrix genes (Fbln1, Col3a1, Mmp2, Lum, Col16a1) is lower in SPP1 (effect −0.31, p = 0.097) and correlates with Spp1 within 5xFAD (ρ = 0.57, p = 0.055). Both are suggestive but not significant, and the effect shrinks further after adjusting for osteoblast content (`ecm_score_spp1_robustness.csv`).

## Revision note
The first version of this analysis had two flaws:
1. Two mitochondrial genes were used as fibroblast deconvolution markers. Since mitochondrial reads are the technical axis, this inflated fibroblast estimates in the 5xFAD libraries.
2. Nasal gland genes were misassigned to plasma cells, so 0453 appeared to be 47% plasma cells.

Both are now fixed, and two bone genes (Bglap, Bglap2) that had leaked into the fibroblast set are blacklisted. The fibroblast gene set shrank from 443 to 151 genes. The earlier "ECM/EMT program down with SPP1 (FDR 0.006)" result did not survive, and is withdrawn as a firm finding. The gene-level hits (Cldn11, Thbd, Fbln1, Matn4, Ifit3b) survived.

## Bottom line and next steps
* **Robust:** SPP1 shifts a handful of fibroblast genes, pointing at the arachnoid barrier (Cldn11) and at matrix and endothelial-interface genes (Fbln1, Thbd). These are modest effects at n = 6 per group.
* **Not robust:** a broad fibroblast matrix-program change, and anything about 5xFAD vs WT.
* **Validate with:** a cleanly peeled dura or sorted PDGFRα+ fibroblasts; CLDN11 whole-mount staining; knockdown measured in the targeted cell type.

## Files
| File | Contents |
|---|---|
| `deconvolution_proportions.csv`, `deconvolution_grouped_categories.csv` | Per-sample composition, all classes and grouped |
| `deconvolution_fit_quality.csv`, `deconvolution_markers.csv` | Fit quality per sample and the marker genes used |
| `fibroblast_specific_genes.csv`, `fibroblast_specificity_all_genes.csv.gz` | Fibroblast gene set with signal shares |
| `DE_*`, `DE_summary.csv`, `GSEA_*` | Fibroblast-normalized DESeq2 and GSEA |
| `ecm_score_*.csv` | ECM score, correlation with Spp1, robustness models |
| `F1`–`F7` figures | Reference heatmap, deconvolution (F2, F2b), technical axis, volcanoes, GSEA, heatmaps, ECM score |
