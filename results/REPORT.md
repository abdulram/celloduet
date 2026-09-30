# VT3YXB bulk RNA-seq: 5xFAD vs WT, SPP1-targeted vs non-targeting control

**Design (18 libraries, all sex labels confirmed by Xist / chrY genes):**

| Group | n | Sexes |
|---|---|---|
| WT NTC | 6 | 3 F / 3 M |
| 5xFAD NTC | 6 | 3 F / 3 M |
| 5xFAD SPP1 | 6 | 3 F / 3 M |

There is no WT SPP1 group, so a genotype × treatment interaction can't be estimated.

Code: `analysis/run_analysis.py` (main pipeline) and `analysis/confound_check.py` (technical-confound diagnostics).
Tools: PyDESeq2 (`~ sex + group`, Wald tests, BH FDR) and gseapy preranked GSEA on the Wald statistic, using MSigDB Hallmark v7.0 (human symbols upper-cased) plus custom mouse sets.
19,622 genes pass filtering (≥10 counts in ≥3 samples, small/structural RNAs removed).

---

## Key findings

### 1. The tissue is a brain-border / skull preparation, not brain parenchyma, and its composition varies a lot between samples
The top expressed genes are Ttr, Col1a1/2, Mgp, Ptgds, Dcn and hemoglobins. Brain genes (Snap25, Plp1) are low, while bone and osteoclast genes (Bglap, Ibsp, Ctsk, Acp5, Mmp13) and nasal/olfactory genes (Omp, Cyp2g1, Ugt2a1, Bpifa1) are clearly present.
The profile fits **meninges/dura on skull (calvaria), with variable carry-over of choroid plexus, bone marrow and olfactory mucosa**.

* **PC1 (40% of variance) is mostly olfactory-mucosa contamination** (Cyp2a5, Cyp2g1, Ugt2a1, Scgb1c1; ρ = −0.68 with an olfactory score). 0453 and 613 carry the most, which is why 0453 is the outlier in your correlation heatmap.
* Other large sources of variation: choroid plexus (Ttr ranges from 650 to 30,000 CPM), osteoclast/bone content, blood and skeletal muscle (`figures/04_composition_modules.png`).
* Neither genotype nor treatment drives PC1–PC5 cleanly (`figures/02_pca.png`, `figures/03_pc_covariate_correlation.png`).

### 2. The 5xFAD-vs-WT difference is confounded with a technical library effect. Treat it as not interpretable as-is
The base model gives **6,811 DE genes (padj < 0.05)** for 5xFAD NTC vs WT NTC, but the signature looks technical:

| Median per group | WT NTC | 5xFAD NTC | 5xFAD SPP1 |
|---|---|---|---|
| % uniquely mapped | **87.4** | 80.5 | 81.3 |
| % multi-mapped | **10.3** | 16.8 | 16.5 |
| % mitochondrial | **5.2** | 10.0 | 9.9 |

* Five of six WT libraries (611, 624, 626, 0466 and, partly, 623) share a distinct mapping profile. **WT 616 looks like the 5xFAD libraries and scores like 5xFAD on the "genotype" signature** (`figures/10_genotype_signature_vs_technical.png`; signature vs % mito r = 0.94, vs % multi-mapped r = 0.90).
* The genes "up in 5xFAD" include 253 pseudogenes, 167 ribosomal-protein genes and 12 mitochondrial genes. None of these classes appear among the down genes. The "down in 5xFAD" genes are enriched for long and nuclear-retained transcripts (Neat1, Prrc2a, Hspg2, Srrm1, Kdm6b). This is a classic RNA-quality / extraction-batch pattern.
* Adding % mito as a single covariate cuts the count from 6,811 to **1**; % multi-mapped cuts it to **379**. Composition covariates (osteoclast, choroid plexus, olfactory) change almost nothing (`genotype_DE_sensitivity.csv`).
* The GSEA "5xFAD" hits (up: OXPHOS, MYC targets, mTORC1; down: apical junction, WNT, TGF-β) match that technical axis.

**Caveat on the other side:** % mito is itself 80% correlated with genotype, so adjusting for it removes any real biology that happens to align with it. The data can't separate genotype from the technical factor. **Check your lab records** for whether the WT samples were dissected, extracted or library-prepped separately, and check the RIN/DV200 values.

### 3. SPP1 vs NTC (within 5xFAD): little gene-level signal and no significant Spp1 reduction in this tissue
* **Spp1 is not significantly reduced:** base LFC −1.2, padj 0.67 (Mann-Whitney p = 0.24). The median is ~33% lower, but the three SPP1 **males** (0453, 612, 613: median 1,476 CPM) drive that. The SPP1 females sit at NTC levels (3,891 CPM).
* **Spp1 in this tissue tracks osteoclast/bone content** (ρ = 0.60 with an osteoclast module). Its top genome-wide correlates are Atp6v0d2, Ctsk, Acp5, Mmp13, Ibsp and Snx10. NTC 620's extreme Spp1 (23,000 CPM) comes with the highest osteoclast score. **Bulk Spp1 here mostly reports bone and osteoclasts, not microglia or macrophages,** so a knockdown targeted at CNS myeloid cells would be hard to see in this preparation (`figures/05_spp1_vs_osteoclast.png`).
* Only **4 genes** reach padj < 0.05, and they are modest: **Rgs4** (LFC −0.91), **Cldn11** (−0.44; arachnoid barrier cells/oligodendrocytes), **Prg4** (−0.65) and **Bpifa1** (+5.3). Bpifa1 is a nasal gland gene, so that hit is olfactory contamination. After adjusting for composition and technical covariates, Rgs4, Prg4, Cldn11 and Man1a remain just above the cutoff (padj ≈ 0.056–0.058).
* **GSEA is the most interesting result, and the least certain.** It is consistent between the base and adjusted models. In SPP1 vs NTC:
  * **Up:** interferon/ISG (Ifit1, Ifit3, Isg15), IFN-α/γ response and E2F/G2M.
  * **Down:** EMT/ECM, TGF-β, angiogenesis, coagulation, complement, TNF-α/NF-κB and the DAM microglia set.

  That fits SPP1 loss leading to less fibrotic/ECM and DAM-like activation and more type-I IFN tone. However, preranked GSEA permutes genes, not samples, so with n = 6 vs 6 its q-values are optimistic. The DAM set also contains Spp1 itself. **Treat this as a hypothesis to test**, not a result.

### 4. QC is otherwise fine
Depth is 12–21 M reads per library with 78–88% uniquely mapped and about 18–23k genes detected. Sex genes agree with all labels. The ERCC spike-ins have zero counts, so they aren't usable for normalization. The 5xFAD transgenes (human APP/PSEN1) aren't in the mouse reference, so genotype can't be confirmed from these data.

---

## Recommendations
1. **Resolve the WT batch question first.** If the WT tissue was processed separately, the genotype comparison needs new, batch-balanced samples. Covariate adjustment can't rescue a near-perfect confound.
2. **Standardize the dissection**, or remove contaminating tissue computationally. Olfactory mucosa, choroid plexus and bone marrow content vary 10–100× between animals. If the target is dura/meninges, strip it off the skull. If the target is myeloid cells, sort them (e.g. CD11b+).
3. **Measure knockdown in the target cell type** (qPCR or sorted cells) instead of relying on bulk Spp1, which here is dominated by bone.
4. **Validate the SPP1 → IFN-up / ECM- and DAM-down pattern** with sample-permutation GSEA or a larger n, and in sorted cells.
5. Consider analyzing each sex separately for the SPP1 effect, since the Spp1 drop is only seen in males (n = 3 per cell, so this is exploratory).

## Files
| File | Contents |
|---|---|
| `qc_and_composition.csv` | Per-sample QC, sex check and composition module scores |
| `DE_<contrast>_<base\|adjusted>.csv` | Full DESeq2 results. The adjusted model is `~ sex + osteoclast + choroid plexus + olfactory + %mito + group` |
| `DE_summary.csv` | DE gene counts per contrast and model |
| `GSEA_<contrast>_<model>.csv` | Preranked GSEA results |
| `genotype_DE_sensitivity.csv`, `genotype_signature_gene_classes.csv` | Confound diagnostics |
| `pca_scores.csv`, `pca_loadings_top2000.csv`, `pc_covariate_spearman.csv` | PCA |
| `figures/01`–`10` | QC, PCA, composition, Spp1, volcanoes, heatmap, GSEA, marker panel, confound plot |

---

## Follow-up: fibroblast-specific re-analysis
See **`results/fibroblast/FIBROBLAST_REPORT.md`**. In short:
- Using a CELLxGENE single-cell reference plus deconvolution, 443 genes were selected whose bulk signal comes mostly from fibroblasts, and DE was re-run on them with fibroblast-based normalization.
- The genotype contrast is still driven by the technical axis.
- SPP1 vs NTC shows a modest, sex-consistent reduction in fibroblast matrix/barrier genes (Cldn11, Fbln1, Loxl3, collagens; EMT/ECM GSEA FDR 0.006). Part of this effect overlaps with osteoblast content.
