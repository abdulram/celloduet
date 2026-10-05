# 5xFAD dura bulk RNA-seq: dissection / cell-composition confounding analysis

**Purpose.** Determine which 5xFAD-associated expression changes in this dura dataset are robust to differences in dissected material and cell composition, and which are explained by them. This analysis is deliberately conservative. A negative result is reported as a result.

**One-line answer.** We cannot reliably estimate any 5xFAD effect in this dataset.
- Genotype is associated with the major composition axes, with poor overlap between genotypes.
- Genotype is even more strongly confounded with a library-processing (technical) axis.
- The single WT animal that shares the 5xFAD libraries' technical profile (616) looks "5xFAD" on 84% of the baseline genotype genes.
- No pathway reaches Tier 1 or Tier 2.

---

## Dataset
| | |
|---|---|
| Samples | 18 libraries (facility run VT3YXB), one per animal |
| Groups | WT_NTC n = 6 (3F/3M); 5xFAD_NTC n = 6 (3F/3M); 5xFAD_SPP1 n = 6 (3F/3M) |
| Primary genotype contrast | **5xFAD_NTC vs WT_NTC**. SPP1-treated 5xFAD animals are kept in models for dispersion estimation and shown for context. There is no WT_SPP1 group. |
| Age | **Not available** |
| Batch / dissection date / dissector | **Not available** |
| Pathology (plaques, Aβ, gliosis, behavior, vascular/lymphatic measures) | **Not available** (searched every data file; `09_pathology_associations/pathology_metadata_search.tsv`) |
| Sex | Available; all labels confirmed by Xist / chrY genes |
| Animal IDs | Available. IDs fall in two numbering series, "04xx" and "6xx". This is recorded as an *unverified* possible cohort: WT 1/6 in 04xx vs 5xFAD 6/12. |
| Transgene verification | Not possible: the count matrix is mouse-only, plus ERCC spike-ins, which have 0 reads. |

## Analysis environment
- **R** 4.3.3, with DESeq2 1.42.0, limma 3.58.1, MatchIt 4.5.5, car 3.1.2, ggplot2 3.4.4 and pheatmap 1.0.12.
- **Python** 3.11.15, with numpy 2.4.6, pandas 2.3.3, scipy 1.17.1, statsmodels 0.15.0, matplotlib 3.11.2 and tiledbsoma 2.3.0 (CELLxGENE Census access).
- Full versions: `logs/package_versions.txt` and `logs/sessionInfo_*.txt`. The random seed is 7 everywhere.

## Input files (read-only)
| Path | Content |
|---|---|
| `data/VT3YXB-expression-matrix.tsv` | Gene × sample raw counts and CPM (Ensembl mouse IDs, gene names, biotypes) |
| `data/VT3YXB-mapping-stats.csv` | Uniquely mapped, multi-mapped and unmapped reads per sample |
| `data/VT3YXB-gene-biotype-summary.csv` | Facility biotype QC |
| `data/genesets/h.all.v7.0.symbols.gmt` | MSigDB Hallmark v7.0 (human symbols, mapped to mouse by symbol) |
| CELLxGENE Census 2025-01-30 (public S3) | Single-cell reference: adult mouse dura dataset `58b01044-…` plus other adult/embryonic data for classes the dura dataset lacks (`scripts/00_build_dura_reference.py`) |

No original files were modified. There were no previous DESeq2 objects. Earlier exploratory analyses in this repository (`results/`, `analysis/`) were not reused as inputs; their conclusions are consistent with this analysis.

## QC findings (`01_QC/`)
- **Depth and gene detection are fine:** 12–21 M reads per library, 18–23k genes with ≥5 reads.
- **Two library profiles exist, and they split almost perfectly by genotype:**

  | Median | WT_NTC | 5xFAD_NTC | 5xFAD_SPP1 |
  |---|---|---|---|
  | % uniquely mapped | 87.4 | 80.5 | 81.3 |
  | % multi-mapped | 10.3 | 16.8 | 16.5 |
  | % mitochondrial | 5.2 | 10.0 | 10.0 |

  All 12 5xFAD libraries have the high-multimapping / high-mito profile, as do WT **616** (fully) and **623** (intermediate: multi-mapping 14.3%, mito 6.0%).
- **No samples were removed.** 0453 (heavy olfactory/nasal carry-over) and 620 (heavy osteoclast/bone content) are extreme on composition, but there was no objective QC failure.
- **PCA by age or batch is impossible** (no metadata). It was replaced by PCA colored by ID series and by library profile.

## Major sources of variation (`01_QC/transcriptome_PC_top_genes.tsv`, `tables/txPC_variance_explained_by_composition.tsv`)
| Tx PC | % variance | What it is | R² by Composition_PC1–3 | R² by technical axis | R² by group |
|---|---|---|---|---|---|
| PC1 | 40.2 | Dura purity: olfactory/nasal (Cyp2a5, Cyp2g1, Scgb1c1) vs fibroblast/endothelial | **0.90** | 0.60 | 0.45 |
| PC2 | 22.9 | Bone (Composition_PC3, ρ = −0.89) + multi-mapping (ρ = −0.81) + muscle | **0.83** | 0.50 | 0.34 |
| PC3 | 8.5 | Skeletal-muscle carry-over (Myh4, Acta1, Ckm) | 0.02 | 0.31 | 0.03 |
| PC4 | 5.5 | Sex (Xist, Ddx3y) | 0.46 | 0.01 | 0.01 |
| PC5 | 3.9 | Blood / marrow neutrophils (S100a8/9, Ngp, Camp) | 0.25 | 0.01 | 0.20 |

**About 55% of all transcriptome variance (PC1 + PC2) is explained by measured composition.** The next components are muscle carry-over, sex and blood. Genotype is not a primary axis.

## Cell composition (`02_cell_composition/`, `tables/signature_validation.tsv`)
- **Primary method: marker-signature scoring.** We defined 20 broad signatures: 13 native dura lineages plus 7 adjacent/contaminating tissues. Each has 5–8 canonical markers and is scored as the mean z of the blind VST.
  - **Validation rule:** at least 3 markers detected and median marker-vs-rest r ≥ 0.4.
  - **16 validated:** fibroblast/stromal, blood endothelial, macrophage/BAM, monocyte, dendritic, B, NK, mast, neutrophil, CNS, bone, osteoclast, olfactory/nasal, choroid plexus, erythroid, skeletal muscle.
  - **Not quantifiable** (markers do not co-vary): lymphatic endothelial (r = 0.33), mural/pericyte/SMC (0.09), T cell (−0.03), plasma cell (0.22).
- **Secondary method: reference deconvolution.** We used a 20-class dura-centric CELLxGENE pseudobulk reference with a trimmed NNLS ensemble.
  - **Unstable for most classes.** It assigns a median 44% to choroid plexus (driven by very high Ttr; implausible), puts immune classes near 0%, and agrees poorly with the marker scores for immune cells (ρ ≤ 0.45).
  - Only fibroblast, blood endothelial, osteoblast, osteoclast, erythroid, choroid plexus and muscle are both stable and concordant (`tables/deconvolution_validation.tsv`).
  - **The reference does not adequately represent dura.** Lymphatic and mural cells come from other tissues, and the osteoclast, choroid plexus and olfactory profiles are embryonic. Deconvolution percentages are not used for inference.
- **Signatures show marker-RNA abundance, not cell counts.** `02_cell_composition/lineage_marker_heatmap.png` shows every marker per animal, ordered by the composition axis.

## Genotype–composition confounding (`03_composition_PCA/`, `tables/genotype_vs_composition.tsv`)
**Composition PCA** (16 validated signatures, standardized; explained variance 26.9 / 22.7 / 16.8%):

| Axis | What it represents | 5xFAD_NTC vs WT_NTC | AUC | Common support |
|---|---|---|---|---|
| Composition_PC1 | Dura purity: fibroblast + blood endothelial + DC/B vs olfactory/CNS/muscle | Hedges g −1.62, p = 0.015 | 0.08 | **Poor** (30% range overlap) |
| Composition_PC2 | Myeloid/NK/mast | g 0.23, p = 0.59 | 0.61 | Adequate |
| Composition_PC3 | Bone (osteoblast + osteoclast) | g +1.63, p = 0.009 | 0.94 | **Poor** (12% overlap) |
| % mito (technical) | Library / RNA profile | g +2.26, p = 0.041 | 0.86 | Only WT 616 overlaps |

Composition and the technical axis are correlated across all animals (ρ = −0.55) but **not within 5xFAD** (ρ = 0.10). They are distinct axes, each confounded with genotype.

**Interpretation.** In 5xFAD animals, the dissected material contained relatively less dura stroma and blood endothelium and more bone and adjacent tissue. This could reflect (A) disease remodeling, (B) systematic differences in dissection, or (C) both. **Bulk data cannot distinguish these.** WT 616, the one WT processed like the 5xFAD libraries, also sits with the 5xFAD animals on Composition_PC1 (−0.63), which favors (B).

> **FLAG — common support:** WT and 5xFAD overlap poorly on Composition_PC1 and Composition_PC3, and essentially not at all on the technical axis. Statistical adjustment cannot be assumed to remove these confounds.

## Baseline genotype analysis (`04_baseline_DE/`)
- **Model:** DESeq2 `~ sex + group`, 18 animals, contrast 5xFAD_NTC vs WT_NTC. Age and batch could not be modeled.
- **Genes:** 6,787 at padj < 0.05 (3,401 up, 3,386 down).
- **Pathways** (limma camera on VST, accounting for correlation between genes): no gene set has FDR < 0.9. The top-ranked sets are the technical sentinels (cytosolic ribosome p = 0.08, mitochondrial-encoded p = 0.09, OXPHOS p = 0.11).
- This is the **unadjusted association**, not a disease effect.

## Composition-adjusted analysis (`05_composition_adjusted_DE/`, `tables/adjusted_models_summary.tsv`)
All models are full rank. The variance inflation for genotype is 1.5–2.0, so the models are estimable but the covariates are correlated with genotype (r = −0.69 to 0.80).

| Model | Genes padj < 0.05 | Median effect retained | Baseline genes keeping ≥ 50% of effect |
|---|---|---|---|
| `~ sex + Composition_PC1 + group` | 1,678 | 0.72 | 5,692 / 6,787 |
| `~ sex + Composition_PC1 + Composition_PC2 + group` | 426 | 0.63 | 5,026 |
| `~ sex + Composition_PC1 + Composition_PC3 + group` | **0** | 0.53 | 3,693 |
| `~ sex + pct_mito + group` (technical) | **1** | **0.32** | 2,025 |
| `~ sex + Composition_PC1 + pct_mito + group` | 2 | 0.30 | 1,891 |

- **Dura purity alone explains a minority of the baseline signature.** Adding the bone axis removes all significance.
- **The technical axis explains about two-thirds of the effect size.**
- 127 genes are "potentially masked" under Composition_PC1 only (`tables/genotype_effects_baseline_vs_adjusted.tsv.gz`). None persist in the other models.

## Common-support / matched analysis (`06_common_support_matching/`)
All algorithms were pre-specified (`scripts/06_common_support_analysis.R`).

| Cohort | Rule | Animals | Result |
|---|---|---|---|
| Restricted | Composition_PC1 and PC3 both inside the shared range | WT 616 vs 5xFAD 451 | Too few for DE |
| Matched | MatchIt 1:1 nearest-neighbour, Mahalanobis on Composition_PC1 + PC3, caliper 1 SD | WT 616, 626 vs 5xFAD 617, 451 | 4 genes significant. 95% direction concordance with baseline but 0 baseline genes significant (n = 2 vs 2). |
| Technical profile B | High-multimapping libraries | WT 616, 623 vs all 6 5xFAD_NTC | 99% direction concordance, 9 genes significant. 623 is only intermediate on the technical axis (mito 6.0%). |
| **Animal-position test (06b)** | Each WT placed between other WT (0) and 5xFAD (1) on all 6,787 baseline genes | — | **WT 616 median 0.96, 5xFAD-like on 84% of genes. Other WT −0.36 to −0.16 (2–5%).** |

**Interpretation.** The baseline 5xFAD signature follows the processing/technical axis, not the genotype label. There are two possibilities, and both make genotype non-estimable here:
1. 616 is truly WT, so the "genotype" signature is technical or batch.
2. 616 is a mis-genotyped 5xFAD animal, so genotype and library profile are perfectly confounded.

**Genotyping 616 is the single most informative check available.**

## Cell-state analysis (`07_cell_state_analysis/`, `tables/cell_state_vs_abundance.tsv`)
- **Model:** `state ~ lineage_abundance + sex + genotype` (n = 12). Any gene that is an abundance marker is removed from the state module. p-values come from OLS and from an exact within-sex permutation (400 relabelings).
- **Nominal hits:**
  - Macrophage MHC-II lower in 5xFAD at equal macrophage abundance (β −0.83 SD; permutation p = 0.02). Lost with the % mito covariate (p = 0.66).
  - Meningeal-trafficking module given endothelial abundance: permutation p = 0.04, OLS p = 0.19.
  - DAM and fibroblast-activation: permutation p ≈ 0.08–0.10.
- **None survive** correction across 15 tests or the technical adjustment.
- Endothelial angiogenesis differs before adjustment (p = 0.008) but is almost entirely abundance (r = 0.95 with endothelial abundance).
- T-cell, lymphatic, plasma and mural *state* cannot be analyzed, because their abundance cannot be measured.
- **There is no evidence for altered cell state beyond abundance.** This is an inferred analysis, not single-cell evidence.

## 5xFAD-relevant findings (`08_pathways/`, `tables/pathway_robustness_matrix.tsv`)
Module scores were tested under baseline, +CompPC1, +CompPC1+PC3, +% mito and matched analyses.

- **Inflammation, interferon, complement, antigen presentation, chemokines, monocytes, BAM programs:** no baseline genotype association (all p > 0.1).
- **Myeloid / DAM:** DAM module +1.09 SD, p = 0.068; the trend persists across adjustments but is not significant. Its genes are incoherent: Itgax, Axl and Apoe go down while Spp1, Fabp5 and Tyrobp go up. Spp1 here is mostly bone/osteoclast-derived.
- **Lipid handling:** +0.74 SD at baseline (p = 0.15), +1.6 SD after composition adjustment (p = 0.015–0.04), +1.5 SD with % mito (p = 0.08). Driven mostly by Fabp5 (and Abca1), with Apoe and Lipa opposite. Camera FDR 0.99. **Potentially masked, suggestive at most.**
- **Lymphocytes:** T-cell and plasma signatures cannot be measured. B-cell abundance is nominally lower (AUC 0.19, p = 0.09).
- **Vascular:** blood-endothelial abundance is lower in 5xFAD (p = 0.015; composition, Tier 3). The angiogenesis module is Tier 4 (616 position 1.13).
- **Fibroblast / ECM:** fibroblast activation +1.15 SD (p = 0.049) is explained by the bone axis (5% of the effect retained), so Tier 3. ECM remodeling shows no difference.
- **Lymphatics:** not measurable (incoherent markers).
- **Technical sentinels** (ribosome, mitochondrial-encoded, OXPHOS) are the strongest "genotype" modules (+1.4 to 1.5 SD, p ≤ 0.01), with 616 positioned at about 1.1. **This is the clearest evidence that the baseline genotype contrast carries technical signal.**

## Pathology relationships
None were possible: no pathology, severity or age data exist (`09_pathology_associations/`).

## Robust findings (Tier 1 and Tier 2)
**None.**

## Composition-sensitive findings (Tier 3)
- Genotype is associated with lower Composition_PC1 (dura purity) and higher Composition_PC3 (bone).
- Blood-endothelial, bone, dendritic-cell and B-cell signatures differ by genotype.
- Fibroblast activation module.
- These may be real biology, but this experiment cannot separate them from dissection.

## Uninterpretable findings (Tier 4)
- The genome-wide baseline 5xFAD signature (6,787 genes).
- Angiogenesis-endothelial and ROS modules.
- The nominal cell-state hits.
- All statements about lymphatic, mural, T-cell or plasma-cell biology.

## Limitations
- **Genotype is confounded with library profile** (all 5xFAD vs 1 WT) and poorly overlaps on two composition axes.
- **There is no batch, age, dissection or pathology metadata.** Genotype cannot be verified from the count matrix.
- **n = 6 per primary group.** Matched and restricted cohorts have n = 2 vs 2 or smaller.
- **Deconvolution is unstable** and the reference does not fully represent dura. Four lineages cannot be measured.
- **Gene sets are mostly small, hand-curated modules.** Hallmark symbols were mapped from human by case.
- **Module-score p-values are not corrected across the 38 modules.** Camera, which is corrected and correlation-aware, finds nothing.

## Conclusions
Once plausible differences in dissection, composition and library processing are accounted for, **this dataset supports no conclusion about 5xFAD-associated dural biology.**
- **What the baseline genotype contrast reflects:** mainly a technical/processing axis and a dura-purity/bone composition axis, both strongly confounded with genotype. The decisive observation is that the only WT animal sharing the 5xFAD library profile (616) looks "5xFAD" genome-wide.
- **Suggestive at most:** lipid-handling (Fabp5-driven) and DAM-like trends.

## Recommended follow-up experiments
1. **Genotype WT 616** (and confirm all animals), and retrieve processing records (RNA extraction and library-prep dates, RIN/DV200).
2. **Balanced re-processing.** Process WT and 5xFAD in the same extraction and library batches, interleaved, and record dissector and date.
3. **Standardized region-specific dura sampling**, for example peeling the dura from the calvaria and excluding the olfactory/cribriform region and bone. Spike a marker check into QC.
4. **Flow cytometry** of the dura (CD45, CD206/LYVE1, CD3, CD19, CD31, PDGFRα, PROX1) to measure composition directly.
5. **Whole-mount immunofluorescence** (LYVE1/PROX1 lymphatics, CD31, CD206, CD3, Aβ) to test remodeling in situ.
6. **scRNA-seq/snRNA-seq** of dura (≥3 per genotype) to separate abundance from state, especially for lymphatic and T cells.
7. **Spatial transcriptomics** around sinuses and lymphatics.
8. **Collect pathology per animal** (plaque load, Aβ40/42, Iba1/GFAP) and age, so dural features can be related to severity within 5xFAD.
9. **An independent replication cohort.**

## File index
| What | Path |
|---|---|
| Scripts (run in order) | `scripts/00_build_dura_reference.py`, `01_load_and_qc.R`, `02_cell_signatures.py`, `03_composition_analysis.py`, `03b_define_gene_sets.py`, `04_baseline_deseq2.R`, `05_composition_adjusted_deseq2.R`, `06_common_support_analysis.R`, `06b_animal_position_test.py`, `07_cell_state_analysis.py`, `08_pathway_analysis.py`, `09_pathology_analysis.py`, `10_generate_summary.py` (helpers: `common.R`, `common.py`, `de_helpers.R`) |
| Sample table (metadata, QC, Tx PCs, composition PCs) | `tables/sample_metadata_with_composition.tsv` |
| Composition | `tables/composition_signature_scores*.tsv`, `signature_validation.tsv`, `deconvolution_*.tsv`, `genotype_vs_composition.tsv` |
| DE results | `tables/DE_baseline_*.tsv`, `DE_comp1_*`, `DE_comp12_*`, `DE_comp13_*`, `DE_tech_*`, `DE_comp1tech_*`, `DE_matched_*`, `DE_tech_profile_B_*` |
| Baseline vs adjusted per gene | `tables/genotype_effects_baseline_vs_adjusted.tsv.gz` |
| Pathways | `tables/pathways_camera_*.tsv`, **`tables/pathway_robustness_matrix.tsv`** |
| Cell state | `tables/cell_state_vs_abundance.tsv` |
| Technical-axis test | `tables/technical_axis_animal_position_test.tsv` |
| Final figures | `10_final_figures/Fig01`–`Fig12` |
| Summaries | `analysis_summary.tsv`, `EXECUTIVE_SUMMARY.md` |
| Logs and versions | `logs/` |
