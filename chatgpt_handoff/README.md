# Handoff: bulk RNA-seq of mouse dura, 5xFAD vs WT and SPP1-targeted vs non-targeting control (run VT3YXB)

This package has everything needed to understand, re-plot and critique the analyses done on this dataset:
- **Study summary and design** (Sections 1–2).
- **Data processing** (Section 3).
- **Each analysis**, with methods, results and caveats (Sections 4–11).
- **Analysis history**, including mistakes that were found and corrected (Section 12).
- **Overall conclusions, open questions and suggested experiments** (Sections 13–14).
- **Plotting guide**: which data file reproduces each figure (Section 15).
- **Data dictionary and conventions** (Sections 16–17).

`FILE_MANIFEST.csv` lists every data file with its row count and column names. `figures/` contains the original figures as visual references.

> **Note for the reader (human or AI):** Section 13, the overall conclusions, is the authoritative interpretation. Several early results were later revised or withdrawn (Section 12). Do not report a withdrawn result as a finding. The most important caveat: **genotype (5xFAD vs WT) is confounded with library processing and dissection composition**, so 5xFAD-vs-WT differences generally cannot be interpreted as disease biology.

---

## 0. PRIMARY QUESTION: 5xFAD_SPP1 vs 5xFAD_NTC (`data/10_SPP1_vs_NTC_PRIMARY/`, figures F_*)
The main interest of this study is the **SPP1-targeted vs non-targeting-control comparison within 5xFAD mice** (n = 6 vs 6, 3F/3M each).

**This is the cleanest comparison in the dataset:**
- SPP1 and NTC libraries share the same technical profile: % mito 9.9 vs 10.0, % multi-mapped 16.5 vs 16.8.
- So the processing confound that ruins 5xFAD vs WT does not apply here.

### 0.1 Methods specific to this contrast
- **DESeq2 genome-wide** (19,622 genes) under 5 models:

  | Model | Formula | Animals |
  |---|---|---|
  | baseline | `~ sex + group` | all 18 (contrast SPP1 vs NTC) |
  | comp1 | `~ sex + Composition_PC1 + group` | all 18 |
  | comp13 | `~ sex + Composition_PC1 + Composition_PC3 + group` | all 18 |
  | tech | `~ sex + pct_mito + group` | all 18 |
  | 5xFADonly | `~ sex + treatment` | 12 5xFAD |

  Output: `SPP1_vs_NTC_all_genes_all_models.csv.gz` has LFC, p and padj for every model, plus `n_models_padj05` and `sign_consistent`.
- **Composition balance:** does SPP1 vs NTC differ in dissected tissue? Each composition signature, composition PC and technical metric is tested with an exact within-sex permutation (all 20 × 20 = 400 relabelings of SPP1/NTC that keep 3 F and 3 M per group) (`SPP1_vs_NTC_composition_balance.csv`).
- **Module scores** (38 gene sets, mean z of VST): SPP1 effect in SD units, adjusted for (i) nothing, (ii) Composition PC1 + PC3, (iii) % mito, (iv) fibroblast abundance (`signature__fibroblast_stromal`). Plus the exact permutation p and limma camera (`SPP1_vs_NTC_module_effects.csv`, `SPP1_vs_NTC_camera_pathways.csv`).
- **Key genes:** the same set of adjustments, plus per-sex median differences and the permutation p (`SPP1_vs_NTC_key_genes.csv`).

### 0.2 Results
**Composition balance (figure F_SPP1_A).**
- SPP1 and NTC are balanced on Composition PC1/2/3, % mito, % multi-mapped, bone, osteoclast, olfactory, blood and all immune signatures (permutation p ≥ 0.15).
- **Exception: the fibroblast/stromal signature is lower in SPP1** (−1.11 SD, p = 0.025; AUC 0.14).
- So SPP1 samples contain relatively less fibroblast RNA. This could be a treatment effect (fewer or less active fibroblasts) or dissection; bulk data cannot distinguish them.

**Genes robust in every model** (padj < 0.05 in ≥ 4 of 5 DESeq2 models; also robust to fibroblast-abundance adjustment; permutation p = 0.005; down in both sexes; figure F_SPP1_B):

| Gene | log2FC | What it is |
|---|---|---|
| **Cldn11** | −0.44 | Claudin-11, tight junction of arachnoid barrier cells |
| **Prg4** | −0.65 | Lubricin / proteoglycan 4 |
| **Rgs4** | −0.91 | Regulator of G-protein signaling 4 |

**Robust to composition and fibroblast adjustment, but significant in fewer DESeq2 models:**
- Wfdc17 down (−0.77; padj < 0.05 in tech and 5xFAD-only models).
- **Ifit3b up** (+0.47) and **Matn4 up** (+0.45). Both have permutation p ≤ 0.05 and stay significant after fibroblast adjustment (p = 0.017 and 0.006).

**Matrix / mechanotransduction genes: down, but about half the effect is explained by lower fibroblast signal.**
- Fn1, Fbln1, Thbd, Loxl3, Col12a1 and Wwtr1 (TAZ): log2FC −0.39 to −0.58, permutation p = 0.01–0.03, mostly stronger in males.
- They reach padj < 0.05 only in the % mito model.
- After fibroblast-abundance adjustment the effects shrink by about 50% (p = 0.05–0.15).
- Fbln1 is itself a marker in the fibroblast signature, so that adjustment is partly circular for Fbln1.

**Pathways (figure F_SPP1_C).**
- Unadjusted module scores: SPP1 lowers Hallmark TGF-β (−1.10 SD, permutation p = 0.035), fibroblast activation (−1.14, p = 0.045), DAM (−1.13, p = 0.06), EMT, angiogenesis, hypoxia and lipid handling (p = 0.07–0.09). It **raises type-I interferon/ISG** (+1.05, p = 0.055).
- Camera FDR is 0.51 for all of these.
- **After fibroblast-abundance adjustment,** the TGF-β/EMT/ECM/DAM decreases mostly disappear (p 0.3–0.8). The **interferon increase strengthens** (+1.63 SD, p = 0.019).
- Composition PC1 + PC3 or % mito adjustment leaves most effects in place.

**Spp1 itself:** not significantly lower (log2FC −1.2, permutation p = 0.115). The drop is only in males (VST −1.40 in M vs +0.01 in F). Bulk Spp1 mostly reflects osteoclasts (Section 4), so knockdown in the targeted cells cannot be assessed.

**SPP1 and 5xFAD disease state:** SPP1 does **not** move 5xFAD animals toward or away from WT on any 5xFAD signature (Section 7). The TGF-β ligand Tgfb2, the one robust 5xFAD-vs-WT gene, shifts only slightly and non-significantly toward WT with SPP1 (Section 8).

### 0.3 Interpretation (conservative)
- **Robust SPP1-associated changes** in 5xFAD dura:
  - lower **Cldn11** (arachnoid barrier), **Prg4** and **Rgs4**;
  - higher **Ifit3b** and **Matn4**;
  - a likely increase in **type-I interferon** response.
- **Matrix / TGF-β / fibroblast-activation decreases** (Fn1, Fbln1, collagens, Loxl3, Wwtr1) go together with lower overall fibroblast signal. Either SPP1 reduces fibroblast abundance or activity (a real effect), or SPP1 dissections contained less fibroblast tissue. Immunostaining is needed to decide.
- Effects are about 20–45% and often stronger in males; n = 3 per sex per group.

---

## 1. Study in one paragraph
Dura (dural meninges, dissected from the skull) was collected from 18 mice:

| Group | n | Sexes |
|---|---|---|
| WT treated with a non-targeting control (WT_NTC) | 6 | 3 F / 3 M |
| 5xFAD (Alzheimer's amyloid model) with non-targeting control (5xFAD_NTC) | 6 | 3 F / 3 M |
| 5xFAD with an SPP1 (osteopontin)-targeting treatment (5xFAD_SPP1) | 6 | 3 F / 3 M |

There is no WT_SPP1 group. Bulk RNA-seq was done by a sequencing facility (run ID VT3YXB). The questions were:
1. Are there 5xFAD-associated changes in the dura, especially in dural fibroblasts?
2. Does SPP1 treatment change them?

**Short answer** (the SPP1 vs NTC comparison is the primary interest; see Section 0):
- The 5xFAD-vs-WT comparison is dominated by a **technical/library-processing difference** and by **differences in dissected tissue composition**. Almost no genotype effect survives scrutiny. The exception is lower **Tgfb2** ligand expression, which is a candidate (Tier 2).
- SPP1-vs-NTC is a clean comparison (same processing). It shows a **small, sex-consistent reduction of a few fibroblast/barrier/matrix genes**: Cldn11, Fbln1, Thbd and Fn1 down; Matn4 and Ifit3b up.
- SPP1 does **not** move 5xFAD animals toward or away from WT on any 5xFAD signature.

## 2. Samples and metadata
| library_id | sample (animal_genotype_sex_treatment) | group |
|---|---|---|
| VT3YXB_1–3 | 0453_5xFAD_M_SPP1, 612_5xFAD_M_SPP1, 613_5xFAD_M_SPP1 | 5xFAD_SPP1 (males) |
| VT3YXB_4–6 | 621_5xFAD_F_SPP1, 0463_5xFAD_F_SPP1, 0455_5xFAD_F_SPP1 | 5xFAD_SPP1 (females) |
| VT3YXB_7–9 | 616_WT_M_NTC, 611_WT_M_NTC, 626_WT_M_NTC | WT_NTC (males) |
| VT3YXB_10–12 | 624_WT_F_NTC, 623_WT_F_NTC, 0466_WT_F_NTC | WT_NTC (females) |
| VT3YXB_13–15 | 617_5xFAD_M_NTC, 0452_5xFAD_M_NTC, 0451_5xFAD_M_NTC | 5xFAD_NTC (males) |
| VT3YXB_16–18 | 629_5xFAD_F_NTC, 0469_5xFAD_F_NTC, 620_5xFAD_F_NTC | 5xFAD_NTC (females) |

- **Sex labels:** all confirmed from Xist vs chrY genes (Ddx3y, Uty, Kdm5d, Eif2s3y).
- **Not available:** age, processing batch, dissection date/person, and any pathology (plaques, Aβ, gliosis, behavior).
- **Animal-ID series:** IDs fall in two numbering series, "04xx" (leading zero) and "6xx". This was recorded as a *possible* cohort; it is unverified. WT has 1/6 in 04xx; 5xFAD has 6/12.
- **Animal 616 (WT) is important.** It is the only WT library with the 5xFAD-like technical profile (Section 6).
- In some tables the `animal` column lost the leading zero (453 = 0453). Use `sample` or `library_id` as the key.

## 3. Data and processing
### 3.1 Raw data (from the facility)
- `data/00_core/counts_raw_all_genes.csv.gz`: gene-level read counts for 78,317 Ensembl mouse genes × 18 samples, with gene name and biotype. Alignment and quantification were done by the facility (STAR-style uniquely mapped reads; exact pipeline not provided).
- `data/00_core/cpm_all_genes.csv.gz`: facility counts-per-million.
- `data/01_qc_pca/mapping_stats_reads.csv`: uniquely mapped, multi-mapped and unmapped reads per sample.
- The 92 ERCC spike-ins had 0 reads in every sample and were dropped.
- The reference is mouse only, so **human APP/PSEN1 transgene reads cannot be counted**. Genotype cannot be verified from this data.

### 3.2 Processing used in the analyses
- **Gene filter for modeling:** ≥ 10 counts in ≥ 3 samples, excluding rRNA, Mt_rRNA, Mt_tRNA, misc_RNA, snRNA, snoRNA, scaRNA, miRNA and ribozyme. 19,622 genes remain.
- **Normalization:** DESeq2 median-of-ratios size factors.
- **Transformation:** DESeq2 variance-stabilizing transform, `vst(blind=TRUE)` → `data/00_core/vst_blind_filtered_genes.csv.gz`. Used for PCA, heatmaps, signature and module scores.
- **Transcriptome PCA:** top 2,000 most variable VST genes, centered, not scaled.
- **Differential expression:**
  - DESeq2 1.42 in R 4.3.3 for the composition-confounding analysis (Sections 6–11).
  - PyDESeq2 0.5.4 for the initial and fibroblast analyses (Sections 4, 5, 7).
  - Wald tests, Benjamini-Hochberg FDR.
  - **log2FC sign: first-named group minus second.** For example, "5xFAD_NTC_vs_WT_NTC" > 0 means higher in 5xFAD.
- **Gene sets:**
  - MSigDB Hallmark v7.0 (human symbols mapped to mouse by case-insensitive symbol match).
  - Curated mouse modules: DAM (Keren-Shaul 2017), BAM, lymphatic, trafficking, ECM and others.
  - Technical sentinel sets: mitochondrial-encoded and cytosolic ribosomal genes.
  - All gene lists: `data/05_pathways_cell_state/gene_sets_long.csv`.
- **Pathway tests:**
  - Initial analyses used gseapy preranked GSEA on the Wald statistic. Its gene-permutation q-values are optimistic at n = 6.
  - Later analyses used limma `camera` on VST (sample-level, accounts for inter-gene correlation), plus per-sample module scores (mean z of VST) tested with OLS and exact within-sex label permutation.
- **Single-cell references:** CELLxGENE Census release 2025-01-30, read directly from its public S3 bucket. Pseudobulk = per-cell CP10K, averaged per dataset, then across datasets; up to 1,200 cells per class.
  - **33-class atlas reference** (`reference_atlas_33class_*`): used for the fibroblast analysis.
  - **20-class dura-centric reference** (`reference_dura_20class_*`): immune/stromal/endothelial classes from the adult mouse dura dataset 58b01044-c5e5-4b0f-8a2d-6ebf951e01ff; lymphatic and mural cells from other adult tissues.
  - In both references, osteoclast, choroid plexus and olfactory classes come from an embryonic atlas, because no adult data existed.
- **Software:** Python 3.11 (pandas 2.3, numpy 2.4, scipy 1.17, statsmodels 0.15, gseapy 1.3, pydeseq2 0.5.4, tiledbsoma 2.3); R 4.3.3 (DESeq2 1.42, limma 3.58, MatchIt 4.5.5, car 3.1). Random seed 7.

## 4. Initial whole-tissue analysis (`data/02_initial_whole_tissue/`, figures A01–A05)
**Methods:**
- DESeq2 `~ sex + group`.
- An "adjusted" model adds composition covariates: osteoclast, choroid plexus and olfactory marker modules, plus % mito.
- Preranked GSEA.

**Results:**
- **What the tissue is:** a dura/meninges-on-skull preparation, not brain. Top genes are Ttr, Col1a1/2, Mgp, Ptgds, Dcn and hemoglobins. There is variable carry-over of bone/osteoclast, olfactory/nasal mucosa, choroid plexus, blood and muscle.
- **What drives the transcriptome PCs:**
  - Transcriptome PC1 (40%) is olfactory/nasal carry-over (Cyp2a5, Cyp2g1, Scgb1c1). Samples 0453 and 613 are extreme.
  - PC2 (23%) is bone plus library profile.
- **5xFAD_NTC vs WT_NTC:** 6,811 genes at padj < 0.05.
  - The "up" genes are dominated by pseudogenes (253), ribosomal genes (167) and mitochondrial genes (12).
  - Adjusting for % mito alone leaves 1 gene.
  - The signature correlates with % mito (r = 0.94) and % multi-mapped reads (r = 0.90) (`genotype_DE_sensitivity.csv`, `genotype_signature_gene_classes.csv`; figure A05).
- **Bulk Spp1 tracks osteoclast content** (Spearman ρ = 0.60; top correlates Atp6v0d2, Ctsk, Acp5, Mmp13). Spp1 is not significantly lower in SPP1 animals (log2FC −1.2, padj 0.67); the drop is seen only in SPP1 males (figure A03).
- **SPP1 vs NTC:** only 4 genes at padj < 0.05: Rgs4, Cldn11 and Prg4 down; Bpifa1 up. Bpifa1 is nasal, i.e. contamination.
- The GSEA "SPP1: IFN up, EMT/ECM down" result came from preranked GSEA and is optimistic; treat it as hypothesis only.

## 5. Fibroblast-specific analysis (`data/06_fibroblast_specific/`, figures B01–B04)
**Goal:** minimize contamination by analyzing only genes whose bulk signal comes from fibroblasts.

**Methods (final version):**
1. **Deconvolution against the 33-class reference:** NNLS on class marker genes (≥ 3–5× specific).
   - Mitochondrial, ribosomal and Ig/TCR-V genes are excluded as markers.
   - Weights are 1/√(expression + 1).
   - Iterative trimming drops markers the fit misses by more than 4–8× (these come from cell types absent from the reference, e.g. nasal glands).
   - The final estimate is the mean of 4 settings.
   - Output: `deconv_fibro_ensemble_proportions.csv` (31 classes) and `deconv_fibro_ensemble_grouped_categories.csv` (9 grouped categories, used for figure B03).
2. **Fibroblast signal share per gene:** the fraction of each gene's predicted bulk signal that comes from fibroblast classes.
   - Strict: median ≥ 0.8 and min ≥ 0.6. Relaxed: median ≥ 0.6 and min ≥ 0.4.
   - Genes whose bulk expression tracks a contaminant module (|r| ≥ 0.5) are dropped, and canonical contaminant markers (Bglap, Ctsk, Ttr, Omp, Hbb…) are blacklisted.
   - Result: **151 genes** (66 strict, 85 relaxed); `fibroblast_specific_genes.csv`. Strongly enriched for Hallmark EMT/ECM (FDR 6e-10).
3. **DESeq2 on those 151 genes only,** so size factors come from fibroblast genes ("fibroblast-normalized"). Two models: `~ sex + group` and `~ sex + pct_mito + group`.

**Results:**
- **Composition (figure B03):** fibroblasts are about 17–38% of marker signal. Osteoclast is high in 620 (46%), 0455 and 0452; these are the three highest-Spp1 samples. Blood is high in 623, 616, 0455 and 0452 (matches hemoglobin reads). Olfactory/nasal content is high in 0453 and 613. Values are marker-RNA shares, not cell counts.
- **5xFAD vs WT:** 57/151 genes at padj < 0.05, but **0** after the % mito adjustment.
- **SPP1 vs NTC:** Cldn11, Thbd and Fbln1 are down; Matn4 and Ifit3b are up (Section 10).
- **ECM score** (Fbln1, Col3a1, Mmp2, Lum, Col16a1): lower with SPP1 (−0.31 SD, p = 0.097). Within 5xFAD it tracks Spp1 (ρ = 0.57, p = 0.055). It weakens after adjusting for osteoblast content (`ecm_score_*.csv`; figure B04).
- **Deconvolution caveat:** estimates depend on settings and on the reference. Choroid plexus and macrophage estimates disagree with independent markers; olfactory, osteoclast and erythrocyte estimates agree well.

## 6. Composition-confounding analysis of genotype (`data/03_composition/`, `data/04_genotype_confounding/`, figures D_*)
This is a structured analysis of whether 5xFAD effects survive differences in dissection, composition and technical factors.
- **Composition measure:** 20 broad marker signatures, each the mean VST z-score of 5–8 canonical markers (`signature_scores_all20.csv`).
  - 16 pass a coherence check (≥ 3 markers detected and median item-rest r ≥ 0.4).
  - **Lymphatic endothelial, mural, T-cell and plasma-cell signatures fail** (their markers do not co-vary), so these lineages cannot be measured (`signature_validation.csv`).
- **Composition PCA** on the 16 validated signatures (`composition_PC_*`):
  - **Composition_PC1** (26.9%): dura purity, fibroblast and blood endothelial vs olfactory/CNS/muscle.
  - **Composition_PC2** (22.7%): myeloid/NK/mast.
  - **Composition_PC3** (16.8%): bone (osteoblast and osteoclast).
- **Composition explains the transcriptome's main axes** (R²: transcriptome PC1 0.90, PC2 0.83), about 55% of total variance (`txPC_variance_explained_by_composition.csv`).
- **Genotype vs composition** (`genotype_vs_composition.csv`):
  - 5xFAD has lower Composition_PC1 (Hedges g −1.62, AUC 0.08, p = 0.015) and higher Composition_PC3 (g +1.63, AUC 0.94, p = 0.009).
  - **Common support between genotypes is poor** (figure D_Fig06).
- **Technical axis:** median % mito is 5.2 (WT) vs 10.0 (5xFAD); % multi-mapped is 10.3 vs 16.8. All 12 5xFAD libraries plus WT 616 share the high profile; 623 is intermediate on multi-mapping only.
- **Adjusted DESeq2** (`adjusted_models_summary.csv`; baseline 6,787 genes at padj < 0.05):

  | Model | Genes padj < 0.05 | Median effect retained |
  |---|---|---|
  | + Composition_PC1 | 1,678 | 72% |
  | + PC1 + PC3 | 0 | 53% |
  | + % mito | 1 | 32% |

  Variance inflation for genotype is 1.5–2.0 (`adjusted_models_collinearity.csv`).
- **Matched cohort:** MatchIt 1:1, Mahalanobis on Composition_PC1 + PC3, caliper 1 SD. It yields only 2 vs 2 (WT 616, 626 vs 5xFAD 617, 451).
- **Decisive single-animal test** (`technical_axis_animal_position_*`; figure D_Fig10c): on the 6,787 baseline genotype genes, each WT animal is placed between the other WT animals (0) and 5xFAD (1).
  - **WT 616's median position is 0.96; it looks 5xFAD on 84% of genes.** Other WT animals: 2–5%.
  - Either the "genotype" signature is technical/batch, or 616 is mis-genotyped. Either way, genotype is not separable in this dataset.
- **Pathways** (`pathway_robustness_matrix.csv`, figures D_Fig08 and D_Fig12):
  - Camera finds no set at FDR < 0.9.
  - The strongest "genotype" modules are the technical sentinels (ribosome, mito, OXPHOS; +1.4 to 1.5 SD).
  - **No Tier 1 or Tier 2 pathway.** Fibroblast activation is Tier 3 (bone axis); angiogenesis-endothelial and ROS are Tier 4 (technical).
  - Suggestive only: lipid handling (Fabp5-driven) and a DAM-like trend (incoherent directions).
- **Cell state vs abundance** (`cell_state_vs_abundance.csv`, figure D_Fig09): no robust state change beyond abundance.
- Full write-up: `composition_confounding_analysis/ANALYSIS_README.md` in the source repo (summarized here).

## 7. SPP1 vs the 5xFAD signature in fibroblast genes (`data/07_spp1_vs_5xfad_signature/`, figures C01–C03)
- **5xFAD signature, defined without bias:** WT vs all 12 5xFAD animals, within the 151 fibroblast genes. 60 genes at padj < 0.05, 0 after % mito adjustment. WT 616 position: median 0.99.
- **SPP1 shift along that signature:** 3–6% of the 5xFAD–WT gap (exact within-sex permutation over 400 relabelings, p ≈ 0.66–0.72). Gene-wide concordance ρ ≈ 0.00.
- **No evidence that SPP1 changes 5xFAD-associated fibroblast changes.**
- An earlier version that defined the signature as WT vs NTC showed an apparent 55% "rescue" (p = 0.01). That was **regression-to-the-mean selection bias** and is withdrawn.

## 8. TGF-β signaling panel (`data/08_tgfb_panel/`, figures E_*)
**Panel and weights:**
- High (weight 1): Serpine1, Smad7, Tgfbi, Pmepa1.
- Medium (0.5): Skil, Ccn2 (= Ctgf), Fn1.
- Low (0.25): Smurf2.
- Context genes: Tgfb1, Tgfb2, Tgfb3, Tgfbr1, Tgfbr2, Ltbp1.

**Results:**
- **Not one pathway.** The readouts co-vary (mean r = 0.49 vs 0.01 for random gene sets), but in two composition-linked modules:
  - Smad7, Tgfbi, Pmepa1 track fibroblast/endothelial content (r 0.7–0.9).
  - Serpine1, Ccn2, Skil, Smurf2 track bone/mural content.
- **5xFAD vs WT targets:** the high-confidence score is lower (−1.19 SD, p = 0.035), driven by Smad7 (46 → 22 CPM) and Pmepa1 (266 → 113 CPM). This **disappears** with composition adjustment (−0.27) or the % mito covariate (+0.28), and WT 616 looks 5xFAD on these genes. The composition-residual score shows no difference (p = 0.62). **Tier 3/4.**
- **Tgfb2 ligand: Tier 2 candidate.** About 35% lower in 5xFAD (75 → 37 CPM; VST ranges do not overlap: WT 8.99–9.40, 5xFAD_NTC 8.36–8.79). It is robust to composition, % mito (p = 0.006), ID series, and all covariates jointly (p = 0.018). **WT 616 is WT-like (position 0.12).** It is a single gene and the source cell type is unclear. Tgfb1 behaves similarly but is weaker.
- **SPP1:** no change in pathway activity (score permutation p = 0.11–0.32). **Fn1 is lower** (−1.39 SD, permutation p = 0.01, composition-robust).

## 9. Genotype verification attempt
- App and Psen1 counts are *lower* in 5xFAD, but WT 616 is equally low, so this is technical.
- Thy1 (the transgene cassette contains Thy1 exons) is not different.
- Neuronal content is too low for the Thy1-driven transgenes to show.
- **Verifying genotype requires FASTQ/BAM files or genotyping records.**

## 10. SPP1 vs NTC "matrix" changes, consolidated (`data/09_spp1_matrix_genes/`)
| Gene | Role | log2FC (fibroblast-normalized) | CPM NTC → SPP1 | Evidence |
|---|---|---|---|---|
| Cldn11 | Arachnoid barrier tight junction | −0.34 | 204 → 160 | padj 0.048 (0.039 mito-adjusted); whole tissue padj 0.001; both sexes |
| Fbln1 | Fibulin-1 (basement membrane) | −0.33 | 286 → 224 | padj 0.039 (mito-adjusted); both sexes |
| Thbd | Thrombomodulin (endothelial/meningeal) | −0.36 | 223 → 181 | padj 0.048 / 0.039; stronger in males |
| Fn1 | Fibronectin | whole tissue −0.45 | 336 → 264 | permutation p = 0.01; composition-robust; mostly males |
| Matn4 | Matrilin-4 | +0.55 | 11 → 15 | padj 0.048 / 0.039 |
| Ifit3b | Interferon-induced | +0.56 | 17 → 23 | padj 0.001 (mito-adjusted) |
| Rgs4, Prg4 | — | whole tissue −0.91 / −0.65 | — | padj 0.001 (about 0.056 after composition adjustment) |

- **Update (Section 0, genome-wide R analysis):** Cldn11, Prg4 and Rgs4 are the only genes significant in ≥ 4 of 5 models. Fn1, Fbln1, Thbd, Loxl3 and Col12a1 are significant only with the % mito covariate, and lose about half their effect after adjusting for the lower fibroblast signal in SPP1 samples.
- Col3a1, Lum, Mmp2, Col16a1 and Loxl3 all trend down but are **not individually significant**.
- Effects are small (about 20%), n = 6 per group, often stronger in males.
- SPP1 knockdown itself is not detectable in this tissue.
- The SPP1 group has slightly less osteoblast signal, which partly confounds the matrix trend.
- In the summary table, `permutation_p_within_sex` was computed only for genes in the TGF-β analysis (Fn1); it is NaN for the others.

## 11. Things that could NOT be assessed
- Lymphatic endothelium, T cells, plasma cells and mural cells (no coherent markers).
- Any pathology or severity association (no data).
- Age effects.
- Batch, as a recorded variable.
- Transgene presence.

## 12. History of revisions (important for correct interpretation)
1. The initial 5xFAD-vs-WT signature (6,811 genes) was later shown to be technical/composition-driven.
2. **Withdrawn:** "SPP1 lowers the fibroblast ECM/EMT program (FDR 0.006)". This came from a deconvolution that (a) used mitochondrial genes as fibroblast markers and (b) misassigned nasal gland genes (Scgb1c1, Agr2) to plasma cells. After the fix the gene set shrank from 443 to 151 genes and the program-level effect fell to p ≈ 0.1. The gene-level hits survived.
3. **Withdrawn:** "SPP1 rescues 55% of the 5xFAD fibroblast signature (p = 0.01)". This was selection bias. The corrected estimate is 3–6%, p ≈ 0.7.
4. An earlier "Spp1 tracks osteoclast content" analysis and the deconvolution both agree that bulk Spp1 reflects bone/osteoclast cells.

## 13. Overall conclusions (authoritative)
**High confidence**
- **Composition dominates the variation:** dura purity vs olfactory/CNS/muscle carry-over, and bone.
- **Genotype is confounded with library processing** (WT 616 test) and with composition, with poor common support.
- **The genome-wide 5xFAD signature is not interpretable as disease biology.**
- **SPP1 does not shift the 5xFAD signature.**

**SPP1 vs NTC (primary comparison; clean of the technical confound)**
- **Robust:** Cldn11, Prg4 and Rgs4 down; Ifit3b and Matn4 up; type-I interferon module up after fibroblast adjustment.
- **Conditional:** matrix / TGF-β / fibroblast-activation decreases (Fn1, Fbln1, Thbd, Loxl3, Col12a1, Wwtr1), accompanied by lower fibroblast signal in SPP1 samples. These are either a real reduction in fibroblast content or activity, or dissection.

**Suggestive / candidate**
- Tgfb2 lower in 5xFAD (Tier 2 candidate).

**Cannot determine**
- Cell-intrinsic 5xFAD effects.
- Lymphatic or T-cell biology.
- Pathology links.
- Whether SPP1 knockdown occurred in target cells.

## 14. Recommended next steps
- **Genotype WT 616** and check processing records (RNA extraction and library batches, RIN).
- **Re-process** WT and 5xFAD interleaved in the same batches.
- **Standardize dissection** (peel dura from the calvaria; exclude bone and the olfactory/cribriform region).
- **Immunostaining of dura whole-mounts:**
  - CLDN11 with E-cadherin (arachnoid barrier);
  - fibronectin and FBLN1 with PDGFRα;
  - THBD with CD31;
  - SPP1 (knockdown check);
  - TGF-β2 with CD31/PDGFRα;
  - pSMAD2/3;
  - Aβ (6E10);
  - composition markers: CD31, LYVE1/PROX1, PDGFRα, CD206, CD3.
  - Image matched regions (e.g. along the sinuses), score blinded, include both sexes, quantify per area.
- Flow cytometry or scRNA-seq of dura; record age and pathology; replicate in an independent cohort.

## 15. Plotting guide: figure → data
Figures are in `figures/`. The `group` palette used throughout: WT_NTC `#2a78d6` (blue), 5xFAD_NTC `#eb6834` (orange), 5xFAD_SPP1 `#1baf7a` (green). Circle = female, triangle = male. Diverging colormap: blue → grey `#f0efec` → red.

| Figure | What it shows | Data |
|---|---|---|
| A01 | Initial PCA (PC1–PC4) | `00_core/sample_master_table.csv` (Tx_PC1–5 come from the R analysis; initial-analysis PCs are reproducible from VST) + `01_qc_pca/transcriptome_PC_variance.csv` |
| A02 | Marker-module composition heatmap | master table, columns `initial_marker_module__*` (z-score across samples) |
| A03 | Spp1 vs osteoclast module; Spp1 by group | master: `Spp1_cpm` (log2(x+1)), `initial_marker_module__osteoclast` |
| A04 | Initial GSEA bar charts | `02_initial_whole_tissue/GSEA_*.csv` (Term, NES, FDR q-val) |
| A05 | Genotype signature score vs % multi / % mito | master: `genotype_signature_score_initial`, `pct_multi`, `pct_mito` |
| User heatmap | Sample-sample correlation | `01_qc_pca/sample_correlation_pearson_top2000_vst.csv` (hierarchical clustering) |
| B01 | Reference expression of the 40 most-expressed strict fibroblast genes | `06_fibroblast_specific/fibroblast_specific_genes.csv` (tier == strict, sort by fib_max_cp10k) + `reference_atlas_33class_pseudobulk_cp10k.csv.gz` (log10(CP10K + 0.01)) |
| B02 | Deconvolution heatmap and fibroblast fraction by group | `06_fibroblast_specific/deconv_fibro_ensemble_proportions.csv`; master `fibroblast_fraction_deconv` |
| B03 | **Per-sample cell representation (stacked bars)** | `06_fibroblast_specific/deconv_fibro_ensemble_grouped_categories.csv` (× 100); fit ρ from `deconv_fibro_ensemble_fit_quality.csv` |
| B04 | SPP1 ECM score by group/sex and vs Spp1 | `06_fibroblast_specific/ecm_score_per_sample.csv` |
| C01 | WT 616 position histogram; adjusted vs unadjusted LFC (fibroblast genes) | `07_spp1_vs_5xfad_signature/fibroblast_genes_genotype_and_spp1.csv` (t616, LFC columns), `t_leave_one_out_WT.csv` |
| C02 | Top genotype genes per animal | `06_fibroblast_specific/fibroblast_normalised_vst.csv` + gene list from C01 data (sort by padj_5xFAD_vs_WT) |
| C03 | SPP1 on the 5xFAD signature | `07_spp1_vs_5xfad_signature/signature_scores_per_sample.csv`, `spp1_signature_permutation_tests.csv`, `spp1_vs_5xFAD_concordance.csv` |
| D_Fig01 | Transcriptome PCA by genotype | master: `Tx_PC1`, `Tx_PC2` + `01_qc_pca/transcriptome_PC_variance.csv` |
| D_Fig02/02b | Signature heatmap / per-marker heatmap | `03_composition/signature_scores_all20.csv`; markers in `signature_marker_genes_item_rest.csv`, values from `00_core/vst_blind_filtered_genes.csv.gz` (row z-score) |
| D_Fig03 | Signature scores by group | `03_composition/signature_scores_all20.csv` (+ `signature_validation.csv` for which are validated) |
| D_Fig04 | Composition PCA + loadings | master: `Composition_PC1–3`; `03_composition/composition_PC_loadings.csv`, `composition_PC_variance.csv` |
| D_Fig05 | Transcriptome PCs vs composition | `03_composition/txPC_vs_composition_spearman_rho.csv` / `_p.csv` |
| D_Fig06 | Common support strip plots | master + `03_composition/genotype_vs_composition.csv` (AUC, common_support) |
| D_Fig07 | Baseline vs adjusted genotype LFC | `04_genotype_confounding/genotype_effects_baseline_vs_adjusted.csv.gz` (`LFC_baseline` vs `LFC_<model>`, `class_<model>`) |
| D_Fig08 / D_Fig12 | Pathway effects by model; robustness tiers | `05_pathways_cell_state/pathway_robustness_matrix.csv` |
| D_Fig09 | State vs abundance | `05_pathways_cell_state/cell_state_scores_per_sample.csv` + `03_composition/signature_scores_all20.csv` + `cell_state_vs_abundance.csv` |
| D_Fig10a/b | Matched cohort selection; cohort LFCs | `04_genotype_confounding/common_support_cohorts_animals.csv`; `DE_matched_*.csv.gz`, `DE_tech_profile_B_*.csv.gz` vs `DE_baseline_*` |
| D_Fig10c | **WT-animal position test** | `04_genotype_confounding/technical_axis_animal_position_per_gene.csv.gz` (histogram per WT animal; clip −1.5 to 2.5) |
| E_TGFb_A | Panel genes per animal | `00_core/vst_blind_filtered_genes.csv.gz` rows for the panel genes; p-values from `08_tgfb_panel/tgfb_panel_effects.csv` |
| E_TGFb_B | Effect heatmap across models | `08_tgfb_panel/tgfb_panel_effects.csv` (`*_SD` columns; `*_p` for asterisks) |
| E_TGFb_C | Panel scores per animal | `08_tgfb_panel/tgfb_panel_scores.csv` |
| E_TGFb_D | Gene-gene correlations; cell origin | `08_tgfb_panel/tgfb_panel_gene_correlations.csv`; `tgfb_panel_reference_expression_cp10k.csv` (log10) |
| — | SPP1 matrix-gene table | `09_spp1_matrix_genes/spp1_vs_ntc_matrix_gene_summary.csv` |
| **F_SPP1_A** | **SPP1 vs NTC composition balance** | master table (5xFAD rows) + `10_SPP1_vs_NTC_PRIMARY/SPP1_vs_NTC_composition_balance.csv` |
| **F_SPP1_B** | **Key SPP1 genes per animal** | `00_core/vst_blind_filtered_genes.csv.gz` (5xFAD columns) + `10_SPP1_vs_NTC_PRIMARY/SPP1_vs_NTC_key_genes.csv` |
| **F_SPP1_C** | **Module effects under 4 adjustments** | `10_SPP1_vs_NTC_PRIMARY/SPP1_vs_NTC_module_effects.csv` (`SPP1_effect_SD*`, `perm_p_within_sex`) |
| **F_volcano / F_MA** | **Genome-wide SPP1 vs NTC (5xFAD-only model)** | `10_SPP1_vs_NTC_PRIMARY/DE_SPP1_vs_NTC_5xFADonly.csv.gz` |

## 16. Key data dictionary
**`00_core/sample_master_table.csv`** has one row per animal (18 rows) and is the main file for per-sample plots:
- **Identifiers:** `library_id`, `sample`, `animal`, `genotype` (WT/5xFAD), `sex` (F/M), `treatment` (NTC/SPP1), `group`, `id_series` (04xx/6xx, possible cohort).
- **QC:** `total_reads`, `pct_unique`, `pct_multi`, `assigned_counts`, `genes_detected`, `pct_mito`, `pct_hemoglobin`, `sex_check`, `tech_profile` (A_low_multimap / B_high_multimap), `size_factor`.
- **Transcriptome PCs:** `Tx_PC1–5` (R analysis; top 2,000 VST genes).
- **Composition:** `Composition_PC1–3`.
- **Spp1:** `Spp1_cpm`.
- **Signatures and modules:**
  - `signature__<lineage>`: marker-signature scores (mean z of VST); see validation for which are reliable.
  - `initial_marker_module__<x>`: initial-analysis modules (mean log2 CPM of markers).
- **Deconvolution:**
  - `deconv_dura_ref__<category>`: 20-class dura-reference deconvolution, grouped (fraction 0–1; unreliable for immune classes).
  - `deconv_fibro_ensemble__<category>`: 33-class ensemble deconvolution, grouped (fraction 0–1; figure B03).
- **Other scores:**
  - `fibroblast_score_strict_genes`, `fibroblast_fraction_deconv`;
  - `fibroblast_ECM_score`: SPP1 ECM leading-edge score;
  - `genotype_signature_score_initial`: initial 5xFAD-vs-WT signature score;
  - `TGFb_weighted_all8`, `TGFb_high_conf_4`, `TGFb_weighted_composition_residual`.

**DE tables:**
- `gene_id`, `gene_name`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`.
- Fibroblast DE adds `tier`, `fib_class` and `specificity_ratio`.

**`04_genotype_confounding/genotype_effects_baseline_vs_adjusted.csv.gz`:**
- For each model in {comp1, comp12, comp13, tech, comp1tech}: `LFC_<model>`, `padj_<model>`, `retention_<model>` (adjusted / baseline LFC) and `class_<model>`.
- Model formulas: comp1 = +Composition_PC1; comp12 = +PC1 + PC2; comp13 = +PC1 + PC3; tech = +pct_mito; comp1tech = +PC1 + pct_mito. All also include `~ sex + group`.

**`05_pathways_cell_state/pathway_robustness_matrix.csv`:**
- `<model>_effect_SD`, `<model>_p` (OLS on module score, n = 12 NTC animals);
- `*_retention`;
- `matched_effect_SD` (2 vs 2);
- `position_616` (0 = WT-like, 1 = 5xFAD-like; unstable when the baseline effect is about 0);
- `camera_FDR_*`;
- `confidence_tier`.

**`08_tgfb_panel/tgfb_panel_effects.csv`:**
- `geno_*` effects (SD units, WT_NTC vs 5xFAD_NTC) under `base`, `compPC1`, `compPC1PC3`, `mito` and `fibroblast` models;
- `position_616`;
- `spp1_*` effects (5xFAD_SPP1 vs 5xFAD_NTC);
- `spp1_perm_p` (exact within-sex permutation).

**Reference pseudobulks:** genes × cell classes, mean CP10K. Class composition and cell counts are in the `*_summary.csv` files.

## 17. Conventions and pitfalls for re-analysis
- **Fold-change direction:** log2FC is always first group minus second. Effects in "SD units" are standardized within the animals in that comparison.
- **The primary genotype contrast is 5xFAD_NTC vs WT_NTC.** SPP1 animals were included in DESeq2 models only to help estimate dispersion.
- **Treat any 5xFAD-vs-WT result as confounded** unless it survives the % mito adjustment **and** WT 616 looks WT-like.
- **Deconvolution numbers are RNA shares, not cell counts.**
- **Gene-permutation GSEA q-values are optimistic;** prefer camera or sample-permutation results.
- **Do not re-derive a 5xFAD signature from WT vs NTC and then test SPP1 on it** (selection bias, Section 12.3).

**`10_SPP1_vs_NTC_PRIMARY/` columns:**
- `SPP1_vs_NTC_all_genes_all_models.csv.gz`: `LFC_<model>`, `p_<model>` and `padj_<model>` for models baseline, comp1, comp13, tech and 5xFADonly; plus `min_padj_any_model`, `n_models_padj05` and `sign_consistent`.
- `SPP1_vs_NTC_key_genes.csv`:
  - `SPP1_effect_SD` (VST, standardized across the 12 5xFAD animals; `~ sex`);
  - `_fibroblast_adj` and `_compPC1PC3` variants, each with its p;
  - DESeq2 padj per model;
  - `VST_diff_females` / `VST_diff_males` (median SPP1 − median NTC);
  - `perm_p_within_sex` (exact, 400 relabelings).
- `SPP1_vs_NTC_module_effects.csv`: module effects with the same adjustments, plus `camera_p` and `camera_FDR`.
- `SPP1_vs_NTC_composition_balance.csv`: SPP1 − NTC difference (SD), `perm_p`, `AUC_SPP1_gt_NTC` per composition or technical measure.
