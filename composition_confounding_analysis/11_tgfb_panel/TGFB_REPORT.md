# TGF-β signaling across groups (WT_NTC, 5xFAD_NTC, 5xFAD_SPP1)

**Scripts:** `scripts/11a_tgfb_panel_deseq2.R` (per-gene DESeq2 under every model) and `scripts/11b_tgfb_panel_analysis.py` (cell origin, coherence, scores, permutation tests, figures).
**Tables:** `tables/tgfb_panel_*.tsv`. **Figures:** `11_tgfb_panel/TGFb_A`–`D`.

The analysis uses the same framework as the rest of `composition_confounding_analysis/`. Every genotype effect is checked against:
- **dura-purity composition** (Composition_PC1) and the **bone axis** (PC3);
- the **technical/library axis** (% mito);
- **fibroblast abundance**;
- **WT animal 616**, the only WT library with the 5xFAD-like technical profile. On the genome-wide "5xFAD signature", 616 looks 5xFAD-like on 84% of genes.

## Panel and weights
| Gene | Confidence | Weight | Main source in the dura single-cell reference | Bulk correlation with composition |
|---|---|---|---|---|
| Serpine1 | High | 1 | mural 35%, skeletal muscle 29%, osteoblast 14% | olfactory −0.62 |
| Smad7 | High | 1 | mast 22%, blood endothelial 18% | **endothelial 0.93, fibroblast 0.85** |
| Tgfbi | High | 1 | monocyte/macrophage 42% | **endothelial 0.86, fibroblast 0.84** |
| Pmepa1 | High | 1 | fibroblast 28%, osteoclast 20% | **fibroblast 0.74, endothelial 0.71** |
| Skil | Medium | 0.5 | mast 28% | osteoclast 0.61, olfactory −0.77 |
| Ccn2 (Ctgf) | Medium | 0.5 | mural 54% | olfactory −0.67 |
| Fn1 | Medium | 0.5 | monocyte 47%, osteoblast 15% | fibroblast 0.74, osteoblast 0.52 |
| Smurf2 | Low | 0.25 | osteoclast / choroid plexus / olfactory | osteoblast 0.60, osteoclast 0.58 |

Ctgf is annotated as **Ccn2** in this Ensembl release.
Context genes (pathway input, not readouts): Tgfb1, Tgfb2, Tgfb3, Tgfbr1, Tgfbr2, Ltbp1.

## 1. Do the readouts behave like one pathway? Partly, but the structure is composition
The mean pairwise correlation of the 8 readouts is r = 0.49, against r = 0.01 for random expression-matched gene sets (p = 0.002). The correlation matrix (`TGFb_D`) splits into two modules:
- **Module A: Smad7, Tgfbi, Pmepa1** (r = 0.69–0.81 with each other). These track the **fibroblast/endothelial (dura-purity) axis**.
- **Module B: Serpine1, Ccn2, Skil, Smurf2** (r = 0.54–0.84). These track **bone/osteoclast/mural content** and fall with olfactory carry-over.
- **Fn1** sits between the two.

So the panel co-varies mainly because the cells expressing it vary between dissections, not necessarily because TGF-β/SMAD activity varies.

## 2. 5xFAD vs WT
Effects are in SD units across the 12 NTC animals. 616 position: 0 = like the other WT animals, 1 = like 5xFAD.

| Readout | Baseline | + Comp PC1 | + PC1 + PC3 | + % mito | + fibroblast abundance | 616 position |
|---|---|---|---|---|---|---|
| High-confidence score (4 genes) | **−1.19 (p = 0.035)** | −0.27 | −0.35 | +0.28 | −0.47 | 1.59 |
| Weighted score (8 genes) | −0.85 (p = 0.16) | +0.09 | −0.22 | +0.65 | −0.07 | 2.09 |
| Composition-residual score | −0.31 (p = 0.62) | — | — | +0.89 | — | 2.97 |
| Smad7 | **−1.52 (p = 0.002)** | −0.68 | −0.49 | −0.33 | −1.03 | 1.07 |
| Pmepa1 | **−1.31 (p = 0.017)** | −0.85 | −1.00 | +0.15 | −0.76 | 1.30 |
| Tgfbi | −1.13 (p = 0.055) | 0.00 | +0.34 | +0.27 | −0.42 | 1.44 |
| Serpine1, Skil, Ccn2, Fn1, Smurf2 | −0.4 to +1.0 (all n.s.) | | | | | |
| **Tgfb2** (ligand) | **−1.73 (p < 0.001)** | **−1.69** | **−1.79** | **−1.66 (p = 0.006)** | **−1.68** | **0.12** |
| Tgfb1 (ligand) | −1.34 (p = 0.015) | −1.52 | −1.77 | −0.91 (p = 0.27) | −1.26 | 0.05 |
| Tgfb3 (ligand) | −1.34 (p = 0.013) | −0.65 | −0.57 | −0.22 | −0.77 | 0.89 |
| Tgfbr1 / Tgfbr2 / Ltbp1 | n.s. | | | | | |

**Readouts are lower in 5xFAD, but this is not robust (Tier 3/4).**
- At baseline, the high-confidence targets Smad7, Pmepa1 and Tgfbi are 30–60% lower in 5xFAD. For example, Pmepa1 is 266 vs 113 CPM and Smad7 is 46 vs 22 CPM.
- These genes track fibroblast/endothelial content. 5xFAD dissections contain less of both, and the effects shrink by 50–100% once composition is adjusted.
- They vanish with the technical covariate, and **WT 616 looks 5xFAD on all three** (positions 1.1–1.4).
- Once each gene's composition-driven part is removed, the panel score shows no genotype difference (−0.31 SD, p = 0.62).

**Ligand Tgfb2 is lower in 5xFAD, and robust to every check run here (Tier 2 candidate).**
- WT ranges 8.99–9.40 (VST) and 5xFAD_NTC ranges 8.36–8.79: no overlap, about 35% lower (DESeq2 log2FC −0.65; 75 vs 37 CPM).
- The effect is unchanged by composition adjustment (PC1, PC1 + PC3, fibroblast abundance), by % mito (p = 0.006), by the animal-ID series, and by all covariates fitted together (p = 0.018).
- **WT 616 sits with the WT animals (position 0.12).** This is unlike almost every other genotype-associated gene in the dataset.
- It passes Bonferroni across the 14 genes examined. It is not significant at genome-wide FDR once covariates are added (padj 0.48–0.96), because the genome-wide test has little power.
- Tgfb1 shows the same pattern (lower in 5xFAD, robust to composition, 616 WT-like at 0.05), but it weakens with % mito (p = 0.27) and is lowly expressed (2–5 CPM).
- The source cell is uncertain. In the reference, Tgfb2 is spread across olfactory, osteoclast and mural cells. In bulk it tracks endothelial content (r = 0.52) and runs opposite to bone (r = −0.53).

**Interpretation for 5xFAD.** There is no evidence for increased TGF-β pathway activity in 5xFAD dura. The apparent decrease in SMAD target genes is composition and technical. The one signal that survives is **lower TGF-β2 ligand expression** in 5xFAD. Since the downstream targets show no composition-independent change, this would point to reduced ligand supply without a detectable change in signaling output. That is a hypothesis, not a finding: it is a single gene, n = 6 per group, the source cell is unknown, and genotype can't be verified from this data.

## 3. SPP1 vs NTC (within 5xFAD; same library profile, so no technical confound)
| Readout | SPP1 − NTC (SD) | Permutation p (400, within sex) | + PC1 + PC3 | DESeq2 log2FC (padj) |
|---|---|---|---|---|
| Weighted score (8) | −0.87 | 0.11 | −0.40 | — |
| High-confidence score (4) | −0.61 | 0.32 | −0.08 | — |
| **Fn1** | **−1.39** | **0.010** | **−1.05** | −0.45 (genome-wide padj 0.27; 0.049 with % mito) |
| Skil | −0.96 | 0.085 | −0.71 | −0.17 |
| Smurf2 | −1.05 | 0.075 | −0.67 | — |
| Serpine1, Smad7, Tgfbi, Pmepa1, Ccn2 | −0.2 to −0.8 | 0.17–0.88 | — | — |
| Tgfb1 / Tgfb2 | +0.73 / +0.67 | 0.27 / 0.23 | — | — |

- SPP1 animals show **lower Fn1** (−21% CPM; strongest in males), robust to composition adjustment. This fits the earlier SPP1 matrix findings (Fbln1, Col3a1, Thbd down).
- Module-B readouts (Skil, Smurf2) trend down. The canonical SMAD-feedback targets (Smad7, Pmepa1, Tgfbi, Serpine1) do not change.
- Ligands move slightly *toward WT* (not significant): Tgfb2 median 8.70 in SPP1 vs 8.52 in NTC vs 9.11 in WT. That is at most a partial, non-significant shift.
- **No significant change in TGF-β pathway activity with SPP1.** There is a specific reduction in Fn1, an ECM / TGF-β-associated output, which is better read as part of the SPP1 matrix effect than as SMAD signaling.

## Bottom line
| Claim | Confidence |
|---|---|
| TGF-β/SMAD target activity differs between 5xFAD and WT | **Not supported**: composition/technical-driven (Tier 3/4) |
| Tgfb2 ligand expression is lower in 5xFAD | **Tier 2 candidate**: robust to composition, technical axis, ID series and the 616 test; single gene; genotype unverified |
| Tgfb1 lower in 5xFAD | Suggestive (composition-robust, weakened by the technical covariate) |
| SPP1 changes TGF-β pathway activity | **Not supported** (score permutation p = 0.11–0.32) |
| SPP1 lowers Fn1 | Suggestive-to-moderate (permutation p = 0.01; robust to composition; part of the broader SPP1 ECM effect) |

**Validation:** TGF-β2 immunostaining or ELISA in dura whole-mounts (WT vs 5xFAD); pSMAD2/3 staining to test pathway activity directly; Fn1 immunostaining for the SPP1 effect; qPCR of Tgfb2, Pmepa1 and Smad7 on cleanly peeled dura.
