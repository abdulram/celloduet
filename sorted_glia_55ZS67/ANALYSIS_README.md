# 55ZS67: bulk RNA-seq of sorted CD11b+ and ACSA-2+ cells, SPP1 vs NTC

The brief from the user is in `BRIEF_FROM_USER.md`. The one-page results are in `EXECUTIVE_SUMMARY.md`. Per-finding tiers are in `tables/analysis_summary.tsv`.

## Inputs (`data/`)
- **Count and CPM matrix:** `55ZS67-expression-matrix.tsv.gz`, from the sequencing provider, Ensembl gene IDs.
  - Counts are fractional (about 9% of entries), consistent with multimappers being split between loci. They are rounded for DESeq2.
- **Mapping stats** (unique / multi / unmapped reads) and **biotype summaries:** provider files.
- **ERCC spike-ins** are in the annotation, but they total only about 100 reads per library. They are uninformative and excluded.
- **Not available:** FASTQ, BAM, batch information, RIN, sort records.

## Design (inferred from sample names; see open questions)
| sample | facility name | fraction | treatment | replicate |
|---|---|---|---|---|
| 55ZS67_1–3 | rerun_NTC_microglia_1–3 | CD11b | NTC | 1–3 |
| 55ZS67_4–6 | rerun_Spp1_microglia_1–3 | CD11b | SPP1 | 1–3 |
| 55ZS67_7–9 | rerun_NTC_astrocyte_1–3 | ACSA2 | NTC | 1–3 |
| 55ZS67_10–12 | rerun_Spp1_astrocyte_1–3 | ACSA2 | SPP1 | 1–3 |

- **Sex:** all male (Xist < 1 CPM; chrY genes 13–130 CPM).
- **Genotype:** all 5xFAD (confirmed by the user). Age not stated. No WT arm.
- **Sorting:** CD11b+ depleted before ACSA-2 enrichment; warm dissociation; whole brain. The libraries are a re-run after poor RNA quality in the first attempt.

## Workflow (run in order from `sorted_glia_55ZS67/`)
| script | what it does |
|---|---|
| `00_gene_sets.py` | Writes curated lineage panels and state modules to `data/modules.gmt` (brief §5.3/5.7). |
| `01_intake_qc.py` | Sample sheet, mapping/biotype/mito/ribo metrics, sex calls → `01_QC/`, `tables/sample_sheet.tsv`. |
| `02a_marker_table.py` | Raw CPM of lineage, state and Spp1-axis genes. |
| `02_purity_artifacts.py` | Panel scores and validation, cross-fraction carry-over, IEG dissociation score, and exact permutation (20 relabelings) of every metric vs treatment within fraction → `02_purity_celltyping/`, `03_artifacts/`. |
| `04_deseq2_per_fraction.R` | Per fraction `~ treatment` DESeq2 (Wald; normal-shrunk LFC). One-covariate sensitivity models. Leave-one-out refits. VST. limma-voom `camera` on curated modules and MSigDB Hallmark (human → mouse by symbol case). |
| `05_pca_variance.py` | PCA per fraction and PC–covariate R². |
| `06_carryover_explainability.py` | For each gene, the log2FC that carry-over alone would produce (intrinsic + a·other-fraction profile, with a = measured carry-over share). Flags DE genes it explains. |
| `07_pairing_test.py` | Tests whether replicate *k* is the same animal in both fractions (residual-profile correlation). |
| `08a_build_brain_reference.py`, `08b_deconvolution.py` | CELLxGENE Census (2025-01-30) adult-brain pseudobulk with 15 classes (no choroid-plexus class was available). NNLS with trimming on specific markers (mt, ribosomal, Ig/TCR-V and Hb genes excluded), 20/30/50 markers per class, median reported. |
| `09_figures_key_genes.py` | Per-sample plots. |
| `10_pathology_programs.py` | 20 literature 5xFAD pathology programs. Spp1 removed from all sets; dissociation genes excluded. Per-animal scores; Welch CI and exact permutation; minimal detectable difference; Spp1/DAM knockdown index; ACSA-2 sets normalised to host-lineage content → `09_pathology/`, `figures/fig10`. |
| `11_rna_quality_checks.py` | Gene-length vs log2FC and a long/short-gene index per sample (degradation proxy without BAMs). |

**Software:** R 4.3 / DESeq2 1.42 / limma / edgeR (`logs/sessionInfo_R.txt`); Python 3 with pandas, numpy, scipy and tiledbsoma. Seeds are fixed (1 and 7).

## Results

### 1. Technical QC (`01_QC/qc_metrics.tsv`, `03_artifacts/metric_vs_treatment_tests.tsv`)
- **All libraries are comparable overall:** 19–21 M reads, 82–89% uniquely mapped, 9–11% multi-mapped.
- **% mitochondrial differs by fraction:** CD11b 4.5–5.5%, ACSA2 10.7–12.7%. Genes detected (≥ 5 counts): 17–21k.
- **Assigned share of uniquely mapped reads is low** (42–66%), consistent with intronic/intergenic reads not being counted.
- **Confound: within each fraction, a library-quality axis is perfectly separated by treatment.**
  - CD11b: SPP1 has higher % unmapped (4.1 vs 2.8) and a lower assigned rate (51 vs 60%).
  - ACSA2: NTC has higher % unmapped (7.0 vs 4.1) and more genes detected.
  - With n = 3 the exact p is 0.1, which is the smallest possible.
  - The direction is opposite in the two fractions, so it is not a treatment biology effect. More likely each fraction × treatment was handled as a unit.
  - Processed pseudogenes and unnamed loci (about 20% of DE calls) follow this axis. They are treated as Tier 4.

### 2. Sort identity and purity
- **CD11b+:**
  - 95–99% microglial RNA. Microglia markers (P2ry12 ≈ 1,300, Hexb ≈ 12,000 CPM) are 50–100× higher than in ACSA2.
  - BAM 0.1–2.7%, monocyte < 0.6%, endothelial ~1%.
  - Sample 6 (SPP1) has visible neutrophil/blood RNA (S100a8/9, Ngp, Camp, Lcn2; still < 1% share).
  - Astrocyte/oligodendrocyte markers are about 1% of ACSA2 levels, slightly higher in SPP1 (complete separation).
  - Choroid plexus Ttr is 500–2,100 CPM in every sample. It is likely ambient (Ttr is extremely abundant in the CP), and the reference has no CP class.
- **ACSA-2+:**
  - Astrocyte markers present (Atp1b2/ACSA-2 antigen 700–1,100 CPM, Slc1a3, Aldoc, Gja1, Aqp4), **but the oligodendrocyte signal dominates** (Plp1 18,000–24,000 CPM; Mbp, Mog, Mobp).
  - Deconvolution: oligodendrocyte 43–63%, astrocyte 25–48%, ependymal 4–10%, microglia 0.4–1.8%.
  - The oligo/astro split varies by up to 20 points across marker settings. It is not stable at the decimal level, but every setting agrees oligodendrocyte RNA ≥ astrocyte RNA.
  - Mature-oligodendrocyte (or myelin-debris) RNA is a known ACSA-2 MACS contaminant.
- **Panel validation** (item-rest r ≥ 0.4, ≥ 3 genes) fails for low-level contaminant panels in CD11b (OPC, pericyte, fibroblast, B cell) and for neuron, CP and erythrocyte panels in ACSA2. Those panel scores are not interpreted.

### 3. Dissociation / ex vivo activation
- The IEG score is high in both fractions (microglia: Fos 2,000–4,100 CPM; Ccl3/Ccl4 1,700–4,300 CPM), indicating dissociation activation.
- It does not separate the arms (no complete separation; PC2 in CD11b, PC2 in ACSA2). It was used only as a sensitivity covariate.

### 4. Exploration (`04_exploration/`)
- **CD11b:** PC1 (38%) correlates with treatment (r = 0.95) and with carry-over scores (astro, pericyte, endothelial; r ≈ 0.9). The two cannot be separated at the PC level. PC2 (22%) is Ttr and IEG; PC3 (18%) is the mito/ribo/quality axis.
- **ACSA2:** PC1 (42%) is genes detected, treatment, and oligo/ependymal content. PC2 (27%) is IEG and astrocyte share.

### 5. DE: CD11b+, SPP1 vs NTC (`05_DE_CD11b/`)
- 15,493 genes tested; 70 at padj < 0.05 and 116 at padj < 0.1. Of the 116:
  - 61 are ≥ 4× higher in ACSA2 (carry-over);
  - 23 are pseudogenes or unnamed (quality axis);
  - about 15 are microglia-enriched and not explained by carry-over.
- **Spp1:** log2FC −0.81 (shrunk −0.52), padj 1.9e-7.
  - CPM: NTC 1,608 / 2,385 / 2,156 vs SPP1 1,187 / 1,024 / 1,311.
  - Leave-one-out max p = 8e-4.
  - Adjusted: for assigned rate −0.55, for astro carry-over −1.21, for IEG −0.68.
  - Among 1,497 genes with baseMean > 1,000, only Plp1 (carry-over) and Tsc22d3 (stress, n.s.) have larger |log2FC|.
- **Co-varying microglial genes (Tier 2):**

  | gene | log2FC | padj |
  |---|---|---|
  | Gpnmb | −0.74 | 0.10 |
  | Ch25h | −0.51 | 0.013 |
  | Fabp5 | −0.31 | 0.055 |
  | Rgs16 | −1.05 | 0.001 |
  | Chst1 | −0.54 | 0.005 |
  | Mrc2 | −0.41 | 0.013 |
  | Lat2 | −0.39 | 0.016 |
  | Hpse | −0.40 | 0.019 |
  | Mir142hg | −0.34 | 0.016 |
  | Lmo2 | −0.30 | 0.07 |
  | Mcoln3 | −1.40 | 0.05 |
  | Lgi2 | −0.63 | 0.09 |
  | Gpr165 (up) | +1.32 | 2e-4 |
  | Ccl7 (up) | +0.90 | 0.007 |

- **Unchanged:**
  - DAM core: Cst7, Itgax, Clec7a, Lpl, Cd9, Trem2, Tyrobp.
  - Homeostatic: P2ry12, Tmem119, Cx3cr1.
  - Spp1 receptors.
  - MHC-II (n.s., driven by sample 5), IFN, Tgfb1/2.
- **Camera (curated):**
  - Down: Lipid (FDR 1e-5), DAM2 (6e-5, includes Spp1), MGnD (1e-4), DAM2-without-Spp1 (1e-3), lysosome (2e-3).
  - Up: astrocyte, oligodendrocyte and BAM panels (carry-over).
- **Hallmark:**
  - EMT and myogenesis up (carry-over / vascular).
  - TNFα/NF-κB up. This matches the IEG/dissociation set, which is also up (FDR 1e-3) even though the summary IEG score was not separated. Treat it as Tier 3–4.

### 6. DE: ACSA-2+, SPP1 vs NTC (`06_DE_ACSA2/`)
- 17,529 genes tested; 16 at padj < 0.05 and 61 at padj < 0.1. Of the 61, 27 are pseudogenes or unnamed, and most of the rest are ependymal/cilia genes (down) or neuronal genes (Kif5a/c, Kcnc3, up).
- **Spp1:** log2FC −0.35 (n.s.); −0.47 after microglial carry-over adjustment (n.s.). Spp1/Hexb is unchanged.
- **Camera:**
  - Up: cholesterol synthesis and astrocyte panel (both consistent with the higher astrocyte share); IFN (A_IFN FDR 8e-4, M_IRM_IFN 2e-3).
  - Down: microglia, DAM and complement (lower microglial carry-over).
  - OXPHOS and fatty-acid Hallmark up, aligned with the composition and quality axes.

### 7. Glucocorticoid program (`figures/fig08`)
- Fkbp5, Tsc22d3, Sgk1, Klf9, Per1, Mt2, Cdkn1a and Sult1a1 are higher in SPP1 in both fractions, most strongly in SPP1 replicate 2 (samples 5 and 11).
- Core clock genes are flat, which argues against time of day.
- Classified as Tier 4 (harvest stress or a systemic treatment effect; not separable).

### 8. Pairing and cross-talk (`08_crosstalk/`)
- The matched replicate pairing gives the highest cross-fraction residual correlation in both arms (NTC mean r 0.27, SPP1 0.09; each ranks 1st of 6 pairings).
- Blood (Hb) and stress genes co-vary within the same replicate numbers.
- Formal ligand–receptor analysis was **not** run. The ACSA2 fraction is too mixed to attribute ligands to astrocytes, and pairing is not confirmed.

### 9. Pathology programs (round 2; `09_pathology/`)
**Gene sets.** Sets were fixed from the literature before looking:
- **Microglia:** DAM stage 1/2 (Keren-Shaul 2017), MGnD (Krasemann 2017), ARM and IRM (Sala Frigerio 2019), plaque-induced genes (Chen 2020, microglia-expressed subset), plus curated programs.
- **Astrocytes:** DAA (Habib 2020), pan-reactive (Zamanian 2012), A1 (Liddelow 2017).
- **Oligodendrocytes:** DOL (Kenigsbuch 2022, approximate list).
- Spp1 was removed from every set. Ch25h and Fabp5 were kept out of the lipid set because round 1 flagged them in these same animals (selection-bias rule).

**Microglia (CD11b):** no program changes beyond its 95% CI.

| program | change | 95% CI |
|---|---|---|
| DAM2 | −3% | −15 to +10 |
| PIG | −4% | −11 to +3 |
| MGnD | −10% | −23 to +6 |
| homeostatic | −1% | −11 to +11 |
| lysosome | −7% | −17 to +4 |
| proliferation | −5% | |
| IRM | 0% | |
| MHC-II | +8% | −28 to +63 (driven by sample 5) |

- Minimal detectable differences at 80% power: DAM2 ±17%, PIG ±10%, homeostatic ±14%.

**Spp1 knockdown index** (log2 Spp1 − mean log2 DAM2): NTC 0.96 / 1.68 / 1.37 vs SPP1 0.62 / 0.46 / 0.70. The groups separate completely. Spp1 output per unit of DAM response is roughly halved.

**ACSA-2:**
- Raw: DAA −3%, pan-reactive 0%, A1 +7%, DOL +1%, cholesterol +11%, IFN +18%.
- The astrocyte lineage score is +0.24 log2 in SPP1 (composition). Normalised per astrocyte: DAA −18%, pan-reactive −16%, A1 −10%, cholesterol −6%, IFN −0.1%.
- Normalised per oligodendrocyte: DOL −3%.
- Nothing is significant, and the detectable difference is about ±30%. **The round-1 ACSA-2 IFN signal is composition (Tier 3).**

**RNA quality.**
- CD11b: log2FC is unrelated to gene length (Spearman ρ = 0.02), and the long/short index is the same in both arms.
- ACSA-2: a mild length trend (ρ = −0.09; short-gene quintile +0.04 vs long −0.01 median log2FC), with the index the same in both arms.
- Degradation does not explain the Spp1 result.

**Interpretation.** SPP1 treatment produced a partial, specific Spp1 knockdown in microglia without a detectable change in the amyloid-response programs of microglia, astrocytes or oligodendrocytes. Transcriptomic proxies therefore give no evidence that plaque burden changed. Small decreases (≤ 20%) in the lipid/lysosomal DAM axis (Apoe, Gpnmb, Ctsb) cannot be excluded or confirmed.

## Deviations from the brief and reasons
- **No 5xFAD-vs-WT contrast, no sex term, no 400-relabeling permutation:** none is possible with this design (single genotype, all male, 3 vs 3 → 20 relabelings, minimum p = 0.1).
- **No transgene check:** no FASTQ, and no human APP/PSEN1 in the annotation.
- **No histology:** none provided. `09_pathology/` therefore holds transcriptomic pathology proxies only.
- **No choroid-plexus class in the deconvolution reference:** none in the Census adult brain data. Ttr-derived RNA is therefore unassigned.

## Failure conditions present
- Library-quality metrics separate treatment within each fraction.
- Purity (carry-over) differs by treatment in both fractions.
- n = 3 per arm.
- Genotype cannot be verified.
- The ACSA-2 fraction is not astrocyte-dominant.
