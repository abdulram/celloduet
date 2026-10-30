# Brief for a new analysis agent: bulk RNA-seq of sorted CD11b+ and ACSA-2+ cells

You (the agent) are being handed a **new dataset**: bulk RNA-seq of cell populations sorted on **CD11b** (microglia, macrophages and other myeloid cells) and **ACSA-2** (astrocytes).

Your job is to **decipher what the data mean**. Verify the sorts, find and control artifacts, then extract robust biology. Report conservatively.

This brief contains:
- context from a previous, related project by the same user;
- the questions to ask before analyzing;
- a step-by-step analysis plan;
- the pitfalls that caused false positives last time;
- environment notes.

Read all of it before starting.

---

## 1. Context: the previous project (same lab, same research program)
- **Previous dataset:** bulk RNA-seq of **mouse dura**. Groups: WT_NTC, 5xFAD_NTC and 5xFAD_SPP1 (n = 6 each, 3F/3M).
  - **5xFAD** is an Alzheimer's amyloid model (human APP + PSEN1 transgenes, Thy1 promoter, neuronal).
  - **SPP1** = an osteopontin-targeting treatment. **NTC** = non-targeting control.
- **The user's main interest is 5xFAD_SPP1 vs 5xFAD_NTC.** Expect the same in this dataset. Also expect 5xFAD vs WT as a reference contrast.
- **Lessons from the dura analysis** (all of these happened; check for each one here):
  1. **Genotype was confounded with library processing.** All 5xFAD libraries had about 10% mitochondrial and 17% multi-mapped reads; 5 of 6 WT had about 5% and 10%. A naive 5xFAD-vs-WT test gave about 6,800 "DE genes" that were technical. One WT animal (616) processed like the 5xFAD libraries looked "5xFAD" on 84% of those genes. **Always test whether technical metrics separate your groups before interpreting DE.**
  2. **Tissue composition dominated the variance.** Dissection carry-over (bone, olfactory mucosa, blood, muscle) drove the top PCs. In sorted cells the equivalent is **sort purity / contamination** and **dissociation-induced activation**.
  3. **Bulk Spp1 did not show knockdown** in dura, because there Spp1 mostly came from osteoclasts. **In CD11b+ brain cells, Spp1 is a microglial disease-associated (DAM) gene, so this dataset may be the first real test of whether SPP1 treatment works.** Prioritize it.
  4. **Selection-bias trap:** a "5xFAD signature" built from WT vs NTC and then used to test SPP1 vs NTC gave a fake "55% rescue" (p = 0.01). It disappeared when the signature was built from WT vs all 5xFAD. Never select genes with the same animals you then test.
  5. **Deconvolution artifacts:** markers absent from the reference (nasal gland genes) were misassigned to plasma cells, and mitochondrial genes used as "markers" tracked the technical axis. Exclude mitochondrial, ribosomal and Ig/TCR-variable genes from markers. Check every deconvolved class against independent canonical markers.
  6. **Gene-permutation GSEA was optimistic** at n = 6. Prefer limma `camera` or sample-label permutation (exact within-sex relabeling).
  7. **Adjusting for a covariate that is itself a marker of the outcome is circular** (e.g. adjusting Fbln1 for a fibroblast score that contains Fbln1).
  8. **The prior dura SPP1 findings** to compare against: SPP1 lowered Cldn11, Prg4 and Rgs4 and raised Ifit3b, Matn4 and a type-I interferon module. Matrix/TGF-β decreases were partly explained by lower fibroblast content. Tgfb2 was lower in 5xFAD dura (candidate).
- **Previous outputs:** branch `claude/rnaseq-data-analysis-krbtjk` of repo `abdulram/celloduet`. See `chatgpt_handoff/README.md` and `composition_confounding_analysis/`. The scripts there (DESeq2 helpers, signature validation, permutation tests, CELLxGENE reference building) are reusable.

## 2. Terminology: "demultiplexing" and "cell typing" in BULK sorted data
- **No single cells exist in bulk data**, so there is no per-cell demultiplexing or clustering.
- **"Demultiplexing" can only mean sample demultiplexing:** splitting pooled sequencing reads into samples by index barcode (bcl2fastq / BCL Convert, usually done by the facility). If the user provides BCL files or undetermined reads, handle that first. Otherwise confirm the sample sheet maps each FASTQ/count column to the right animal and fraction.
  - **Exception:** if samples were pooled with hashtags or the user actually has single-cell data, stop and ask. That is a different workflow (HTO/genotype demultiplexing, e.g. via Seurat HTODemux or demuxlet/vireo).
- **"Cell typing" in bulk means three things:**
  - (a) **Sort identity:** does each sample express the expected lineage?
  - (b) **Purity and contamination:** what fraction of RNA comes from other cell types?
  - (c) **Subpopulation composition** within the sort: e.g. microglia vs border-associated macrophages vs monocytes within CD11b+.
- **Methods:** coherence-validated marker signatures (primary) and reference-based deconvolution (secondary). See Section 5.3.

## 3. Ask the user first (do not guess these)
1. **Tissue and region:** whole brain, cortex, hippocampus? Or dura/meninges again?
2. **Groups and n:** genotype (WT/5xFAD), treatment (SPP1/NTC), sex and age per animal.
   - Is there a WT_SPP1 group?
   - Are **both fractions from the same animals** (a paired design)?
3. **Sorting method:** MACS or FACS? Antibody clones? Was CD11b+ depleted before ACSA-2 enrichment? Viability gating? Any digestion (enzymatic, 37 °C) or cold mechanical dissociation? Were transcription/translation inhibitors (actinomycin D, triptolide, anisomycin) used during dissociation?
4. **Confirm "ASCA2" means ACSA-2** (astrocyte cell surface antigen-2, encoded by Atp1b2).
5. **Batches:** processing days, sort days, RNA extraction and library-prep batches, sequencing lanes, operator. Ask for RIN/DV200, cell numbers sorted and RNA input per sample. **Make a batch × group table immediately.**
6. **Library type and files:** full-length (e.g. SMART-seq) or 3' (QuantSeq)? Stranded? Raw FASTQ/BAM, or only a count matrix? Which genome and annotation?
7. **Genotyping records**, and whether FASTQs are available for transgene confirmation.
8. **Any pathology** per animal (plaque load, Aβ, Iba1/GFAP histology, behavior).
9. **The primary question.** Expected: SPP1 vs NTC within 5xFAD, per cell type. Also any hypotheses or genes of interest.

## 4. Environment notes (cloud container used previously)
- **Outbound network goes through a proxy.** Blocked last time: Enrichr, NCBI/GEO, Ensembl, Bioconductor, CRAN, CELLxGENE discovery.
- **Reachable:** PyPI, apt (run `apt-get update` first), GitHub raw, conda-forge/bioconda, and the CELLxGENE S3 bucket.
- **R:** `apt-get update && apt-get install -y --no-install-recommends r-base-core r-bioc-deseq2 r-bioc-limma r-cran-matchit r-cran-car r-cran-ggplot2 r-cran-pheatmap` (R 4.3, DESeq2 1.42).
- **Python:** `pip install pandas numpy scipy statsmodels matplotlib pydeseq2 gseapy cellxgene-census`. If pip fails on Debian-managed packages, add `--ignore-installed packaging`.
- **MSigDB Hallmark:** `https://raw.githubusercontent.com/zqfang/GSEApy/master/tests/extdata/h.all.v7.0.symbols.gmt` (human symbols; map to mouse by case-insensitive symbol match).
- **CELLxGENE Census without the discovery server:** reuse `analysis/reference/census_ctx.py` from the previous repo. It opens `s3://cellxgene-census-public-us-west-2/cell-census/2025-01-30/soma/` via tiledbsoma with the proxy settings. Mouse brain atlases there have microglia, BAMs, astrocytes, oligodendrocytes, OPCs, neurons, endothelium and pericytes, which is suitable for purity deconvolution.

## 5. Analysis plan
Put everything in a numbered, scripted workflow with seeds and session info. Analyze **each sorted fraction separately** for DE. Combine fractions only for QC, purity and crosstalk analyses.

### 5.1 Intake and sanity
- **Build the sample sheet:** sample ID, animal, fraction (CD11b / ACSA2), genotype, treatment, sex, age, batch fields, RIN, cells sorted.
- **If FASTQs exist:** sample demultiplexing check, FastQC/MultiQC, and alignment with STAR or salmon to GRCm39 + GENCODE vM-series. Also align, or use k-mer matching, against **human APP and PSEN1 cDNA** to confirm genotype; the 5xFAD transgenes carry Swedish/Florida/London APP and M146L/L286V PSEN1 mutations.
  - In sorted glia, transgene reads will be sparse (Thy1-driven, neuronal), but neuronal RNA contamination usually gives some.
  - A strong positive control for 5xFAD is the microglial **DAM program** (below). If 5xFAD CD11b+ samples lack DAM at ≥ 4 months, suspect genotype, age or sort issues.
- **Sex check:** Xist vs chrY genes (Ddx3y, Uty, Kdm5d, Eif2s3y).

### 5.2 Technical QC (the first big lesson)
- **Metrics:** library size, genes detected, % uniquely mapped, % multi-mapped, % mitochondrial, % rRNA, % intronic/intergenic, and gene-body 3'/5' bias (RNA degradation; needs BAMs).
- **Test whether any technical metric separates groups** (AUC, permutation). If a group is technically distinct, flag it prominently before any DE.
- **Single-animal anchor test:** if one animal of group A shares group B's technical profile, check whether it looks like B on the "group" genes. That is how the dura confound was proven.
- **Technical sentinel gene sets:** mitochondrial-encoded and cytosolic ribosomal genes. If these are your top "DE pathways", suspect technical signal.

### 5.3 Sort identity, purity and cell typing
**Marker panels.** Mean z of VST per panel. Validate each panel: ≥ 3 markers detected and median item-rest r ≥ 0.4.

| Population | Markers |
|---|---|
| Microglia (homeostatic) | P2ry12, Tmem119, Sall1, Cx3cr1, Hexb, Fcrls, Siglech, Gpr34, Olfml3 |
| Border-associated macrophages | Mrc1, Lyve1, Cd163, Pf4, F13a1, Ms4a7 |
| Monocytes | Ly6c2, Ccr2, Plac8, Chil3 |
| Neutrophils | S100a8, S100a9, Ly6g, Retnlg, Ngp, Camp |
| Dendritic cells | Flt3, Xcr1, Clec9a, Cd209a |
| T cells | Cd3e, Cd3g, Cd2 |
| B cells | Cd79a, Ms4a1, Cd19 |
| NK cells | Nkg7, Ncr1, Klrb1c |
| Astrocytes | Aldh1l1, Slc1a3, Slc1a2, Gja1, Aqp4, **Atp1b2 (ACSA-2 antigen)**, Sox9, Gfap, Aldoc |
| Neurons | Snap25, Rbfox3, Syt1, Stmn2 |
| Oligodendrocytes | Mbp, Plp1, Mog, Mobp |
| OPCs | Pdgfra, Cspg4, Olig2 |
| Endothelial | Cldn5, Pecam1, Flt1, Ly6c1 |
| Pericytes / SMC | Pdgfrb, Rgs5, Kcnj8, Acta2 |
| Fibroblasts / meningeal | Dcn, Col1a2, Slc47a1 |
| Ependymal / choroid plexus | Foxj1, Ttr |
| Erythrocytes | Hbb-bs, Hba-a1 |

Lineage-identity genes above should **not** be reused in state modules.

**Expected patterns:**
- CD11b+ from brain is about 90% or more microglia, with BAMs, monocytes and neutrophils as minor populations. It **contains no astrocyte RNA unless contaminated**.
- ACSA-2+ sorts often carry microglia, oligodendrocyte or endothelial contamination. Quantify **cross-fraction contamination**: microglial markers in ACSA-2 samples and astrocyte markers in CD11b samples.

**Deconvolution (secondary):**
- Build a pseudobulk reference from adult mouse brain atlases in the CELLxGENE Census, matching the region if possible.
- Fit NNLS on marker genes with iterative trimming. Exclude mitochondrial, ribosomal and Ig-V genes as markers.
- Report **RNA shares, not cell counts**. Validate each class against the marker panels and report stability across settings.

**Purity per sample:**
- Report purity and test whether it differs by group.
- If it does, assess common support (do groups overlap in purity?) before adjusting.

### 5.4 Dissociation / sorting artifact
- Enzymatic dissociation at 37 °C induces an **ex vivo activation signature** in microglia (Marsh et al. 2022 Nat Neurosci; van Hove/van den Brink 2017): Fos, Fosb, Jun, Junb, Egr1, Egr2, Atf3, Ier2, Ier3, Dusp1, Hspa1a, Hspa1b, Klf2, Klf4, Rgs1, Ccl3, Ccl4, Nfkbia, Zfp36.
- Astrocytes show similar immediate-early and heat-shock responses.
- **Score it per sample** and test whether it differs by group (batch or sort-day confound?).
- Consider it as a covariate only if it is not collinear with group, and report results with and without it.

### 5.5 Exploration
- **Per fraction:** VST, PCA, and sample-sample correlation.
- **Correlate PCs** with group, sex, technical metrics, purity scores, dissociation score and batch.
- **Quantify variance explained** by each (R²). In dura, composition explained about 55% of variance; find out what dominates here.

### 5.6 Differential expression (per fraction)
- **Primary contrast:** `5xFAD_SPP1 vs 5xFAD_NTC`, using `~ sex + treatment` within 5xFAD, or `~ sex + group` across all samples for dispersion.
- **Reference contrast:** `5xFAD_NTC vs WT_NTC`, only after passing the technical and purity checks.
- **Paired design:** if both fractions come from the same animals, joint models are possible (`~ animal + fraction:treatment`). Keep per-fraction DE as primary.
- **Covariates:** parsimonious only (e.g. one purity axis, one technical axis). Check rank and VIF. Report baseline vs adjusted effect retention per gene.
- **Robustness:** at n ≈ 6 per group, add an exact within-sex label permutation (20 × 20 = 400 relabelings for 3F + 3M per group) for key genes and modules. Use limma `camera` for gene sets.
- **Classify findings in tiers:**

  | Tier | Meaning |
  |---|---|
  | 1 | Consistent across baseline, adjusted and permutation; both sexes |
  | 2 | Robust to adjustment, but weaker in some checks |
  | 3 | Explained by purity or composition |
  | 4 | Confounded or uninterpretable |

### 5.7 Biology to examine
**SPP1 efficacy (top priority):**
- **Spp1 mRNA in CD11b+ cells, SPP1 vs NTC.** Spp1 is a DAM / MGnD / plaque-associated microglia gene, high in 5xFAD.
- **Spp1 receptors and downstream:** Cd44, Itgav, Itgb1, Itgb3, Itgb5 in both fractions.
- **Is the knockdown fraction-specific?** Compare knockdown in microglia vs astrocytes.

**Microglia (CD11b+):**

| Program | Genes / source |
|---|---|
| Homeostatic | P2ry12, Tmem119, Cx3cr1, Sall1 |
| DAM stage 1 | Tyrobp, Apoe, B2m, Ctsb, Ctsd, Lyz2 |
| DAM stage 2 | Cst7, Itgax, Clec7a, Lpl, Spp1, Gpnmb, Axl, Cd9, Igf1, Csf1 (Keren-Shaul 2017) |
| MGnD | Krasemann 2017 |
| ARM / IRM | Sala Frigerio 2019 |
| Interferon-responsive | Ifit1, Ifit3, Isg15, Irf7, Rtp4 |
| MHC-II | H2-Aa, H2-Ab1, Cd74 |
| Proliferation | Mki67, Top2a |
| Lipid / cholesterol handling | Apoe, Lpl, Abca1, Ch25h |
| Phagocytosis / lysosome | — |
| Complement | C1qa/b/c, C3ar1 |
| Inflammatory | Tnf, Il1b, Nfkbia |

- **BAM vs microglia shifts** within CD11b+.

**Astrocytes (ACSA-2+):**

| Program | Genes / source |
|---|---|
| Disease-associated astrocytes (DAA) | Gfap, Vim, Serpina3n, Ctss, C4b, Cd9, Apoe (Habib 2020) |
| Neurotoxic / "A1" | C3, Serping1, H2-D1, Gbp2, Fkbp5, Ugt1a1, Psmb8, Srgn, Amigo2 (Liddelow 2017) |
| "A2" / pan-reactive | Lcn2, Steap4, S1pr3, Timp1, Hspb1, Cxcl10, Cd44, Osmr, Cp, Serpina3n, Aspg, Vim, Gfap (Zamanian 2012; A1/A2 is now considered an oversimplification, so report as modules) |
| Interferon | — |
| Synapse support | Sparcl1, Gpc4, Gpc6, Thbs1, Thbs2 |
| Glutamate / K+ homeostasis | Slc1a2, Slc1a3, Kcnj10, Glul |
| Cholesterol synthesis | Hmgcr, Hmgcs1, Sqle, Fdft1, Srebf2 |
| Barrier / endfoot | Aqp4, Gja1 |

**Cross-talk (if paired fractions):** compare ligand expression in one fraction against receptor expression in the other:
- C3 (astro) – C3ar1 (microglia)
- Il33 – Il1rl1
- Spp1 – Cd44 / integrins
- Csf1 / Il34 – Csf1r
- Apoe – Trem2 / Lrp1
- Tgfb1/2 – Tgfbr1/2
- Cx3cl1 – Cx3cr1
- Gas6 / Pros1 – Mertk / Axl

Treat these as hypotheses. Bulk RNA cannot show signaling.

**Compare with the dura results:**
- Is interferon up with SPP1 here too?
- Is Tgfb2 lower in 5xFAD here too?

### 5.8 Pathology / severity (if available)
- **Within 5xFAD:** relate DAM, DAA and Spp1 levels to plaque load. Within-genotype associations avoid the genotype-contrast confound.

## 6. Deliverables
- **Folders:** `01_QC/`, `02_purity_celltyping/`, `03_artifacts/`, `04_exploration/`, `05_DE_CD11b/`, `06_DE_ACSA2/`, `07_modules/`, `08_crosstalk/`, `09_pathology/`, `tables/`, `figures/`, `scripts/`, `logs/`.
- **Files:**
  - `ANALYSIS_README.md`: full methods and results;
  - `EXECUTIVE_SUMMARY.md`: 1–2 pages, split into high confidence / suggestive / likely artifact / cannot determine;
  - `analysis_summary.tsv`: one row per finding with tier and supporting file;
  - a figure-to-data map.
- **Figures:** show individual animals (circle F, triangle M). Use consistent group colors. The previous palette was WT_NTC `#2a78d6`, 5xFAD_NTC `#eb6834`, 5xFAD_SPP1 `#1baf7a`.

## 7. Reporting rules
- **State failure conditions explicitly:** group × batch confounding, technical separation, purity confounding, too few replicates, unverifiable genotype, failed marker validation.
- **Never remove samples** to strengthen a result. Document any exclusion with an objective, pre-specified reason.
- **Correlation, adjustment and purity:** a difference that disappears after adjustment is not proven false. A difference that persists is not proven cell-intrinsic.
- **Effect sizes and confidence matter more than p-values.** With n = 6 per group, report effect sizes and say what can and cannot be ruled out.
- **Ask the user** whenever a design detail changes the analysis (Section 3).
