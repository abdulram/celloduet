# Executive summary: 55ZS67 sorted CD11b+ / ACSA-2+ bulk RNA-seq (SPP1 vs NTC)

**Status:** round 2 (pathology deep-dive added). Round 1 was a first pass, built from the count matrix only. Several design facts are **inferred, not confirmed** (see "Open questions"). Treat every tier as provisional until they are answered.

## Design facts confirmed by the user (round 2)
- All 12 animals are **5xFAD**. Brain regions: **all (whole brain)**.
- **Warm (enzymatic) dissociation.** Ex vivo activation is expected, and it is present.
- **CD11b+ cells were depleted before ACSA-2 enrichment.** The ~1% microglial RNA in ACSA-2 is residual carry-over that escaped depletion.
- **"rerun":** the first libraries had poor RNA quality and were re-sequenced. A gene-length check finds no degradation bias in CD11b. ACSA-2 shows a mild bias, the same in both arms (see `03_artifacts/gene_length_quality_index.tsv`).
- **Still unconfirmed:** pairing of fractions by replicate number; SPP1 modality and timing; age.

## Does 5xFAD pathology differ between SPP1 and NTC? (round 2, `09_pathology/`)
**No pathology-associated program differs detectably.**
- Scored 20 published gene sets, with Spp1 removed from each set and ex vivo dissociation genes excluded:
  - DAM 1/2, MGnD, ARM;
  - plaque-induced genes (PIG);
  - homeostatic, lysosome, phagocytosis, complement, MHC-II, IFN, proliferation and lipid-handling programs;
  - astrocyte DAA / pan-reactive / A1 programs;
  - disease-associated oligodendrocytes (DOL).
- All 95% CIs include 0. With 3 vs 3 the study can detect about **±17% (DAM2)** and **±10% (PIG)** in microglia, but only about **±30%** for astrocyte and oligodendrocyte programs.
- **Microglial plaque-response magnitude is unchanged:**
  - DAM2 −3% (CI −15 to +10%);
  - PIG −4% (CI −11 to +3%);
  - homeostatic −1%.
  - DAM intensity scales with plaque exposure in 5xFAD, so this argues **against a large change in amyloid burden or microglia–plaque engagement**. That is an indirect inference; histology is needed to confirm it.
- **Spp1 falls per unit of DAM response.** Spp1 relative to the DAM2 score is 0.96–1.68 (log2) in NTC vs 0.46–0.70 in SPP1, with no overlap. The knockdown is specific, not a weaker plaque response.
- **Small, coherent but non-significant microglial dips sit in the Spp1-linked lipid/lysosome axis:** Apoe −22%, Gpnmb −40%, Ctsb −11%, lysosome set −7%, MGnD −10%, DAM1 −6% (all six animals separate on DAM1, but the effect is tiny).
- **Astrocyte and oligodendrocyte disease states are unchanged.** After normalising to astrocyte content (astrocyte RNA is ~20% higher in SPP1 ACSA-2 sorts):
  - DAA −18% (p = 0.26), pan-reactive −16%, A1 −10%;
  - DOL −3%.
  - The ACSA-2 **interferon increase from round 1 disappears (−0.1%)**. It was composition, and is now Tier 3.

## MMPs, phagocytosis, lipid handling in microglia (round 3)
- **MMPs: no microglial change.**
  - Only Mmp14, Mmp12 and Mmp2 are expressed. Mmp12 +26% is nominal and tracks BAM content.
  - Mmp9/Mmp8 come from neutrophil contamination in one sample.
- **Phagocytosis: recognition and engulfment genes are unchanged** (Trem2, Tyrobp, Mertk, C1q, Cd68, NADPH oxidase).
  - Several lysosomal/phagocytic genes are ~10–15% lower, each only nominal: Axl, Fcgr2b, Lamp2, Ctss/Ctsh/Ctsc/Ctsb.
  - ~3% of that is dilution from the purity difference.
- **Lipid handling: a coherent mild decrease in the lipid-associated arm.**
  - Apoe −22%, Gpnmb −40%, Ch25h −29%, Fabp5 −19%, Fabp3 −18%, Soat1 −12%.
  - Cholesterol-synthesis (SREBP2) targets −12%, with 18/19 genes down even though astrocyte carry-over would push them up.
  - Lipid uptake (Trem2, Lrp1, Cd36, Lpl) and efflux (Abca1) are unchanged.
  - Tier 2 hypothesis: fewer lipid-laden microglia, or Spp1 upstream of part of the lipid program. Needs lipid-droplet staining or single-cell data.

## What the data are (inferred from file and sample names)
- 12 libraries, provider order 55ZS67 (all names carry the prefix `rerun_`). No FASTQ/BAM files.
- 2 sorted fractions × 2 treatments × n = 3:
  - CD11b+ ("microglia"): samples 1–3 NTC, 4–6 SPP1.
  - ACSA-2+ ("astrocyte"): samples 7–9 NTC, 10–12 SPP1.
- **No WT group and no 5xFAD-vs-WT contrast is possible.** Genotype is not stated anywhere in the files.
- **All 12 samples are male** (Xist ≈ 0 CPM; Ddx3y/Uty/Kdm5d/Eif2s3y expressed).
- The facility annotation contains only mouse genes plus 92 ERCC spike-ins. The ERCC spike-ins carry essentially zero reads, so they were probably not added. **Human APP/PSEN1 transgenes cannot be checked from this matrix.**
- **Likely paired design.** Residual expression profiles link CD11b replicate *k* with ACSA-2 replicate *k* in both arms (matched pairing ranks 1st of 6 in both arms; chance probability ≈ 1/36). Replicate *k* is probably the same animal or processing slot. Please confirm.

## High confidence
1. **The CD11b+ sort is pure microglia:** about 96–99% microglial RNA by deconvolution. Astrocyte, oligodendrocyte and neuron markers are about 1% of their ACSA-2 levels.
2. **The CD11b+ microglia carry a strong DAM program** (Cst7 ≈ 1,900 CPM, Itgax ≈ 250, Clec7a ≈ 480, Spp1 ≈ 2,000). This fits an amyloid model such as 5xFAD, but without a WT arm it is not proof of genotype.
3. **Spp1 mRNA is lower in SPP1-treated microglia:** 1,174 vs 2,049 CPM, about −43%.
   - DESeq2 log2FC −0.81, padj 2e-7. Spp1 is the #1 gene genome-wide.
   - The 3 SPP1 and 3 NTC samples do not overlap (exact permutation p = 0.1, which is the smallest p possible with 3 vs 3).
   - Survives every leave-one-out fit and every covariate model (log2FC −0.55 to −1.21).
   - Carry-over predicts 0 change for Spp1.
   - Other DAM genes are unchanged (Cst7 −0.02, Itgax −0.05, Lpl −0.13, Trem2 −0.03), so this is not a loss of the DAM state. **This is the first direct evidence that the SPP1 treatment lowers Spp1 in the cells that make it.** Effect size is moderate (partial knockdown).

## Suggestive (Tier 2: consistent and robust to adjustment, but n = 3 and correlated with the Spp1 drop)
- **A small set of lipid/plaque-associated microglial genes falls with Spp1:**
  - Gpnmb −40%, Ch25h −29%, Fabp5 −18%, Apoe −22% (n.s.), Igf1 −16%;
  - also Rgs16, Lgi2, Mcoln3, Hpse, Mrc2, Lat2.
  - All are microglia-enriched (carry-over cannot explain them) and robust to leave-one-out.
  - Camera: lipid-handling module down (FDR 1e-5); DAM2 without Spp1 down (FDR 1e-3, driven by Gpnmb/Igf1/Axl). The core DAM genes are flat.
  - **Interpretation:** either a downstream effect of Spp1 loss or a shift in a Spp1^hi/Gpnmb^hi plaque-associated subset. Bulk data cannot separate the two.
- ~~ACSA-2+ interferon module up with SPP1~~: **downgraded to Tier 3 in round 2.** After normalising to astrocyte content it is −0.1%, so it is explained by the higher astrocyte share in SPP1 sorts.

## Likely artifact (Tier 3–4)
- **About 60% of the 116 CD11b "DE genes" (padj < 0.1) are carry-over** (astrocyte/oligodendrocyte genes such as Gpr37l1, Atp1a2, Slc1a2, Ptn, Sparcl1, Plp1, Mobp).
  - SPP1 CD11b sorts carry slightly more non-microglial RNA: astrocyte marker share 1.5% vs 0.95%, with complete separation.
  - They also carry more border-associated macrophage (BAM) RNA (1.8% vs 0.35%; Pf4, Cd163, Mrc1 up).
- **About 20% of DE genes in each fraction are processed pseudogenes or unnamed loci.** They track a library-quality axis (% unmapped, % assigned), and that axis is **perfectly separated by treatment within each fraction**:
  - CD11b: SPP1 has more unmapped reads (4.1% vs 2.8%) and fewer assigned (51% vs 60%).
  - ACSA-2: NTC has more unmapped reads (7.0% vs 4.1%).
  - This is the same kind of confound as in the dura project, but milder. These genes are Tier 4.
- **ACSA-2+ DE is dominated by composition:**
  - Ependymal/cilia genes (Dnah6/12, Spag17, Hydin, Sspo, Cfap43) are lower in SPP1 (ependymal 5.3% vs 7.3%).
  - Microglial carry-over is lower in SPP1 (1.0% vs 1.6%).
  - Astrocyte share is higher in SPP1 (38% vs 32%).
  - OXPHOS/fatty-acid Hallmark "up" sets line up with this composition and quality axis.
- **A glucocorticoid/stress program is up in SPP1 in both fractions** (Fkbp5, Tsc22d3, Sgk1, Klf9, Per1, Mt2, Cdkn1a):
  - Clock genes (Dbp, Nr1d1, Per3) are flat, so this is not time-of-day.
  - It is driven mostly by SPP1 replicate 2 (samples 5 and 11), which supports pairing.
  - **Possible causes:** stress or handling at harvest, or a systemic treatment effect. It cannot be separated from treatment here.
- **Dissociation immediate-early genes (IEGs) are high in both fractions** (Fos ~3,000 CPM, Ccl3/Ccl4 ~2,000–3,500 CPM in microglia). The level is similar between arms, so it is not a group confound, but ex vivo activation is present.

## Cannot determine
- **Astrocyte knockdown of Spp1.** ACSA-2 Spp1 is 27–48 CPM, about 2% of the microglial level.
  - Most of it is consistent with microglial carry-over.
  - Spp1/Hexb is unchanged (0.22 vs 0.22), so there is no detectable knockdown in that fraction.
- **Any astrocyte-intrinsic SPP1 effect.** The ACSA-2+ fraction is only about 30–40% astrocyte RNA; about 55% is oligodendrocyte RNA (Plp1 ≈ 20,000 CPM) and about 6% ependymal. ACSA-2 DE is therefore a mixed-glia contrast.
- **5xFAD vs WT effects, and the Tgfb2 question from dura.** There is no WT arm. Tgfb2 is unchanged by SPP1 in either fraction.
- **Sex effects.** All animals are male.
- **Spp1 receptors** (Cd44, Itgav/b1/b3/b5): no change in either fraction.
- **Dura SPP1 genes:**
  - Rgs4, Cldn11 and Prg4 are not changed in CD11b. The apparent Cldn11 "up" is oligodendrocyte carry-over.
  - Matn4 is not expressed in CD11b.

## Power statement
With n = 3 vs 3:
- the smallest exact permutation p is 0.1;
- DESeq2 p-values rest on the dispersion prior;
- an effect smaller than about 30% (|log2FC| < 0.5) for a typical gene cannot be ruled out or confirmed.

## Open questions for you (answers change the analysis)
1. ~~Genotype~~: answered, all 5xFAD. Still open: age; any WT sorted data?
2. **SPP1 modality.** Is it an ASO, AAV-shRNA, CRISPRi or siRNA? What is the route and the time from dosing to harvest? Is there a vector sequence (for reads in unmapped FASTQs)?
3. **Pairing.** Do replicate numbers mean the same animal across fractions?
4. **Sorting** (partly answered: CD11b depleted first, warm dissociation, whole brain). Still open: MACS or FACS? Transcription inhibitors?
5. ~~rerun~~: answered (poor RNA the first time). Still open: RIN/DV200 of the rerun, and whether arms were extracted or prepped on different days.
6. **Raw data.** Can you provide FASTQs for transgene confirmation and gene-body coverage?
7. **Pathology.** Are plaque load or Iba1 data available per animal?
