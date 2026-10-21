# Executive summary: 55ZS67 sorted CD11b+ / ACSA-2+ bulk RNA-seq (SPP1 vs NTC)

**Status:** first pass, built from the count matrix only. Several design facts are **inferred, not confirmed** (see "Open questions"). Treat every tier as provisional until they are answered.

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
- **ACSA-2+ interferon module up with SPP1** (camera FDR 8e-4; Ifit3 +0.52, Ifit3b +0.58, each n.s. alone). This matches the dura finding (type-I IFN up with SPP1), but here it is small and set-level only.

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
1. **Genotype.** Are all 12 animals 5xFAD? What is their age? Is any WT tissue sequenced elsewhere?
2. **SPP1 modality.** Is it an ASO, AAV-shRNA, CRISPRi or siRNA? What is the route and the time from dosing to harvest? Is there a vector sequence (for reads in unmapped FASTQs)?
3. **Pairing.** Do replicate numbers mean the same animal across fractions?
4. **Sorting.** MACS or FACS? Was CD11b depleted before ACSA-2? Enzymatic at 37 °C, or cold? Were inhibitors used? Region?
5. **Why "rerun".** Re-sequenced or re-prepped libraries? Are there original runs? What are the batches, RIN and cell numbers?
6. **Raw data.** Can you provide FASTQs for transgene confirmation and gene-body coverage?
7. **Pathology.** Are plaque load or Iba1 data available per animal?
