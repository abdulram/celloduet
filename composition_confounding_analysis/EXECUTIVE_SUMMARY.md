# Executive summary: 5xFAD dura bulk RNA-seq, after accounting for dissection and composition

**Question.** After accounting for probable differences in dura dissection and cellular composition, what can we conclude about 5xFAD-associated dural biology?

**Answer.** Essentially nothing, with confidence. The experiment cannot separate genotype from two other things:
1. **How the tissue was dissected.** 5xFAD samples contain less dura stroma and blood endothelium and more bone and adjacent tissue.
2. **How the libraries were processed.** All 12 5xFAD libraries share a distinct mapping/mitochondrial profile that only one WT library (616) shares.

The apparent 5xFAD signature (6,787 genes) is best explained by these factors. This is a useful negative result: it shows what the next experiment must control.

Design: WT n = 6 vs 5xFAD n = 6 (plus 6 SPP1-treated 5xFAD for context), 3F/3M each. No age, batch, dissection or pathology metadata.

---

## HIGH CONFIDENCE
- **Most transcriptome variation is composition.** Transcriptome PC1 (40%) is a dura-purity axis: fibroblast/endothelial vs olfactory/nasal, CNS and muscle carry-over. PC2 (23%) is bone plus library profile. Measured composition explains 90% and 83% of them.
- **Genotype is associated with the major composition axes, with poor overlap:**
  - dura purity: AUC 0.08, p = 0.015;
  - bone: AUC 0.94, p = 0.009.
  WT and 5xFAD animals barely overlap, so regression adjustment cannot be trusted to fix this.
- **Genotype is almost perfectly confounded with library processing:** % mito 5% vs 10%, multi-mapping 10% vs 17%. The technical sentinel modules (ribosome, mitochondrial, OXPHOS) are the strongest "genotype" modules.
- **The decisive check:** WT 616, the only WT library processed like the 5xFAD libraries, looks 5xFAD-like on 84% of the 6,787 "genotype" genes. The other WT animals do so on 2–5%. Either the signature is technical, or 616 is mis-genotyped. In both cases genotype is not separable.
- **Adjustment collapses the signature.** Adjusting for % mito leaves 1 significant gene; adding the bone axis leaves 0. Pathway tests that account for correlation between genes (camera) find nothing at any FDR.

## SUGGESTIVE (hypotheses only; not significant after correction)
- **Lipid-handling module:** higher in 5xFAD once composition is adjusted (+1.6 SD, p ≈ 0.02–0.04 in module models; camera FDR 0.99). Driven mainly by Fabp5; Apoe goes the opposite way.
- **DAM-like trend:** +1.1 SD, p = 0.07. Its genes are incoherent: Itgax, Axl and Apoe go down; Spp1, Fabp5 and Tyrobp go up. Spp1 is mostly bone-derived here.
- **Lower macrophage MHC-II at equal macrophage abundance** (permutation p = 0.02). This is lost after technical adjustment.

## LIKELY COMPOSITION / DISSECTION-DRIVEN
- **Lower blood-endothelial, dendritic and B-cell signal, and higher bone signal, in 5xFAD samples.** These could be disease remodeling or dissection, and are indistinguishable here.
- **"Fibroblast activation" up and "angiogenesis" down.** These are explained by the bone axis and by endothelial abundance / technical axis.
- **The genome-wide 5xFAD signature as a whole.**

## CANNOT DETERMINE FROM THIS DATASET
- Any cell-intrinsic 5xFAD effect in dural fibroblasts, endothelium, macrophages or lymphocytes.
- **Lymphatic endothelial, T-cell, plasma-cell and mural biology.** Their marker sets do not behave coherently in these samples, so neither abundance nor state can be estimated.
- Relationships with disease severity or age (no pathology or age data), and whether the 5xFAD animals are transgenic (not visible in the count matrix).

## What would make this answerable
1. **Genotype animal 616** and check processing records (RNA-prep and library batches, RIN).
2. **Process WT and 5xFAD interleaved** in the same batches.
3. **Standardize dura dissection** (peel from calvaria; exclude cribriform/olfactory region and bone).
4. **Measure composition directly** with flow cytometry or whole-mount immunofluorescence.
5. **Use scRNA-seq/snRNA-seq** for state vs abundance.
6. **Record per-animal pathology and age.**
7. **Replicate in an independent cohort.**

Full methods and numbers: `ANALYSIS_README.md`. Machine-readable findings: `analysis_summary.tsv`. Central figure: `10_final_figures/Fig12_pathway_robustness_summary.png`.
