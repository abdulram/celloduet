# Does SPP1 change 5xFAD-associated changes in dura fibroblasts?

**Approach:** use only fibroblast-specific genes (151 genes from `results/fibroblast/`, fibroblast-normalized). First find real 5xFAD vs WT differences in those genes, then ask whether SPP1 shifts 5xFAD animals along them.

Code: `analysis/fibroblast_5xfad_spp1.py`.

## TL;DR
1. **Genotype can't be verified from this data.** The expression matrix contains mouse genes plus ERCC spike-ins only. Reads from the human APP/PSEN1 transgenes have no feature to be counted against, and no indirect trace exists (details below). The FASTQ or BAM files are needed.
2. **No fibroblast gene shows a 5xFAD vs WT difference that can be separated from the technical axis.**
   - 60 of 151 genes differ at padj < 0.05.
   - For those genes, WT animal 616 falls with the 5xFAD mice (median position 0.99 on a WT = 0 → 5xFAD = 1 scale). 616 has the same library/technical profile as the 5xFAD animals. Ordinary WT animals, left out one at a time, sit at about 0.
   - **None** survive adjustment for % mito reads.
   - The 6 genes where 616 looks WT-like (Slc23a2, Fmod, Cilp2, Fmo1, P3h3, Matn4) are about what chance predicts: 10%, vs 2–15% for ordinary WT animals.
3. **SPP1 does not shift 5xFAD fibroblasts toward or away from WT** on any version of the 5xFAD signature:

   | Signature | SPP1 − NTC (score units) | Share of the 5xFAD–WT gap | Permutation p |
   |---|---|---|---|
   | All 60 genotype genes | −0.11 | 3% | 0.66 |
   | 6 technically robust candidates | −0.14 | 6% | 0.72 |
   | All 151 genes, gene-wide concordance with the 5xFAD effect (% mito adjusted) | ρ = −0.00 | — | 1.00 |

   The p-values come from all 400 relabelings of SPP1/NTC within sex. Only one gene hints at a rescue: Matn4 is lower in 5xFAD (−0.41 log2FC unadjusted, not robust to the technical adjustment) and higher with SPP1 (+0.55, padj 0.048).

**Bottom line:** in this data set there is no detectable 5xFAD signature in dura fibroblasts, so there is nothing for SPP1 to modify. There is also no evidence that SPP1 moves the 5xFAD dura-fibroblast state. The SPP1 gene-level changes reported earlier (Cldn11, Thbd, Fbln1 down; Matn4, Ifit3b up) look like SPP1 effects that are independent of the genotype axis.

## A pitfall that was caught
The first version defined the 5xFAD signature as **WT vs 5xFAD-NTC** and then tested **SPP1 vs NTC**. That design is biased toward finding a "rescue". Genes are selected where the NTC animals happen to be extreme, so any independent group (here SPP1) regresses back toward WT. It produced an apparent 55% rescue with p = 0.01.

The signature is now defined as **WT vs all 12 5xFAD animals pooled**, which treats SPP1 and NTC symmetrically. The apparent rescue then vanishes (3–6%, p ≈ 0.7).

## Genotype verification
* **No human features in the matrix.** It contains only mouse (ENSMUSG) genes and ERCC spike-ins, so human APP/PSEN1 reads are either unmapped or partially mismapped.
* **Mouse App and Psen1 are *lower* in 5xFAD** (App 339 vs 530 CPM). But WT 616, the technically 5xFAD-like WT library, is equally low (334). This is the technical axis, not a transgene effect.
* **Thy1 doesn't differ.** The 5xFAD cassette carries mouse Thy1 exons, so Thy1 is the other possible trace, but 53 vs 41 CPM is noise.
* **Neurons are too sparse here** (Rbfox3 2–14 CPM) for the neuron-specific transgenes to be visible anyway.
* **To verify:**
  1. Genotyping records.
  2. FASTQ/BAM files: count reads matching human APP/PSEN1 sequence, especially the Swedish (K670N/M671L), Florida (I716V) and London (V717I) sites. With FASTQs this is quick.
  3. Brain tissue from the same animals, where the Thy1-driven transgenes are highly expressed.

## Caveats
* This design can only detect genotype effects larger than what one technically 5xFAD-like WT animal (616) can distinguish. Small real effects may be hidden by the batch confound.
* The absence of a 5xFAD signature could also be biological: dura fibroblasts may not respond at this age or disease stage. Animal age is not in the metadata.
* The fibroblast gene set depends on the deconvolution, which is approximate (see `results/fibroblast/FIBROBLAST_REPORT.md`).

## Files
| File | Contents |
|---|---|
| `fibroblast_genes_genotype_and_spp1.csv` | Per gene: 5xFAD vs WT (unadjusted and mito-adjusted), 616 position, SPP1 vs NTC, call |
| `t_leave_one_out_WT.csv` | Null for the 616 test (each ordinary WT animal treated as the test animal) |
| `signature_scores_per_sample.csv`, `spp1_signature_permutation_tests.csv`, `spp1_vs_5xFAD_concordance.csv` | SPP1 tests |
| `G1`–`G3` figures | 616 test and effect shrinkage, top genes per sample, SPP1 on the signature |
