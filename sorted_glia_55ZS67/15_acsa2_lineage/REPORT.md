# ACSA-2+ fraction: astrocyte-only and oligodendrocyte-only genes, and the link to each mouse's microglial Spp1 knockdown

Script: `scripts/15_acsa2_lineage_specific.py` (log: `logs/15.log`).

**Design.** All mice are 5xFAD males, 3 per arm, given L2-siRNA (Spp1-targeting or non-targeting), 15 nmol by a single ICV injection at 6 months, and harvested at 8 months. CD11b+ and ACSA-2+ (MACS) fractions come from the same mouse.

## How "lineage-only" was defined
The ACSA-2 fraction is about 55% oligodendrocyte, 35% astrocyte and 6% ependymal RNA (reference deconvolution). For each gene, the share of its ACSA-2 signal expected from each cell type is: measured RNA share × that cell type's expression in the brain reference.

- **Astrocyte-only:** ≥ 90% of the signal from astrocytes; **131 genes**, e.g. Agt, Aldoc, Slc1a2, Fgfr3, Btbd17, Lcat, Dio2, Fam107a.
- **Oligodendrocyte-only:** ≥ 90% from oligodendrocytes; **296 genes**, e.g. Mobp, Mog, Mag, Plp1, Opalin, Ermn, Enpp6, Trf.
- Genes also needed mean CPM ≥ 5 and had to be protein-coding. Pseudogenes and Gm/Rik loci were dropped because they track a library-quality difference between arms.
- An 80% cutoff (341 / 614 genes) was run as a sensitivity check and gave the same conclusions.

**DE method.** Each lineage was tested on its own genes, normalised to those genes (PyDESeq2 control-gene size factors). A change in how much astrocyte RNA a sort captured therefore cannot look like regulation. The flip side: a uniform up- or down-shift of a whole lineage is absorbed by this normalisation, so it is reported separately as "content".

## Results

### 1. How much of each lineage the sorts captured (15a)

| Lineage | SPP1 / NTC | Welch p |
|---|---|---|
| Astrocyte | 1.23× | 0.10 |
| Oligodendrocyte | 1.02× | 0.89 |
| Ependymal | 0.68× | 0.29 (one NTC sort low) |

SPP1 sorts captured somewhat more astrocyte RNA. That is why earlier library-normalised ACSA-2 results looked like astrocyte genes going up.

### 2. Astrocyte-only genes, per astrocyte (15b)
**No change.** 0 of 131 genes reach q < 0.1, the median log2FC is −0.01, and 47% go up.
- Nominal p < 0.05, unreliable at 3 vs 3:
  - Down: Trim9 (−26%), Agt (angiotensinogen, −26%), Rpe65.
  - Up: Ppp1r3g, Fam107a (+27%).
- Fam107a (also called DRR1) is glucocorticoid-induced. It fits the stress signature seen in both fractions, which is driven mainly by SPP1 mouse 2.
- Astrocyte programs per astrocyte (only genes ≥ 50% astrocyte-derived): interferon +15% (p 0.09), synapse −18% (p 0.12), pan-reactive +4% (p 0.81). All q ≥ 0.24.
- DAA, A1 and endfoot programs could not be scored: too few of their genes are astrocyte-specific in this mixed fraction.

### 3. Oligodendrocyte-only genes, per oligodendrocyte (15b)
**Essentially no change.** 1 of 296 genes reaches q < 0.05: Polr3e, up 2-fold (q 0.039). Polr3e is an RNA polymerase III subunit with no obvious link to SPP1, so this single gene needs replication.
- **Sgk1 +6-fold** (p 0.002, q 0.30). Sgk1 is glucocorticoid-induced: the stress signature again (mouse 2).
- Myelin structural genes +9% (p 0.098); Mobp +27% (p 0.04).
- Disease-associated oligodendrocyte (DOL) program +2% (p 0.77).
- Oligodendrocyte cholesterol synthesis +16% (p 0.30). This is where the earlier "sterol program up in ACSA-2" comes from, so it is opposite in direction to microglia.

### 4. Do astrocyte/oligodendrocyte changes track each mouse's own microglial knockdown? (15c)
The test asks whether the gene-wide SPP1 effect in a lineage (treatment log2FC) lines up with how genes vary *within* each arm according to that mouse's microglial Spp1. A negative ρ is the knockdown-consistent direction: the less Spp1 a mouse's microglia retained, the more SPP1-like its glia look. The null comes from all 36 within-arm shuffles of the mouse labels.

| Lineage | ρ | perm p | Leave-one-mouse-out ρ (any mouse dropped) | Drop stress genes |
|---|---|---|---|---|
| Astrocyte-only | −0.44 | 0.083 | −0.28 to −0.55 | −0.44 |
| Oligodendrocyte-only | −0.52 | **0.028** (minimum possible) | −0.19 to −0.55 | −0.52 |

**Controls:**
- **Knockdown index.** Spp1 relative to the microglial DAM2 program gives the same result (oligodendrocyte ρ −0.47, p 0.028). Within arm, microglial Spp1 tracks the knockdown index (r 0.96), not DAM2 or disease severity (r −0.18).
- **Negative controls:** using microglial DAM2 as the per-mouse predictor gives ρ ≈ 0. Using ACSA-2 library quality (% unmapped) also gives ρ ≈ 0, so the quality difference between arms doesn't explain it.

**Interpretation.**
- Oligodendrocyte-only, and more weakly astrocyte-only, expression differences scale with each mouse's degree of Spp1 knockdown. That is a dose-response pattern, from only 3 mice per arm.
- It is a gene-wide *coherence* signal, not a list of significant genes. It says the small, individually non-significant shifts are not random noise.
- **Important caveat.** Microglial knockdown is also a read-out of how much siRNA each mouse's brain received: ICV distribution varies between animals. So the link may mean "mice with more siRNA exposure show more glial change". That could act directly in oligodendrocytes and astrocytes, not via microglial SPP1 signalling. These data cannot separate the two.

## Bottom line
- **Astrocytes and oligodendrocytes show no detectable per-cell transcriptional change** with SPP1 knockdown at 2 months. At most, single genes or small programs move by about 10–25%, and none survive correction except Polr3e.
- **Their small shifts are not random:** they scale with each mouse's microglial knockdown, most clearly in oligodendrocytes. That is worth following up with more animals, or with a dose series.
- **Remaining caveats:** the glucocorticoid/stress genes (Sgk1, Fam107a) are a confound tied to one SPP1 mouse. With n = 3 vs 3, effects smaller than about 25–30% per gene cannot be ruled out.
