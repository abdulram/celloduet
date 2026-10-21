# Brief for making a visual summary: SPP1 knockdown in 5xFAD glia

Use this file to build **one simple visual summary**, such as a single figure, a poster panel or a slide.

- Everything you need is here; you do not need the raw data.
- Use only the numbers below and do not add new claims.
- Keep the wording conservative: these are mRNA results from sorted cells, with no histology or protein data.

---

## 1. One-sentence takeaway (use as the title or headline)
**SPP1 treatment lowers microglial Spp1 by about 43% without changing the microglial, astrocyte or oligodendrocyte response to amyloid pathology.**

## 2. Study design (show as small chips or a caption)
| Item | Value |
|---|---|
| Model | 5xFAD mice (all animals) |
| Sex | all male (checked from Xist and Y-chromosome genes) |
| Groups | NTC (non-targeting control) n = 3; SPP1 treatment n = 3 |
| Tissue | whole brain, warm (enzymatic) dissociation |
| Sorts | CD11b⁺ (microglia); ACSA-2⁺ (astrocyte-enriched), sorted after CD11b⁺ depletion |
| Assay | bulk RNA-seq, about 20 M reads per sample; libraries re-run after poor RNA quality the first time |

Sort quality, which is useful context for the figure:
- **CD11b⁺:** 95–99% microglial RNA (clean).
- **ACSA-2⁺:** only about 35% astrocyte RNA, about 55% oligodendrocyte and about 6% ependymal. ACSA-2 results therefore describe mixed glia.

---

## 3. Suggested layout: 4 panels, read left to right or top to bottom

### Panel A: "The knockdown works" (dot plot, one dot per animal)
Spp1 in CD11b⁺ microglia, counts per million (CPM):

| Group | Animal 1 | Animal 2 | Animal 3 | Mean |
|---|---|---|---|---|
| NTC | 1608 | 2385 | 2156 | 2049 |
| SPP1 | 1187 | 1024 | 1311 | 1174 |

- **Annotation:** "−43%; no overlap between groups; adjusted p = 2×10⁻⁷ (DESeq2)."
- **Footnote:** "Robust to dropping any animal and to adjusting for sort purity, library quality and dissociation."

### Panel B: "…and it is specific" (horizontal bars, % change SPP1 vs NTC)
Other disease-associated microglia (DAM) genes did not change:

| Gene | % change |
|---|---|
| **Spp1** | **−43** (highlight) |
| Cst7 | 0 |
| Itgax | −3 |
| Clec7a | −10 |
| Lpl | −8 |
| Trem2 | −2 |
| P2ry12 (homeostatic) | 0 |

**Optional note:** "Spp1 relative to each animal's own DAM level is roughly halved."

### Panel C: "No change in 5xFAD pathology programs" (forest plot: % change with 95% CI)
Scores use published gene sets, with Spp1 removed from every set.

| Program (cell type) | % change | 95% CI |
|---|---|---|
| DAM stage 2 (microglia) | −3 | −15 to +10 |
| Plaque-induced genes (microglia) | −4 | −11 to +3 |
| MGnD (microglia) | −10 | −23 to +6 |
| Homeostatic (microglia) | −1 | −11 to +11 |
| Phagocytosis genes (microglia) | +3 | −13 to +22 |
| Reactive astrocytes, DAA | −3 | −24 to +24 |
| Disease-associated oligodendrocytes | +1 | −22 to +30 |

- Draw a vertical line at 0. All intervals cross it.
- **Caption:** "The DAM response scales with plaque exposure, so a large change in plaque load is unlikely (histology needed to confirm). With n = 3 per group, the smallest detectable change is about ±10–17% for microglial programs and about ±30% for astrocyte/oligodendrocyte programs."

**Also unchanged.** A text box or a row of grey "no change" chips is enough for these:
- MMPs (Mmp14 0%; Mmp9/Mmp8 essentially not expressed);
- phagocytic machinery (Trem2, Tyrobp, Mertk, C1q, Cd68);
- ECM production by microglia (−1% as a group);
- Spp1 receptors (Cd44, integrins Itgav/b1/b3/b5).

### Panel D: "One pattern worth testing" (horizontal bars; visually mark as a hypothesis, e.g. amber or a dashed border)
Lipid-processing genes are lower, while lipid uptake and efflux are not:

| Gene / set | Role | % change |
|---|---|---|
| Gpnmb | lipid/lysosomal stress | −40 |
| Ch25h | oxysterol (25-hydroxycholesterol) | −29 |
| Apoe | lipid transport | −22 |
| Fabp5 | fatty-acid binding | −19 |
| Cholesterol-synthesis genes (19-gene set) | sterol synthesis | −12 (18 of 19 genes down) |
| Trem2 | lipid sensing/uptake | −2 (grey) |
| Lrp1 | uptake | −3 (grey) |
| Abca1 | cholesterol efflux | −3 (grey) |

**Caption:**
- Two readings are possible: fewer lipid-laden microglia, or Spp1 acting upstream of part of the lipid response.
- Test with BODIPY/Plin2 or filipin staining near plaques, or with single-cell data.
- Individually these changes are modest or nominal.

---

## 4. Optional small strip: "Set aside (explained by sample handling, not treatment)"
Use a grey footer strip with one line each:
- Most microglia "DE genes" were contaminating astrocyte, oligodendrocyte or vascular RNA (SPP1 sorts are 95.5% vs 98.2% microglial).
- ACSA-2⁺ differences are mostly cell mix. An apparent interferon increase disappears when expressed per astrocyte.
- Collagen changes follow fibroblast/vascular carry-over (up in one sort, down in the other).
- Stress-hormone genes (Fkbp5, Tsc22d3, Sgk1) are driven by one SPP1 animal.

---

## 5. Style guidance
- **Group colours** (from the earlier project, keep consistent): NTC `#eb6834` (orange), SPP1 `#1baf7a` (green).
- **Unchanged genes and programs:** grey. **Hypothesis panel:** amber (`#b7791f`) or a dashed outline.
- **Show individual animals** wherever a per-animal value exists. The earlier project's convention was circles for females and triangles for males; all animals here are male.
- **Keep it simple:** one headline, four panels, one caption line each, and no p-value clutter beyond the Spp1 panel.
- **Wording:** say "lower/higher mRNA", not "reduced function". Say "no detectable change", not "no change".

## 6. Do not claim
- Effects on plaque load, cognition or protein levels. None were measured.
- Astrocyte knockdown of Spp1. Astrocyte Spp1 is about 2% of the microglial level and mostly carry-over, so it cannot be assessed.
- Anything about WT vs 5xFAD or sex differences. All animals are 5xFAD males.

## Source
Analysis in repo `abdulram/celloduet`, branch `claude/new-session-2e6idk`, folder `sorted_glia_55ZS67/`. Full tables are in `tables/analysis_summary.tsv`, `09_pathology/` and `07_modules/`.
