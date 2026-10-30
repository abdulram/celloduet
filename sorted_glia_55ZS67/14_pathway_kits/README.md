# Upload kits: PANTHER and QIAGEN IPA (sorted CD11b+ brain microglia, 5xFAD SPP1-targeted vs NTC)

Made by `scripts/14c_cd11b_links_crosstissue_kits.py`. Comparison: **SPP1 vs NTC in CD11b+ cells**, n = 3 vs 3, all male 5xFAD, whole brain. The step-by-step clicks are the same as for the dura kit (`../../composition_confounding_analysis/13_pathway_kits/README.md`); only the files and gene universes differ.

## Why there are two universes
SPP1 CD11b sorts carry slightly more non-microglial RNA than NTC sorts:

| Carry-over | SPP1 | NTC |
|---|---|---|
| Astrocyte markers | 1.5% | 0.95% |
| Border-associated macrophages (BAM) | 1.8% | 0.35% |
| Endothelial | ~1.3% | ~0.7% |

Genes made by those cells therefore shift with treatment for compositional reasons (Col3a1, Ccn1/2, Cldn11, Gpr37l1, Pf4, Mrc1…). Step 14a flags a gene as **microglia-intrinsic** when the measured carry-over cannot shift it by more than ~5%, and it is not BAM- or vascular/blood-enriched in the brain reference. This keeps 13,769 of 15,492 tested genes (`../tables/CD11b_microglia_intrinsic_genes.tsv`).

| Universe | Use it for | Reference list |
|---|---|---|
| `all` | Everything, comparable to the dura kit | `REFERENCE_all_15492_EnsemblIDs.txt` |
| `mg_intrinsic` | **Microglial biology (preferred)** | `REFERENCE_mg_intrinsic_13769_EnsemblIDs.txt` |

Always pair a list with the reference list from its **own** universe.

## PANTHER: `panther/`

| File | Contents |
|---|---|
| `{universe}_UP_p05` / `_DOWN_p05` | DESeq2 `~ treatment` p < 0.05 |
| `{universe}_UP_p01` / `_DOWN_p01` | p < 0.01 |
| `{universe}_RANKED_wald_stat.txt` | All genes with their Wald statistic (positive = higher in SPP1), for the *Statistical enrichment test* |
| `list_sizes.csv` | Counts |

List sizes:

| Universe | p < 0.01 up / down | p < 0.05 up / down |
|---|---|---|
| all | 275 / 62 | 634 / 260 |
| mg_intrinsic | 87 / 53 | 277 / 238 |

Most of the excess "up" genes in `all` are carry-over, which is why the intrinsic up list is much shorter.

**Suggested runs** (Mus musculus; Fisher test with FDR; GO BP complete, PANTHER Pathways, Reactome):
1. `mg_intrinsic_DOWN_p05` against the intrinsic reference. Expected themes: lipid/sterol metabolism (Insig1, Hmgcr, Ldlr, Cyp51), lysosome, DAM-lipid genes (Gpnmb, Apoe, Ch25h, Fabp5).
2. `mg_intrinsic_UP_p05` against the intrinsic reference. Includes the glucocorticoid/stress genes (Fkbp5, Tsc22d3, Sgk1), driven mostly by one SPP1 replicate.
3. `mg_intrinsic_RANKED_wald_stat.txt` with the enrichment test. Prefer this one: with 3 vs 3, a p-value cut-off list is very noisy.
4. Optional: `all_UP_p05`, to see the carry-over signal (ECM, angiogenesis) for comparison. Do not interpret it as microglial.

**PANTHER pathway diagrams to open:** *Cholesterol biosynthesis*, *Integrin signalling*, *Interferon-gamma signaling*, *Inflammation mediated by chemokine and cytokine signaling*. Comparing them side by side with the dura diagrams is the quickest cross-tissue view.

## IPA: `ipa/IPA_upload_SPP1_vs_NTC_5xFAD_CD11b_microglia.txt`

| Column | Meaning |
|---|---|
| `Ensembl_ID` | Identifier (mouse). Use this |
| `Symbol` | For reference only |
| `log2FC_unadjusted`, `pvalue_unadjusted`, `padj_unadjusted` | `~ treatment` (**use as the main observation**) |
| `log2FC_carryover_adjusted`, `pvalue_carryover_adjusted` | `~ carry_astro + treatment`. The covariate correlates r = 0.84 with treatment, so effects are inflated and unstable; use as a sensitivity check only |
| `microglia_intrinsic` | 1 = intrinsic (see above) |
| `baseMean` | Mean normalised count |

### Steps
1. Upload with Identifier = Ensembl, species Mouse.
2. Set one observation: `log2FC_unadjusted` as *Expr Log Ratio* and `pvalue_unadjusted` as *Expr p-value*.
3. Core Analysis settings:
   - cut-off p < 0.05 (~900 genes) or p < 0.01 (~340 genes);
   - **Reference set: User Dataset.**
4. To restrict to intrinsic genes, filter the dataset in IPA on `microglia_intrinsic = 1` before the analysis, or delete those rows in Excel and upload a second file. Then run a **Comparison Analysis** with the two dura observations, to see dura and microglia side by side.

**Upstream regulators worth checking:**
- **SPP1**: expected to be predicted inhibited.
- **SREBF2 / SREBF1, INSIG1, SCAP, 25-hydroxycholesterol, LXR (NR1H3)**: the sterol axis.
- **TREM2, APOE**: the lipid-associated microglia arm.
- **NR3C1 / dexamethasone**: the stress signature.
- **IFNB1 / IFNG**: expected **not** activated in microglia, unlike dura.

## Caveats
- **n = 3 vs 3, all male.** The smallest exact permutation p is 0.1. Only Spp1 is solid at the single-gene level.
- **Library quality differs by arm:** SPP1 has more unmapped reads and fewer assigned reads. About 20% of DE hits (pseudogenes and unnamed loci) track this; they are mostly outside the pathways above.
- **Brain microglia are not dural macrophages.** A negative result here (e.g. interferon) does not rule out the same change in dural myeloid cells.
