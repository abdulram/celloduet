# Upload kits: PANTHER and QIAGEN IPA (dura, 5xFAD SPP1-targeted vs NTC)

All files come from `scripts/13c_pathway_links_and_kits.py`. Comparison: **5xFAD_SPP1 vs 5xFAD_NTC**, n = 6 vs 6 (3F/3M each). Two DESeq2 models:

| Model | Design | What it asks |
|---|---|---|
| `unadjusted` | `~ sex + treatment` | Everything that differs between SPP1 and NTC dura, including changes in cell composition |
| `fibro_adj` | `~ sex + fibroblast_stromal + treatment` | Changes beyond the lower fibroblast content of the SPP1 samples (step 12b, permutation p = 0.025) |

Run both models. A pathway found only in `unadjusted` is probably carried by fibroblast abundance. A pathway that survives in `fibro_adj` is more likely a change in cell state.

**Read this before interpreting anything.** No gene set passes FDR in either model, and only 5 genes pass genome-wide (in `fibro_adj` only). The input lists use a nominal p threshold, so on their own they contain many false positives. Use PANTHER/IPA to **generate and organise hypotheses**, not to confirm them.

---

## PANTHER (pantherdb.org): `panther/`

| File | Contents |
|---|---|
| `REFERENCE_expressed_genes_19622_EnsemblIDs.txt` | All genes tested by DESeq2, i.e. the expressed background. **Always use this as the reference list.** The whole mouse genome inflates enrichment for any tissue-expressed pathway. |
| `{model}_UP_p05_EnsemblIDs.txt`, `{model}_DOWN_p05_EnsemblIDs.txt` | Genes with p < 0.05, split by direction |
| `{model}_UP_p01_EnsemblIDs.txt`, `{model}_DOWN_p01_EnsemblIDs.txt` | Stricter lists with p < 0.01 |
| `{model}_RANKED_wald_stat.txt` | All genes with their DESeq2 Wald statistic (positive = higher in SPP1), for the enrichment test |
| `list_sizes.csv` | Number of genes in each list |

List sizes:

| Model | p < 0.01 up / down | p < 0.05 up / down |
|---|---|---|
| unadjusted | 129 / 189 | 616 / 720 |
| fibro_adj | 101 / 142 | 401 / 589 |

### A. Overrepresentation test (one run per list)
1. Go to pantherdb.org. Paste or upload a list (e.g. `unadjusted_DOWN_p05_EnsemblIDs.txt`) and set **Organism = Mus musculus**.
2. Choose **Statistical overrepresentation test**, then Submit.
3. **Reference list:** click *Change* and upload `REFERENCE_expressed_genes_19622_EnsemblIDs.txt` (organism Mus musculus).
4. **Annotation data set:** run several times, once each with:
   - *GO biological process complete* (detailed; GO-Slim BP is easier to read)
   - *PANTHER Pathways*: these have the clickable **pathway diagrams**
   - *Reactome pathways*
   - optionally *GO cellular component* or *GO molecular function*
5. **Test:** Fisher's exact. **Correction:** FDR. Then *Launch analysis*.
6. Export each result with *Export results → Table* and save it as e.g. `panther_out/unadjusted_DOWN_p05_GOBP.txt`.

Suggested order, from most to least informative:
1. `unadjusted_DOWN_p05`: the matrix/YAP-TAZ/TGF-β axis
2. `fibro_adj_UP_p05`: the interferon/eosinophil axis
3. `fibro_adj_DOWN_p05`: cholesterol/SREBP
4. `unadjusted_UP_p05`
5. The p < 0.01 versions, to see which terms survive a tighter list

### B. Enrichment test (ranked, no cut-off; closest to GSEA)
1. Upload `{model}_RANKED_wald_stat.txt` (two columns: Ensembl ID, Wald statistic) and set Mus musculus.
2. Choose **Statistical enrichment test**. Pick the annotation set as above.
3. In the results, a positive shift means a pathway trends higher in SPP1. This test uses all genes, so it avoids the arbitrary p < 0.05 cut-off. Prefer it when it disagrees with the overrepresentation test.

### C. PANTHER pathway diagrams
1. From an overrepresentation or enrichment result using *PANTHER Pathways*, click the pathway name. Useful candidates: *Integrin signalling pathway*, *TGF-beta signaling pathway*, *Interferon-gamma signaling pathway*, *JAK/STAT signaling pathway*, *Cholesterol biosynthesis*, *Inflammation mediated by chemokine and cytokine signaling*, *Wnt/Notch signaling*.
2. The diagram highlights the genes from your list. To colour by direction, run the UP list and the DOWN list separately, or upload the ranked file with the enrichment test.
3. You can also browse any pathway directly at pantherdb.org → *Pathways* and click **"Upload and display your data"** to overlay the ranked file.
4. Expected coverage: Integrin signalling shows Itgav, Fn1 and collagens (down, unadjusted). Interferon-gamma and JAK/STAT show Stat1, Irf1 and Irf9 (up, fibro_adj). Cholesterol biosynthesis shows Sqle, Msmo1 and Idi1 (down, fibro_adj).

---

## QIAGEN IPA: `ipa/IPA_upload_SPP1_vs_NTC_5xFAD_dura.txt`

Tab-separated, one row per tested gene (19,622 rows).

| Column | Meaning |
|---|---|
| `Ensembl_ID` | Mouse Ensembl gene ID. **Use this as the identifier**; Symbol is for your reference only |
| `Symbol` | Gene symbol |
| `log2FC_unadjusted`, `pvalue_unadjusted`, `padj_unadjusted` | Model 1 (positive = higher in SPP1) |
| `log2FC_fibroblast_adjusted`, `pvalue_fibroblast_adjusted`, `padj_fibroblast_adjusted` | Model 2 |
| `baseMean` | Mean normalised count |

### Steps
1. *New → Upload Dataset*. Format: flexible/tab-delimited. **Identifier type: Ensembl. Species: Mouse.** Array platform: *Not specified / RNA-seq*.
2. Column roles. **Observation 1** (`unadjusted`): `log2FC_unadjusted` = *Expr Log Ratio*, `pvalue_unadjusted` = *Expr p-value*. **Observation 2** (`fibroblast_adjusted`): the matching columns. Set all other columns to *Ignore*.
3. Check the mapping summary. About 90% of Ensembl IDs should map; unmapped IDs are mostly predicted genes (Gm*) and lncRNAs.
4. *New → Core Analysis → Expression Analysis*:
   - **Cut-offs:** Expr p-value < 0.05 gives about 1,300 genes (unadjusted) or about 990 (fibro_adj), which IPA handles well. Repeat with p < 0.01 (about 320 / 240 genes) as a sensitivity check. Do not use padj: almost nothing passes it.
   - **Reference set: User Dataset.** This is the expressed background and is very important. The default "Ingenuity Knowledge Base" background inflates every tissue-expressed pathway.
   - Species: Mouse. Tissues/cells: leave all, or restrict to relevant ones; restricting can drop real hits.
5. Run one Core Analysis per observation, then **Comparison Analysis** (*New → Comparison Analysis*) across the two. This shows directly which pathways and regulators depend on fibroblast content.

### What to look at in IPA
- **Canonical Pathways:** compare with PANTHER and the camera results below. In IPA the z-score has a direction; |z| ≥ 2 is IPA's usual call.
- **Upstream Regulators**, the most useful part here. Questions to ask:
  - Is **SPP1** itself predicted inhibited? That would be consistent with SPP1 targeting. Spp1 mRNA is −1.27 log2 (p = 0.018), mostly in males.
  - Are **TGFB1**, **YAP1/WWTR1** and **SMAD3** predicted inhibited in `unadjusted` and not in `fibro_adj`?
  - Are **IFNA/IFNB1/IFNG**, **IRF7**, **STAT1** activated?
  - Is **SREBF2** inhibited, or **25-hydroxycholesterol** activated?
  - Is **CSF2** or **IL5** activated? That would fit the eosinophil signal.
- **Causal / Mechanistic Networks:** check whether IPA can connect SPP1 → CD44/ITGAV → TGF-β/YAP; this is the hypothesis in `../13_pathway_links/REPORT.md`.
- **Diseases & Functions:** expect fibrosis, cell migration and antiviral response terms. These are descriptive only.

### Caveats specific to IPA
- With n = 6 per group and no FDR-significant genes, every z-score is built on nominal hits. Trust patterns that agree across PANTHER, IPA and camera, not single IPA calls.
- IPA's knowledge base is mostly human and cancer/immune biology. Dura fibroblast and arachnoid barrier biology are poorly represented, so Cldn11/Prg4/Matn4 (our most robust single genes) will rarely appear in any pathway.
- The `unadjusted` results partly measure fewer fibroblasts and more myeloid cells in the SPP1 dura. IPA cannot tell that apart from regulation inside a cell. The `fibro_adj` model helps, but monocyte and neutrophil content are not adjusted for (see REPORT.md).

---

## Where our own pathway results are

These were run locally because PANTHER/IPA could not be reached from the analysis environment:

- `../tables/spp1_pathway_camera_both_models.tsv`: limma camera on 50 Hallmark sets + curated + mechanistic modules, both models. Use the conservative columns (`p_conservative`, `FDR_conservative`).
- `../13_pathway_links/`: link map, mechanistic gene table, and IFN/eosinophil composition checks. Interpretation is in `REPORT.md`.

When your PANTHER/IPA results are back, compare them with those tables. Agreement across three methods is the strongest evidence this design can give.

**GO enrichment was not run locally.** The only mouse GO/Reactome source reachable here was an R data file from a third-party GitHub package (msigdbdf), and loading it was blocked by the session's safety policy. PANTHER covers GO properly, so use its GO BP output for that.
