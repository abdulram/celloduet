"""14a - Flag microglia-intrinsic genes in the CD11b+ fraction (for the SPP1 vs NTC pathway analysis, step 14).

SPP1 CD11b sorts carry more non-microglial RNA than NTC: astrocyte markers 1.5% vs 0.95%, BAM 1.8% vs 0.35%,
endothelial ~1.3% vs ~0.7%. Any gene made mostly by those cells shifts with treatment for compositional reasons.
A gene is called microglia-intrinsic when all three hold:
  1. the measured glial carry-over difference cannot shift it by more than 5% (|carryover_pred_log2FC| from step 06
     < 0.07; step 06 models CD11b cpm = intrinsic + carry-over share x ACSA-2 cpm);
  2. brain reference BAM / microglia expression < 5 (a 1.5-point BAM-share rise adds < ~8%);
  3. max(endothelial, pericyte, smooth muscle, VLMC fibroblast, neutrophil, monocyte) / microglia < 10.
Genes absent from the reference are kept only if rule 1 holds and they are not obviously from other cell types
(they are flagged 'no_reference').
Output: tables/CD11b_microglia_intrinsic_genes.tsv
"""
import os
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
de = pd.read_csv(os.path.join(ROOT, "05_DE_CD11b", "DE_Spp1_vs_NTC_CD11b_annotated.tsv"), sep="\t", index_col=0)
ref = pd.read_csv(os.path.join(ROOT, "data", "reference", "brain_pseudobulk_cp10k.csv.gz"), index_col=0).reindex(de.index)
eps = 0.01
mg = ref.microglia + eps
bam = (ref.BAM + eps) / mg
vasc = (ref[["endothelial", "pericyte", "smooth_muscle", "VLMC_fibroblast", "neutrophil", "monocyte"]].max(1) + eps) / mg
out = pd.DataFrame({"gene": de.gene, "baseMean": de.baseMean,
                    "ACSA2_over_CD11b": de.ratio_other_over_this,
                    "glial_carryover_pred_log2FC": de.carryover_pred_log2FC, "ref_BAM_over_microglia": bam,
                    "ref_vascular_blood_over_microglia": vasc})
r1 = out.glial_carryover_pred_log2FC.abs().fillna(np.inf) < 0.07
r2 = out.ref_BAM_over_microglia.fillna(0) < 5
r3 = out.ref_vascular_blood_over_microglia.fillna(0) < 10
out["no_reference"] = ref.microglia.isna()
out["microglia_intrinsic"] = r1 & r2 & r3
out["reason_excluded"] = np.select([~r1, ~r2, ~r3], ["glial carry-over (>5% predicted shift)", "BAM-enriched", "vascular/blood-enriched"], "")
out.to_csv(os.path.join(ROOT, "tables", "CD11b_microglia_intrinsic_genes.tsv"), sep="\t", float_format="%.4g")
print(out.microglia_intrinsic.value_counts(), out.reason_excluded.value_counts(), sep="\n")
for g in ["Spp1", "Gpnmb", "Apoe", "Isg15", "Rsad2", "Ifit3", "Irf7", "Cd44", "Itgav", "Wwtr1", "Fn1", "Insig1", "Sqle", "Gpr37l1", "Pf4", "Mrc1"]:
    print(g, out.loc[out.gene == g, ["microglia_intrinsic", "reason_excluded"]].values.tolist())
