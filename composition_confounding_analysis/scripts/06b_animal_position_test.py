"""06b - Single-animal position test on the technical axis (technical common support has n_WT = 1).

On % mitochondrial reads (the technical metric most associated with genotype, step 03) WT animal 616 is the ONLY
WT animal inside the 5xFAD range; on % multi-mapped reads 623 is intermediate. For every baseline genotype gene
(step 04, padj<0.05) we place each WT animal between the mean of the remaining WT animals (0) and the 5xFAD_NTC
mean (1):  t = (x_animal - mean(other WT)) / (mean(5xFAD_NTC) - mean(other WT)).
If the baseline signature is genotype-driven, every WT animal - including 616 - should sit near 0.
If it follows the technical axis, 616 should sit near 1 while technically typical WT animals sit near 0.
Leave-one-out positions of the other WT animals give the null distribution.
"""
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CA = os.path.dirname(HERE)
OUT = os.path.join(CA, "06_common_support_matching")
TAB = os.path.join(CA, "tables")
exec(open(os.path.join(HERE, "common.py")).read())

meta = pd.read_csv(os.path.join(TAB, "sample_metadata_with_composition.tsv"), sep="\t", index_col=0)
vst = pd.read_csv(os.path.join(CA, "01_QC", "vst_blind.tsv.gz"), sep="\t", index_col=0)
base = pd.read_csv(os.path.join(TAB, "DE_baseline_5xFAD_vs_WT.tsv"), sep="\t")
sig = base.gene_id[base.padj < 0.05].values
WT = meta.index[meta.group == "WT_NTC"]
FX = meta.index[meta.group == "5xFAD_NTC"]
pos = {}
for a in WT:
    rest = WT.drop(a)
    d = vst.loc[sig, FX].mean(axis=1) - vst.loc[sig, rest].mean(axis=1)
    pos[str(meta.animal[a])] = (vst.loc[sig, a] - vst.loc[sig, rest].mean(axis=1)) / d
pos = pd.DataFrame(pos)
summ = pd.DataFrame({"animal": pos.columns, "pct_mito": [meta.set_index(meta.animal.astype(str)).loc[c, "pct_mito"] for c in pos.columns],
                     "pct_multi": [meta.set_index(meta.animal.astype(str)).loc[c, "pct_multi"] for c in pos.columns],
                     "median_position": pos.median().values, "frac_genes_5xFAD_like_(t>0.5)": (pos > 0.5).mean().values})
summ.to_csv(os.path.join(TAB, "technical_axis_animal_position_test.tsv"), sep="\t", index=False)
pos.describe().to_csv(os.path.join(TAB, "technical_axis_animal_position_distribution.tsv"), sep="\t")
print(summ.round(3).to_string())

fig, ax = plt.subplots(figsize=(9, 5))
bins = np.linspace(-1.5, 2.5, 41)
for c in pos.columns:
    hi = c == "616"
    ax.hist(pos[c].clip(-1.5, 2.5), bins=bins, histtype="stepfilled" if hi else "step", alpha=0.55 if hi else 0.9,
            color="#eb6834" if hi else ("#4a3aa7" if c == "623" else "#2a78d6"), lw=1.4,
            label=f"WT {c} (mito {summ.set_index('animal').loc[c, 'pct_mito']:.1f}%, multi {summ.set_index('animal').loc[c, 'pct_multi']:.1f}%)")
ax.axvline(0, ls="--", color="#2a78d6")
ax.axvline(1, ls="--", color="#eb6834")
ax.set_xlabel("position of the WT animal between other WT (0) and 5xFAD_NTC (1)")
ax.set_ylabel(f"genes (of {len(sig)} baseline genotype genes)")
ax.set_title("Technical common support: where does each WT animal fall on the baseline 'genotype' genes?")
ax.legend(fontsize=8)
save(fig, os.path.join(OUT, "Fig10c_animal_position_test.png"))
