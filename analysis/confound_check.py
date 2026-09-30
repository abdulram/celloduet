"""Is the 5xFAD-vs-WT signature biological or technical?

Run after run_analysis.py. Adds:
  * a per-sample score of the base-model genotype signature vs. library mapping profile
  * genotype DE counts under single-covariate adjustments
"""
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

import run_analysis as ra  # re-uses loaded data, styling and helpers (re-runs the base pipeline)

warnings.filterwarnings("ignore")
OUT, FIG, meta, qc, vst, cf = ra.OUT, ra.FIG, ra.meta, ra.qc, ra.vst, ra.cf

base = pd.read_csv(os.path.join(OUT, "DE_5xFAD_NTC_vs_WT_NTC_base.csv"), index_col=0)
sig = base[(base.padj < 0.01) & (base.baseMean > 50)]
up, dn = sig.index[sig.log2FoldChange > 0.5], sig.index[sig.log2FoldChange < -0.5]
z = vst.sub(vst.mean(1), axis=0).div(vst.std(1), axis=0)
score = z.loc[up].mean() - z.loc[dn].mean()
qc["genotype_signature_score"] = score

fig, axes = plt.subplots(1, 2, figsize=(12.5, 5))
for ax, col, lab in [(axes[0], "pct_multi", "% multi-mapped reads"), (axes[1], "pct_mito", "% mitochondrial reads")]:
    ra.scatter_groups(ax, qc[col], score)
    r = np.corrcoef(qc[col], score)[0, 1]
    ax.set_xlabel(lab)
    ax.set_ylabel("5xFAD-vs-WT signature score\n(mean z of up genes − down genes)")
    ax.set_title(f"Signature vs {lab} (Pearson r = {r:.2f})")
ra.group_legend(axes[0], loc="upper left")
fig.suptitle("The 'genotype' signature follows a library-level technical axis: WT 616 (processed like the 5xFAD libraries) scores like 5xFAD",
             fontweight="bold", fontsize=10)
fig.tight_layout()
ra.savefig(fig, "10_genotype_signature_vs_technical.png")

# Biotype / gene-class make-up of the signature
cls = pd.Series("other", index=sig.index)
names = sig.gene_name.astype(str)
cls[sig.biotype.str.contains("pseudogene", na=False).values] = "pseudogene"
cls[names.str.startswith("mt-").values] = "mitochondrial"
cls[names.str.match(r"^(Rpl|Rps)\d").values] = "ribosomal protein"
cls[names.str.match(r"^(Myh[1-8]|Acta1|Ckm|Tnn[ict]\d|Myl\d|Mylpf|Csrp3|Myoz\d|Smpx|Lmod\d)$").values] = "skeletal/cardiac muscle"
make_up = pd.crosstab(cls, np.where(sig.log2FoldChange > 0, "up_in_5xFAD", "down_in_5xFAD"))
make_up.to_csv(os.path.join(OUT, "genotype_signature_gene_classes.csv"))
print(make_up)

# Single-covariate sensitivity models for the genotype contrast
rows = []
for cov in [None, "pct_multi", "pct_mito", "mod_osteoclast", "mod_choroid_plexus", "mod_olfactory_mucosa", "mod_skeletal_muscle"]:
    md = ra.design_meta.copy()
    design = "~sex + group"
    if cov:
        md[cov] = ((qc[cov] - qc[cov].mean()) / qc[cov].std()).values
        design = f"~sex + {cov} + group"
    dds = DeseqDataSet(counts=cf.T, metadata=md, design=design, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["group", "5xFAD_NTC", "WT_NTC"], quiet=True)
    st.summary()
    res = st.results_df
    rows.append({"covariate": cov or "(none)", "genotype_DE_padj<0.05": int((res.padj < 0.05).sum()),
                 "rho(covariate, genotype)": (np.corrcoef(qc[cov], (meta.genotype == "5xFAD").astype(int))[0, 1] if cov else np.nan)})
sens = pd.DataFrame(rows)
sens.to_csv(os.path.join(OUT, "genotype_DE_sensitivity.csv"), index=False)
print(sens.to_string())
qc.to_csv(os.path.join(OUT, "qc_and_composition.csv"))
