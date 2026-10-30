"""05: per-fraction PCA on VST (top 2000 variable genes), PC ~ covariate R^2, and figures."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ss = sample_sheet()
q = pd.read_csv(os.path.join(ROOT, "01_QC", "qc_metrics.tsv"), sep="\t", index_col=0)
pur = pd.read_csv(os.path.join(ROOT, "02_purity_celltyping", "purity_carryover_per_sample.tsv"), sep="\t", index_col=0)
sc = pd.read_csv(os.path.join(ROOT, "02_purity_celltyping", "panel_scores_z_within_fraction.tsv"), sep="\t", header=[0, 1], index_col=0)
rows = []
fig, axes = plt.subplots(1, 2, figsize=(10, 4.4))
for ax, frac in zip(axes, ["CD11b", "ACSA2"]):
    v = pd.read_csv(os.path.join(ROOT, "04_exploration", f"vst_{frac}.tsv.gz"), sep="\t", index_col=0)
    top = v.var(1).sort_values(ascending=False).index[:2000]
    X = v.loc[top].T; X = X - X.mean()
    U, S, Vt = np.linalg.svd(X.values, full_matrices=False)
    ve = S**2 / (S**2).sum(); pcs = pd.DataFrame(U * S, index=X.index, columns=[f"PC{i+1}" for i in range(len(S))])
    s = pcs.index
    cov = pd.DataFrame(index=s)
    cov["treatment"] = (ss.loc[s, "treatment"] == "Spp1").astype(float)
    for c in ["pct_unmapped", "pct_assigned_of_unique", "pct_mito", "pct_ribo", "pct_multi", "genes_ge5"]:
        cov[c] = q.loc[s, c]
    for c in ["microglial_RNA_pct_vs_CD11b", "astro_marker_pct_vs_ACSA2", "oligo_marker_pct_vs_ACSA2", "Ttr_cpm",
              "neutrophil_cpm_sum", "hemoglobin_cpm", "dissociation_IEG_log2cpm_mean"]:
        cov[c] = pur.loc[s, c]
    for k in ["L_Endothelial", "L_Pericyte_SMC", "L_Ependymal_CP", "L_Neutrophil", "L_Oligodendrocyte", "L_Astrocyte", "L_Microglia"]:
        if (frac, k) in sc.columns: cov["score_" + k] = sc[(frac, k)].loc[s]
    for i in range(4):
        pc = pcs.iloc[:, i]
        for c in cov:
            r = np.corrcoef(pc, cov[c])[0, 1]
            rows.append(dict(fraction=frac, PC=f"PC{i+1}", var_explained=ve[i], covariate=c, r=r, R2=r**2))
    for t, mk in [("NTC", "o"), ("Spp1", "o")]:
        m = ss.loc[s, "treatment"] == t
        ax.scatter(pcs.loc[m, "PC1"], pcs.loc[m, "PC2"], c=COLORS[t], s=70, marker=mk, edgecolor="k", label=t)
    for i in s: ax.annotate(i.split("_")[1], (pcs.loc[i, "PC1"], pcs.loc[i, "PC2"]), xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.set_xlabel(f"PC1 ({ve[0]*100:.0f}%)"); ax.set_ylabel(f"PC2 ({ve[1]*100:.0f}%)"); ax.set_title(f"{frac}+ (all male)")
    ax.legend(fontsize=8)
    pcs.iloc[:, :5].to_csv(os.path.join(ROOT, "04_exploration", f"pca_scores_{frac}.tsv"), sep="\t", float_format="%.3f")
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "figures", "fig04_pca_per_fraction.png"), dpi=150)
r = pd.DataFrame(rows)
r.to_csv(os.path.join(ROOT, "04_exploration", "pc_covariate_correlations.tsv"), sep="\t", index=False, float_format="%.3f")
pd.set_option("display.width", 200)
for (f, pc), g in r.groupby(["fraction", "PC"], sort=False):
    g = g.sort_values("R2", ascending=False).head(5)
    print(f, pc, f"{g.var_explained.iloc[0]*100:.1f}%", "; ".join(f"{a}={b:+.2f}" for a, b in zip(g.covariate, g.r)))
