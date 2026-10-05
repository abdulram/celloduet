"""03 - Composition PCA, transcriptome-PC vs composition, genotype vs composition, COMMON SUPPORT.

Composition matrix = validated marker-signature scores (step 02), each standardised across the 18 animals.
(Signature scores are not proportions, so no CLR is applied; deconvolution proportions are used only as a
cross-check, CLR-transformed.) PCA -> Composition_PC1..3 (NOT transcriptome PCs).

Primary genotype contrast: WT_NTC (n=6) vs 5xFAD_NTC (n=6). 5xFAD_SPP1 animals are shown for context.
Technical axis (% mito, % multi-mapped; from step 01) is analysed alongside composition because it was the
dominant genotype-associated factor found earlier in this project.

Common support per axis: overlap of the WT and 5xFAD ranges, and AUC (probability that a random 5xFAD animal
has a higher value than a random WT animal; 0.5 = full overlap, 0 or 1 = complete separation).
"""
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr
import statsmodels.formula.api as smf

HERE = os.path.dirname(os.path.abspath(__file__))
CA = os.path.dirname(HERE)
OUT = os.path.join(CA, "03_composition_PCA")
TAB = os.path.join(CA, "tables")
exec(open(os.path.join(HERE, "common.py")).read())

meta = pd.read_csv(os.path.join(TAB, "sample_metadata_qc_pcs.tsv"), sep="\t", index_col=0)
S = pd.read_csv(os.path.join(TAB, "composition_signature_scores.tsv"), sep="\t", index_col=0).loc[meta.index]
Sall = pd.read_csv(os.path.join(TAB, "composition_signature_scores_all.tsv"), sep="\t", index_col=0).loc[meta.index]
dec = pd.read_csv(os.path.join(TAB, "deconvolution_proportions.tsv"), sep="\t", index_col=0).loc[meta.index]
dval = pd.read_csv(os.path.join(TAB, "deconvolution_validation.tsv"), sep="\t")
txvar = pd.read_csv(os.path.join(TAB, "transcriptome_PC_variance.tsv"), sep="\t")

# ---------------------------------------------------------------- composition PCA
Z = (S - S.mean()) / S.std()
U, sv, Vt = np.linalg.svd(Z.values, full_matrices=False)
cvar = sv**2 / (sv**2).sum() * 100
cpcs = pd.DataFrame(U[:, :3] * sv[:3], index=meta.index, columns=["Composition_PC1", "Composition_PC2", "Composition_PC3"])
cload = pd.DataFrame(Vt[:3].T, index=S.columns, columns=cpcs.columns)
# sign convention: PC1 positive = more fibroblast/stromal signal
for c in cpcs:
    lead = cload[c].abs().idxmax()
    if c == "Composition_PC1" and "fibroblast_stromal" in cload.index:
        lead = "fibroblast_stromal"
    if cload.loc[lead, c] < 0:
        cpcs[c] *= -1
        cload[c] *= -1
cload.to_csv(os.path.join(TAB, "composition_PC_loadings.tsv"), sep="\t")
pd.DataFrame({"axis": cpcs.columns, "pct_composition_variance": cvar[:3].round(1)}).to_csv(
    os.path.join(TAB, "composition_PC_variance.tsv"), sep="\t", index=False)
meta = meta.join(cpcs)
meta.to_csv(os.path.join(TAB, "sample_metadata_with_composition.tsv"), sep="\t")
print("Composition PC variance %:", cvar[:5].round(1))
print(cload.round(2).to_string())

# CLR of reliable deconvolution classes (cross-check only)
rel = dval.reference_class[dval.reliable].tolist()
D = dec[rel] + 1e-3
clr = np.log(D).sub(np.log(D).mean(axis=1), axis=0)
clr.to_csv(os.path.join(TAB, "deconvolution_reliable_classes_CLR.tsv"), sep="\t")

# ---------------------------------------------------------------- transcriptome PCs vs composition
tx = meta[[f"Tx_PC{i}" for i in range(1, 6)]]
feats = pd.concat([cpcs, S, meta[["pct_mito", "pct_multi", "pct_unique"]]], axis=1)
rho = pd.DataFrame({t: [spearmanr(tx[t], feats[f])[0] for f in feats] for t in tx}, index=feats.columns)
pv = pd.DataFrame({t: [spearmanr(tx[t], feats[f])[1] for f in feats] for t in tx}, index=feats.columns)
rho.to_csv(os.path.join(TAB, "txPC_vs_composition_spearman_rho.tsv"), sep="\t")
pv.to_csv(os.path.join(TAB, "txPC_vs_composition_spearman_p.tsv"), sep="\t")

# variance of each Tx PC explained (R^2) by composition PCs 1-3 jointly vs technical axis jointly
r2 = []
for t in tx:
    d = meta.assign(y=tx[t])
    r_comp = smf.ols("y ~ Composition_PC1 + Composition_PC2 + Composition_PC3", d).fit().rsquared
    r_tech = smf.ols("y ~ pct_mito + pct_multi", d).fit().rsquared
    r_both = smf.ols("y ~ Composition_PC1 + Composition_PC2 + Composition_PC3 + pct_mito + pct_multi", d).fit().rsquared
    r_geno = smf.ols("y ~ C(group)", d).fit().rsquared
    r2.append({"Tx_PC": t, "pct_transcriptome_variance": float(txvar.set_index("PC").loc[t, "pct_variance"]),
               "R2_composition_PC1to3": r_comp, "R2_technical_mito_multimap": r_tech, "R2_composition_plus_technical": r_both,
               "R2_group": r_geno})
r2 = pd.DataFrame(r2)
r2["weighted_share_composition"] = r2.pct_transcriptome_variance * r2.R2_composition_PC1to3
r2.to_csv(os.path.join(TAB, "txPC_variance_explained_by_composition.tsv"), sep="\t", index=False)
print(r2.round(2).to_string())

# ---------------------------------------------------------------- genotype vs composition (+ technical)
prim = meta[meta.group.isin(["WT_NTC", "5xFAD_NTC"])].copy()
prim["g5"] = (prim.genotype == "5xFAD").astype(int)
allg = meta.copy()
allg["g5"] = (allg.genotype == "5xFAD").astype(int)
rows = []
measures = list(cpcs.columns) + list(S.columns) + ["pct_mito", "pct_multi"]
F = pd.concat([cpcs, S, meta[["pct_mito", "pct_multi"]]], axis=1)
for m in measures:
    x = F.loc[prim.index, m]
    wt, fx = x[prim.genotype == "WT"], x[prim.genotype == "5xFAD"]
    auc = mannwhitneyu(fx, wt).statistic / (len(fx) * len(wt))
    fit = smf.ols("y ~ C(sex) + g5", prim.assign(y=x)).fit()
    pooled = np.sqrt(((len(wt) - 1) * wt.var() + (len(fx) - 1) * fx.var()) / (len(wt) + len(fx) - 2))
    lo, hi = max(wt.min(), fx.min()), min(wt.max(), fx.max())
    span = max(wt.max(), fx.max()) - min(wt.min(), fx.min())
    xa = F.loc[allg.index, m]
    auc_all = mannwhitneyu(xa[allg.g5 == 1], xa[allg.g5 == 0]).statistic / ((allg.g5 == 1).sum() * (allg.g5 == 0).sum())
    rows.append({"measure": m, "kind": "composition PC" if m.startswith("Composition") else ("technical" if m.startswith("pct") else "signature"),
                 "WT_mean": wt.mean(), "5xFAD_NTC_mean": fx.mean(), "hedges_g": (fx.mean() - wt.mean()) / pooled * (1 - 3 / (4 * 10 - 1)),
                 "AUC_5xFAD_gt_WT": auc, "MWU_p": mannwhitneyu(fx, wt).pvalue, "OLS_sex_adj_beta": fit.params["g5"], "OLS_p": fit.pvalues["g5"],
                 "overlap_fraction_of_range": max(0, hi - lo) / span if span else 1.0,
                 "AUC_all12_5xFAD_vs_WT": auc_all,
                 "common_support": "none (complete separation)" if (auc in (0, 1)) else ("poor" if (auc <= 0.1 or auc >= 0.9) else "adequate")})
gc = pd.DataFrame(rows)
gc.to_csv(os.path.join(TAB, "genotype_vs_composition.tsv"), sep="\t", index=False)
print(gc[["measure", "hedges_g", "AUC_5xFAD_gt_WT", "MWU_p", "overlap_fraction_of_range", "common_support"]].round(3).to_string())

# ---------------------------------------------------------------- figures
# Fig 4: composition PCA
fig, axes = plt.subplots(1, 2, figsize=(15, 6))
ax = axes[0]
for g in GROUP_ORDER:
    for sx, mk in SEX_MARKER.items():
        idx = meta.index[(meta.group == g) & (meta.sex == sx)]
        ax.scatter(meta.Composition_PC1[idx], meta.Composition_PC2[idx], c=GROUP_COLORS[g], marker=mk, s=80, edgecolor=SURFACE, lw=1.5, zorder=3)
for l in meta.index:
    ax.annotate(str(meta.animal[l]), (meta.Composition_PC1[l], meta.Composition_PC2[l]), xytext=(5, 4), textcoords="offset points", fontsize=7.5, color=INK2)
ax.set_xlabel(f"Composition_PC1 ({cvar[0]:.1f}% of composition variance)")
ax.set_ylabel(f"Composition_PC2 ({cvar[1]:.1f}%)")
ax.set_title("Composition PCA (validated cell/anatomical signatures)")
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], marker="s", ls="", color=c, ms=9, label=g.replace("_", " ")) for g, c in GROUP_COLORS.items()]
          + [Line2D([], [], marker=m, ls="", color=INK2, ms=8, label={"F": "female", "M": "male"}[s]) for s, m in SEX_MARKER.items()], fontsize=8)
ax = axes[1]
L = cload[["Composition_PC1", "Composition_PC2"]]
ax.axhline(0, color=INK2, lw=0.6)
ax.axvline(0, color=INK2, lw=0.6)
ax.scatter(L.iloc[:, 0], L.iloc[:, 1], s=30, c="#2a78d6")
for s_, r in L.iterrows():
    ax.annotate(s_, (r.iloc[0], r.iloc[1]), xytext=(4, 3), textcoords="offset points", fontsize=8)
ax.set_xlabel("loading on Composition_PC1")
ax.set_ylabel("loading on Composition_PC2")
ax.set_title("What the composition axes represent")
fig.tight_layout()
save(fig, os.path.join(OUT, "Fig4_composition_PCA.png"))

# Fig 5: Tx PCs vs composition (heatmap + key scatters)
fig, axes = plt.subplots(1, 2, figsize=(16, 8.5), gridspec_kw={"width_ratios": [1.1, 1]})
ax = axes[0]
im = ax.imshow(rho.values, cmap=DIV, vmin=-1, vmax=1, aspect="auto")
ax.set_yticks(range(len(rho)), rho.index, fontsize=8)
ax.set_xticks(range(5), [f"Tx_PC{i}\n{txvar.pct_variance[i-1]:.0f}%" for i in range(1, 6)], fontsize=8)
for i in range(rho.shape[0]):
    for j in range(rho.shape[1]):
        v = rho.values[i, j]
        ax.text(j, i, f"{v:.2f}{'*' if pv.values[i, j] < 0.05 else ''}", ha="center", va="center", fontsize=6.5, color="white" if abs(v) > 0.6 else INK)
ax.grid(False)
fig.colorbar(im, ax=ax, shrink=0.5, label="Spearman ρ (* p<0.05)")
ax.set_title("Transcriptome PCs vs composition & technical measures")
ax = axes[1]
ax.axis("off")
sub = fig.add_gridspec(2, 2, left=0.57, right=0.98, top=0.92, bottom=0.08, hspace=0.4, wspace=0.35)
pairs = [("Tx_PC1", "Composition_PC1"), ("Tx_PC1", r.measure if False else "olfactory_nasal"), ("Tx_PC2", "pct_mito"), ("Tx_PC2", "Composition_PC1")]
for k, (a, b) in enumerate(pairs):
    axx = fig.add_subplot(sub[k // 2, k % 2])
    for g in GROUP_ORDER:
        idx = meta.index[meta.group == g]
        axx.scatter(feats.loc[idx, b], tx.loc[idx, a], c=GROUP_COLORS[g], s=30, edgecolor=SURFACE, lw=1)
    axx.set_xlabel(b, fontsize=8)
    axx.set_ylabel(a, fontsize=8)
    axx.set_title(f"ρ={rho.loc[b, a]:.2f}", fontsize=9)
save(fig, os.path.join(OUT, "Fig5_txPC_vs_composition.png"))

# Fig 6: common support
show = ["Composition_PC1", "Composition_PC2", "Composition_PC3", "pct_mito", "pct_multi"] + list(
    gc[gc.kind == "signature"].sort_values("AUC_5xFAD_gt_WT", key=lambda s: -abs(s - 0.5)).measure.head(7))
fig, axes = plt.subplots(3, 4, figsize=(16, 11))
for ax, m in zip(axes.flat, show):
    strip(ax, F[m], meta, label_animals=True)
    row = gc.set_index("measure").loc[m]
    ax.set_title(f"{m}\nWT vs 5xFAD-NTC AUC={row.AUC_5xFAD_gt_WT:.2f} ({row.common_support})", fontsize=8.5,
                 color="#d03b3b" if row.common_support != "adequate" else INK)
fig.suptitle("Common support: do WT and 5xFAD animals overlap on composition and technical axes? (AUC 0.5 = full overlap; 0/1 = separation)", fontweight="bold")
fig.tight_layout()
save(fig, os.path.join(OUT, "Fig6_common_support.png"))

fig, ax = plt.subplots(figsize=(7.5, 6))
for g in GROUP_ORDER:
    for sx, mk in SEX_MARKER.items():
        idx = meta.index[(meta.group == g) & (meta.sex == sx)]
        ax.scatter(meta.Composition_PC1[idx], meta.pct_mito[idx], c=GROUP_COLORS[g], marker=mk, s=80, edgecolor=SURFACE, lw=1.5)
for l in meta.index:
    ax.annotate(str(meta.animal[l]), (meta.Composition_PC1[l], meta.pct_mito[l]), xytext=(5, 3), textcoords="offset points", fontsize=7.5, color=INK2)
ax.set_xlabel("Composition_PC1")
ax.set_ylabel("% mitochondrial reads (technical axis)")
ax.set_title("Joint common support: composition vs technical axis")
save(fig, os.path.join(OUT, "common_support_compPC1_vs_mito.png"))
print("Done 03")
