"""Does SPP1 treatment change 5xFAD-associated changes in dura fibroblasts?

Uses only the fibroblast-specific genes from fibroblast_analysis.py (fibroblast-normalised).
  Selection bias guard: the 5xFAD signature is defined as WT vs ALL 12 5xFAD animals (SPP1 + NTC pooled), so the
  selection is symmetric in SPP1/NTC and the later SPP1-vs-NTC test cannot be inflated by regression to the mean
  (defining it as WT vs NTC and then testing SPP1 vs NTC would be biased toward "rescue").
  Step 1  5xFAD (all) vs WT: which fibroblast gene differences are real rather than the technical axis?
          - DESeq2 with and without the technical covariate (% mito)
          - "616 test": WT animal 616 was processed like the 5xFAD libraries (same mapping / mito profile).
            For each gene, t616 = (616 - WT_other) / (5xFAD_NTC - WT_other). t ~ 0: 616 looks WT (difference
            follows genotype). t ~ 1: 616 looks 5xFAD (difference follows the technical axis).
  Step 2  Robust 5xFAD fibroblast signature (genes that pass the technical checks).
  Step 3  Does SPP1 move 5xFAD fibroblasts along that signature (towards WT = rescue, or further = worsening)?
          - per-sample signature score, SPP1 vs NTC
          - genome-wide (fibroblast genes) concordance of SPP1 effect with the 5xFAD effect
          Significance by exact sample-label permutation (SPP1/NTC relabelled within sex: 20 x 20 = 400).
Outputs: results/fibroblast_5xfad_spp1/
"""
import itertools
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIB = os.path.join(ROOT, "results", "fibroblast")
OUT = os.path.join(ROOT, "results", "fibroblast_5xfad_spp1")
os.makedirs(OUT, exist_ok=True)

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
GROUP_COLORS = {"WT_NTC": "#2a78d6", "5xFAD_NTC": "#eb6834", "5xFAD_SPP1": "#1baf7a"}
GROUP_ORDER = list(GROUP_COLORS)
SEX_MARKER = {"F": "o", "M": "^"}
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold", "legend.frameon": False,
})


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=160, bbox_inches="tight")
    plt.close(fig)


qc = pd.read_csv(os.path.join(ROOT, "results", "qc_and_composition.csv"), index_col=0)
meta = qc[["sample", "genotype", "sex", "guide", "group", "pct_mito", "pct_multi"]].copy()
meta["animal"] = meta["sample"].str.split("_").str[0]
vst = pd.read_csv(os.path.join(FIB, "fibroblast_normalised_vst.csv"), index_col=0)[meta.index]
genes = pd.read_csv(os.path.join(FIB, "fibroblast_specific_genes.csv"), index_col=0)
name = genes.gene_name.astype(str)
de_g = pd.read_csv(os.path.join(FIB, "DE_5xFAD_NTC_vs_WT_NTC_fibnorm.csv"), index_col=0)
de_gm = pd.read_csv(os.path.join(FIB, "DE_5xFAD_NTC_vs_WT_NTC_fibnorm_mitoadj.csv"), index_col=0)
de_s = pd.read_csv(os.path.join(FIB, "DE_5xFAD_SPP1_vs_5xFAD_NTC_fibnorm.csv"), index_col=0)
de_sm = pd.read_csv(os.path.join(FIB, "DE_5xFAD_SPP1_vs_5xFAD_NTC_fibnorm_mitoadj.csv"), index_col=0)

lib = {a: meta.index[meta.animal == a][0] for a in meta.animal}
WT_ALL = meta.index[meta.group == "WT_NTC"]
WT_616 = lib["616"]
WT_OTHER = WT_ALL.drop(WT_616)
NTC = meta.index[meta.group == "5xFAD_NTC"]
SPP = meta.index[meta.group == "5xFAD_SPP1"]
FIVE = NTC.union(SPP)

# ---------------------------------------------------------------- Step 1: genotype differences vs technical axis
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

raw = pd.read_csv(os.path.join(ROOT, "data", "VT3YXB-expression-matrix.tsv"), sep="\t").set_index("gene_id")
counts = raw.loc[genes.index, [f"{l}_count" for l in meta.index]].astype(int)
counts.columns = meta.index
md = meta[["sex"]].copy()
md["genotype"] = pd.Categorical(meta.genotype, categories=["WT", "5xFAD"])
md["pct_mito"] = ((meta.pct_mito - meta.pct_mito.mean()) / meta.pct_mito.std()).values


def deseq_genotype(design):
    dds = DeseqDataSet(counts=counts.T, metadata=md, design=design, refit_cooks=True, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["genotype", "5xFAD", "WT"], quiet=True)
    st.summary()
    return st.results_df


de_g = deseq_genotype("~sex + genotype")
de_gm = deseq_genotype("~sex + pct_mito + genotype")

g = pd.DataFrame({"gene": name})
g["LFC_5xFAD_vs_WT"] = de_g.log2FoldChange.reindex(g.index)
g["padj_5xFAD_vs_WT"] = de_g.padj.reindex(g.index)
g["LFC_5xFAD_vs_WT_mitoadj"] = de_gm.log2FoldChange.reindex(g.index)
g["padj_5xFAD_vs_WT_mitoadj"] = de_gm.padj.reindex(g.index)
diff = vst.loc[g.index, FIVE].mean(axis=1) - vst.loc[g.index, WT_OTHER].mean(axis=1)
g["t616"] = (vst.loc[g.index, WT_616] - vst.loc[g.index, WT_OTHER].mean(axis=1)) / diff.replace(0, np.nan)
g["LFC_SPP1_vs_NTC"] = de_s.log2FoldChange.reindex(g.index)
g["padj_SPP1_vs_NTC"] = de_s.padj.reindex(g.index)
g["LFC_SPP1_vs_NTC_mitoadj"] = de_sm.log2FoldChange.reindex(g.index)
g["padj_SPP1_vs_NTC_mitoadj"] = de_sm.padj.reindex(g.index)

cand = g[g.padj_5xFAD_vs_WT < 0.05]
robust = cand[(cand.t616 < 0.5)]
g["genotype_call"] = "none"
g.loc[cand.index, "genotype_call"] = "follows technical axis (616 looks 5xFAD)"
g.loc[robust.index, "genotype_call"] = "candidate genotype effect (616 looks WT)"
g.loc[g.padj_5xFAD_vs_WT_mitoadj < 0.1, "genotype_call"] = "genotype effect (survives % mito adjustment)"
g.sort_values("padj_5xFAD_vs_WT").to_csv(os.path.join(OUT, "fibroblast_genes_genotype_and_spp1.csv"))
print(g.genotype_call.value_counts())
print("t616 among genotype-DE genes: median", round(cand.t616.median(), 2), " IQR",
      cand.t616.quantile([0.25, 0.75]).round(2).tolist())

# Null expectation for t616: if 616 were an ordinary WT animal, how often would a single WT animal look
# 5xFAD-like? Leave-one-out over the other WT animals (each treated as the "test" animal).
loo = []
for a in WT_OTHER:
    rest = WT_OTHER.drop(a)
    d = vst.loc[cand.index, FIVE].mean(axis=1) - vst.loc[cand.index, rest].mean(axis=1)
    t = (vst.loc[cand.index, a] - vst.loc[cand.index, rest].mean(axis=1)) / d
    loo.append(pd.Series(t.values, index=cand.index, name=meta.animal[a]))
loo = pd.concat(loo, axis=1)
loo.to_csv(os.path.join(OUT, "t_leave_one_out_WT.csv"))
print("t for ordinary WT animals (LOO) medians:", loo.median().round(2).to_dict())

# ---------------------------------------------------------------- Step 2/3: signature scores + SPP1
def z(df):
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1).replace(0, np.nan), axis=0)


def signature_score(idx, lfc):
    """Mean z of up-in-5xFAD genes minus mean z of down genes (higher = more 5xFAD-like)."""
    zz = z(vst.loc[idx])
    up, dn = idx[lfc.loc[idx] > 0], idx[lfc.loc[idx] < 0]
    s = pd.Series(0.0, index=vst.columns)
    if len(up):
        s += zz.loc[up].mean()
    if len(dn):
        s -= zz.loc[dn].mean()
    return s


SIGS = {
    "all genotype-DE fibroblast genes (5xFAD vs WT padj<0.05)": cand.index,
    "technically robust candidates (padj<0.05 & 616 WT-like)": robust.index,
}
scores = pd.DataFrame({k: signature_score(v, g.LFC_5xFAD_vs_WT) for k, v in SIGS.items() if len(v)})
scores.to_csv(os.path.join(OUT, "signature_scores_per_sample.csv"))

# exact permutation of SPP1 labels within sex among 5xFAD animals
F5 = meta.loc[FIVE]
perms = []
fem, mal = F5.index[F5.sex == "F"], F5.index[F5.sex == "M"]
for cf in itertools.combinations(fem, 3):
    for cm in itertools.combinations(mal, 3):
        perms.append(set(cf) | set(cm))
obs_set = set(SPP)


def spp1_effect(score, spp_set):
    s = score.loc[FIVE]
    in_s = s.index.isin(list(spp_set))
    # sex-balanced difference: mean over sexes of (SPP1 - NTC)
    eff = []
    for sx in ["F", "M"]:
        m = (F5.sex == sx).values
        eff.append(s[m & in_s].mean() - s[m & ~in_s].mean())
    return np.mean(eff)


perm_res = []
for k, sc in scores.items():
    obs = spp1_effect(sc, obs_set)
    null = np.array([spp1_effect(sc, p) for p in perms])
    p2 = (np.abs(null) >= abs(obs) - 1e-12).mean()
    gap = sc.loc[NTC].mean() - sc.loc[WT_OTHER].mean()  # NTC vs WT gap on the (symmetric) signature
    perm_res.append({"signature": k, "n_genes": len(SIGS[k]), "SPP1_minus_NTC": obs, "perm_p_two_sided": p2,
                     "NTC_minus_WT": gap, "fraction_of_5xFAD_gap_shifted_by_SPP1": -obs / gap if gap else np.nan})
perm_res = pd.DataFrame(perm_res)

# concordance across all fibroblast genes: SPP1 effect vs 5xFAD effect (sample-level permutation)
def spp1_lfc_vec(spp_set):
    vv = vst.loc[g.index, FIVE]
    in_s = vv.columns.isin(list(spp_set))
    eff = []
    for sx in ["F", "M"]:
        m = (F5.sex == sx).values
        eff.append(vv.loc[:, m & in_s].mean(axis=1) - vv.loc[:, m & ~in_s].mean(axis=1))
    return (eff[0] + eff[1]) / 2


conc = []
for label, ref_lfc in [("5xFAD effect, unadjusted", g.LFC_5xFAD_vs_WT), ("5xFAD effect, % mito adjusted", g.LFC_5xFAD_vs_WT_mitoadj)]:
    obs = spearmanr(ref_lfc, spp1_lfc_vec(obs_set))[0]
    null = np.array([spearmanr(ref_lfc, spp1_lfc_vec(p))[0] for p in perms])
    conc.append({"reference": label, "spearman_rho(SPP1 effect, 5xFAD effect)": obs,
                 "perm_p_two_sided": (np.abs(null) >= abs(obs) - 1e-12).mean(),
                 "interpretation": "negative = SPP1 moves fibroblasts towards WT"})
conc = pd.DataFrame(conc)
perm_res.to_csv(os.path.join(OUT, "spp1_signature_permutation_tests.csv"), index=False)
conc.to_csv(os.path.join(OUT, "spp1_vs_5xFAD_concordance.csv"), index=False)
print(perm_res.round(3).to_string())
print(conc.round(3).to_string())

# ---------------------------------------------------------------- figures
def dot_groups(ax, vals, highlight616=True):
    for i, grp in enumerate(GROUP_ORDER):
        for sx, mk in SEX_MARKER.items():
            idx = meta.index[(meta.group == grp) & (meta.sex == sx)]
            xs = np.full(len(idx), i) + (-0.12 if sx == "F" else 0.12)
            ax.scatter(xs, vals[idx], c=GROUP_COLORS[grp], marker=mk, s=45, edgecolor=SURFACE, lw=1.2, zorder=3)
        ax.hlines(vals[meta.index[meta.group == grp]].median(), i - 0.3, i + 0.3, color=INK, lw=1.5)
    if highlight616:
        ax.scatter([0.12], [vals[WT_616]], s=150, facecolor="none", edgecolor=INK, lw=1.3, zorder=4)
    ax.set_xticks(range(3), ["WT\nNTC", "5xFAD\nNTC", "5xFAD\nSPP1"], fontsize=8)
    ax.grid(axis="x", visible=False)


# F1: t616 distribution
fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
ax = axes[0]
bins = np.linspace(-1.5, 2.5, 33)
ax.hist(cand.t616.clip(-1.5, 2.5), bins=bins, color="#eb6834", edgecolor=SURFACE, lw=1.5, label="WT 616 (technically 5xFAD-like)")
for a in loo.columns:
    ax.hist(loo[a].clip(-1.5, 2.5), bins=bins, histtype="step", color="#2a78d6", lw=1.2, alpha=0.8,
            label="ordinary WT animals (leave-one-out)" if a == loo.columns[0] else None)
ax.axvline(0, color="#2a78d6", ls="--", lw=1)
ax.axvline(1, color="#eb6834", ls="--", lw=1)
ax.set_xlabel("position between WT (0) and 5xFAD (1)")
ax.set_ylabel(f"genes (of {len(cand)} genotype-DE fibroblast genes)")
ax.set_title("Where does WT animal 616 fall on each 5xFAD-vs-WT gene?")
ax.legend(fontsize=8, loc="upper left")
ax = axes[1]
ax.scatter(g.LFC_5xFAD_vs_WT, g.LFC_5xFAD_vs_WT_mitoadj, s=14, c="#c9c8c3", lw=0)
ax.scatter(cand.LFC_5xFAD_vs_WT, g.loc[cand.index, "LFC_5xFAD_vs_WT_mitoadj"], s=18, c="#eb6834", lw=0, label="padj<0.05 unadjusted")
ax.axhline(0, color=INK2, lw=0.6)
ax.axvline(0, color=INK2, lw=0.6)
lim = np.nanmax(np.abs(g[["LFC_5xFAD_vs_WT", "LFC_5xFAD_vs_WT_mitoadj"]].values)) * 1.05
ax.plot([-lim, lim], [-lim, lim], color=INK2, lw=0.6, ls=":")
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_xlabel("5xFAD vs WT log2FC (unadjusted)")
ax.set_ylabel("5xFAD vs WT log2FC (% mito adjusted)")
ax.set_title("Genotype effects shrink toward 0 once the technical axis is modelled")
ax.legend(fontsize=8)
fig.tight_layout()
save(fig, "G1_genotype_vs_technical.png")

# F2: top genotype genes, per-sample (616 circled)
top = cand.sort_values("padj_5xFAD_vs_WT").head(12).index
if len(top):
    fig, axes = plt.subplots(3, 4, figsize=(15, 10), sharex=True)
    for ax, gid in zip(axes.flat, top):
        dot_groups(ax, vst.loc[gid])
        ax.set_title(f"{name[gid]}  (t616={g.t616[gid]:.2f})", fontsize=9.5)
    axes[0, 0].set_ylabel("fibroblast-normalised VST")
    fig.suptitle("Top 5xFAD-vs-WT fibroblast genes (5xFAD = SPP1 + NTC pooled for selection); circled = WT 616 (processed like the 5xFAD libraries)", fontweight="bold")
    fig.tight_layout()
    save(fig, "G2_top_genotype_genes_616.png")

# F3: signature scores and SPP1
fig, axes = plt.subplots(1, len(scores.columns) + 1, figsize=(6 * (len(scores.columns) + 1), 5))
axes = np.atleast_1d(axes)
for ax, (k, sc) in zip(axes, scores.items()):
    dot_groups(ax, sc)
    pr = perm_res.set_index("signature").loc[k]
    ax.set_ylabel("5xFAD fibroblast signature score\n(higher = more 5xFAD-like)")
    ax.set_title(f"{k}\nn={int(pr.n_genes)} genes; SPP1-NTC={pr.SPP1_minus_NTC:.2f}, permutation p={pr.perm_p_two_sided:.2f}", fontsize=9)
ax = axes[-1]
ax.scatter(g.LFC_5xFAD_vs_WT_mitoadj, g.LFC_SPP1_vs_NTC, s=16, c="#c9c8c3", lw=0)
lab = g.reindex(g.LFC_SPP1_vs_NTC.abs().sort_values(ascending=False).index[:8])
for gid, r in lab.iterrows():
    ax.annotate(r.gene, (r.LFC_5xFAD_vs_WT_mitoadj, r.LFC_SPP1_vs_NTC), xytext=(3, 2), textcoords="offset points", fontsize=7)
c2 = conc.set_index("reference").loc["5xFAD effect, % mito adjusted"]
ax.axhline(0, color=INK2, lw=0.6)
ax.axvline(0, color=INK2, lw=0.6)
ax.set_xlabel("5xFAD vs WT log2FC (% mito adjusted)")
ax.set_ylabel("SPP1 vs NTC log2FC")
ax.set_title(f"Does SPP1 oppose the 5xFAD effect? (all {len(g)} fibroblast genes)\n"
             f"ρ={c2['spearman_rho(SPP1 effect, 5xFAD effect)']:.2f}, permutation p={c2.perm_p_two_sided:.2f} (negative = toward WT)", fontsize=9)
fig.tight_layout()
save(fig, "G3_spp1_on_5xfad_signature.png")
print("Done:", OUT)
