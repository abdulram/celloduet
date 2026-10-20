"""12b - PRIMARY contrast 5xFAD_SPP1 vs 5xFAD_NTC: composition balance, module effects, robust genes.

1. Composition balance: is SPP1 vs NTC confounded by dissection? For every composition signature, Composition PC
   and technical metric: difference (SD units), AUC, exact within-sex permutation p (400 relabellings).
2. Module scores (all gene sets): SPP1 effect in SD units (~ sex), (~ sex + Composition_PC1 + Composition_PC3),
   and exact within-sex permutation p for the unadjusted effect.
3. Robust genes (padj < 0.05 in >= 3 of 5 DESeq2 models, step 12a) + listed candidates: per-animal values and
   per-sex effects.
Outputs: 12_spp1_vs_ntc/, tables/SPP1_*.tsv
"""
import itertools
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
CA = os.path.dirname(HERE)
OUT = os.path.join(CA, "12_spp1_vs_ntc")
TAB = os.path.join(CA, "tables")
exec(open(os.path.join(HERE, "common.py")).read())

meta = pd.read_csv(os.path.join(TAB, "sample_metadata_with_composition.tsv"), sep="\t", index_col=0)
S = pd.read_csv(os.path.join(TAB, "composition_signature_scores_all.tsv"), sep="\t", index_col=0)
vst = pd.read_csv(os.path.join(CA, "01_QC", "vst_blind.tsv.gz"), sep="\t", index_col=0)
mods = pd.read_csv(os.path.join(TAB, "pathway_module_scores.tsv"), sep="\t", index_col=0)
allm = pd.read_csv(os.path.join(TAB, "SPP1_vs_NTC_all_models.tsv.gz"), sep="\t")
F5 = meta[meta.genotype == "5xFAD"].copy()
F5["sp"] = (F5.treatment == "SPP1").astype(int)
F5["fibro"] = S.loc[F5.index, "fibroblast_stromal"]  # fibroblast abundance (lower in SPP1: see composition balance)
fem, mal = F5.index[F5.sex == "F"], F5.index[F5.sex == "M"]
perms = [set(a) | set(b) for a in itertools.combinations(fem, 3) for b in itertools.combinations(mal, 3)]
obs_set = set(F5.index[F5.sp == 1])


def sexbal(y, spset):
    return np.mean([y[[l for l in grp if l in spset]].mean() - y[[l for l in grp if l not in spset]].mean() for grp in (fem, mal)])


def perm_p(y):
    o = sexbal(y, obs_set)
    null = np.array([sexbal(y, p) for p in perms])
    return o, (np.abs(null) >= abs(o) - 1e-12).mean()


# ---------------------------------------------------------------- 1. composition balance
F = pd.concat([meta[["Composition_PC1", "Composition_PC2", "Composition_PC3", "pct_mito", "pct_multi", "pct_hemoglobin"]], S], axis=1).loc[F5.index]
rows = []
for c in F:
    y = (F[c] - F[c].mean()) / F[c].std()
    o, p = perm_p(y)
    auc = mannwhitneyu(F.loc[F5.sp == 1, c], F.loc[F5.sp == 0, c]).statistic / 36
    rows.append({"measure": c, "SPP1_minus_NTC_SD": o, "perm_p": p, "AUC_SPP1_gt_NTC": auc})
bal = pd.DataFrame(rows)
bal.to_csv(os.path.join(TAB, "SPP1_vs_NTC_composition_balance.tsv"), sep="\t", index=False)
print(bal.round(3).to_string())

# ---------------------------------------------------------------- 2. module scores
rows = []
for s in mods.columns:
    y = (mods.loc[F5.index, s] - mods.loc[F5.index, s].mean()) / mods.loc[F5.index, s].std()
    o, p = perm_p(y)
    d = F5.assign(y=y)
    f0 = smf.ols("y ~ C(sex) + sp", d).fit()
    f1 = smf.ols("y ~ C(sex) + Composition_PC1 + Composition_PC3 + sp", d).fit()
    f2 = smf.ols("y ~ C(sex) + pct_mito + sp", d).fit()
    f3 = smf.ols("y ~ C(sex) + fibro + sp", d).fit()
    rows.append({"SPP1_effect_SD_fibroblast_adj": f3.params["sp"], "p_fibroblast_adj": f3.pvalues["sp"], "module": s, "SPP1_effect_SD": f0.params["sp"], "p_ols": f0.pvalues["sp"], "perm_p_within_sex": p,
                 "SPP1_effect_SD_compPC1PC3": f1.params["sp"], "p_compPC1PC3": f1.pvalues["sp"],
                 "SPP1_effect_SD_mito": f2.params["sp"], "p_mito": f2.pvalues["sp"]})
mres = pd.DataFrame(rows).sort_values("perm_p_within_sex")
cam = pd.read_csv(os.path.join(TAB, "pathways_camera_SPP1_vs_NTC.tsv"), sep="\t").set_index("set")
mres["camera_p"] = mres.module.map(cam.PValue)
mres["camera_FDR"] = mres.module.map(cam.FDR)
mres.to_csv(os.path.join(TAB, "SPP1_vs_NTC_module_effects.tsv"), sep="\t", index=False)
print(mres.round(3).head(15).to_string())

# ---------------------------------------------------------------- 3. robust + candidate genes
robust = allm[allm.n_models_padj05 >= 3].gene_name.tolist()
cands = ["Cldn11", "Prg4", "Rgs4", "Fn1", "Fbln1", "Thbd", "Matn4", "Ifit3b", "Loxl3", "Col12a1", "Wwtr1", "Wfdc17", "Spp1"]
genes = list(dict.fromkeys(robust + cands))
g2i = allm.set_index("gene_name").gene_id
rows = []
for g in genes:
    if g not in g2i.index or g2i[g] not in vst.index:
        continue
    x = vst.loc[g2i[g], F5.index]
    o, p = perm_p((x - x.mean()) / x.std())
    r = allm.set_index("gene_name").loc[g]
    dz = F5.assign(y=(x - x.mean()) / x.std())
    ff = smf.ols("y ~ C(sex) + fibro + sp", dz).fit()
    fc = smf.ols("y ~ C(sex) + Composition_PC1 + Composition_PC3 + sp", dz).fit()
    rows.append({"SPP1_effect_SD": smf.ols("y ~ C(sex) + sp", dz).fit().params["sp"],
                 "SPP1_effect_SD_fibroblast_adj": ff.params["sp"], "p_fibroblast_adj": ff.pvalues["sp"],
                 "SPP1_effect_SD_compPC1PC3": fc.params["sp"], "p_compPC1PC3": fc.pvalues["sp"], "gene": g, "robust_all_models": g in robust, "baseMean": r.baseMean, "LFC_baseline": r.LFC_baseline,
                 "padj_baseline": r.padj_baseline, "padj_comp1": r.padj_comp1, "padj_comp13": r.padj_comp13, "padj_tech": r.padj_tech,
                 "padj_5xFADonly": r.padj_5xFADonly, "n_models_padj05": r.n_models_padj05,
                 "VST_diff_females": x[[l for l in fem if l in obs_set]].median() - x[[l for l in fem if l not in obs_set]].median(),
                 "VST_diff_males": x[[l for l in mal if l in obs_set]].median() - x[[l for l in mal if l not in obs_set]].median(),
                 "perm_p_within_sex": p})
gres = pd.DataFrame(rows)
gres.to_csv(os.path.join(TAB, "SPP1_vs_NTC_key_genes.tsv"), sep="\t", index=False)
print(gres.round(3).to_string())

# ---------------------------------------------------------------- figures
groups = ["5xFAD_NTC", "5xFAD_SPP1"]
show = ["Composition_PC1", "Composition_PC2", "Composition_PC3", "pct_mito", "fibroblast_stromal", "blood_endothelial",
        "bone_osteoblast", "osteoclast", "olfactory_nasal", "macrophage_BAM", "erythroid_blood", "CNS_neural"]
fig, axes = plt.subplots(3, 4, figsize=(15, 10))
for ax, c in zip(axes.flat, show):
    strip(ax, F[c], meta.loc[F5.index], groups=groups, label_animals=True)
    b = bal.set_index("measure").loc[c]
    ax.set_title(f"{c}\nSPP1-NTC {b.SPP1_minus_NTC_SD:+.2f} SD, perm p={b.perm_p:.2f}", fontsize=9,
                 color="#d03b3b" if b.perm_p < 0.05 else INK)
fig.suptitle("Is SPP1 vs NTC balanced on dissection composition and technical metrics? (5xFAD animals only)", fontweight="bold")
fig.tight_layout()
save(fig, os.path.join(OUT, "SPP1_A_composition_balance.png"))

gplot = gres.gene.tolist()
n = len(gplot)
fig, axes = plt.subplots(int(np.ceil(n / 5)), 5, figsize=(17, 3.4 * int(np.ceil(n / 5))))
for ax, g in zip(axes.flat, gplot):
    x = vst.loc[g2i[g]]
    strip(ax, x, meta.loc[F5.index], groups=groups)
    r = gres.set_index("gene").loc[g]
    ax.set_title(f"{g}{' *' if r.robust_all_models else ''}\nLFC {r.LFC_baseline:+.2f}; models padj<.05: {int(r.n_models_padj05)}/5; perm p={r.perm_p_within_sex:.3f}", fontsize=8.5)
for ax in list(axes.flat)[n:]:
    ax.axis("off")
fig.suptitle("Key SPP1-vs-NTC genes per animal (VST; circle F, triangle M; * = padj<0.05 in >=3 of 5 models)", fontweight="bold")
fig.tight_layout()
save(fig, os.path.join(OUT, "SPP1_B_key_genes_per_animal.png"))

mm = mres.set_index("module").sort_values("SPP1_effect_SD")
fig, ax = plt.subplots(figsize=(10, 11))
y = np.arange(len(mm))
ax.scatter(mm.SPP1_effect_SD, y, c="#0b0b0b", s=28, label="~ sex", zorder=3)
ax.scatter(mm.SPP1_effect_SD_compPC1PC3, y + 0.25, c="#2a78d6", marker="s", s=22, label="+ Composition PC1 + PC3", zorder=3)
ax.scatter(mm.SPP1_effect_SD_mito, y - 0.25, c="#eb6834", marker="^", s=22, label="+ % mito", zorder=3)
ax.scatter(mm.SPP1_effect_SD_fibroblast_adj, y, c="#1baf7a", marker="D", s=18, label="+ fibroblast abundance", zorder=2)
for i, (s, r) in enumerate(mm.iterrows()):
    if r.perm_p_within_sex < 0.1:
        ax.text(2.6, i, f"p={r.perm_p_within_sex:.2f}", va="center", fontsize=7, color="#d03b3b")
ax.axvline(0, color=INK2, lw=0.8)
ax.set_yticks(y, [s.replace("MSIGDB_HALLMARK_", "H: ").replace("CURATED_", "").replace("TECHNICAL_", "TECH: ") for s in mm.index], fontsize=7.5)
ax.set_xlim(-3, 3.2)
ax.set_xlabel("SPP1 vs NTC effect on module score (SD units, 5xFAD animals)")
ax.set_title("Pathway/module effects of SPP1 (red = within-sex permutation p < 0.1)")
ax.legend(fontsize=8, loc="upper left")
ax.grid(axis="y", visible=False)
save(fig, os.path.join(OUT, "SPP1_C_module_effects.png"))
