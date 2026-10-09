"""11b - TGF-beta target-gene panel: cell origin, coherence, panel score, genotype and SPP1 effects.

Panel (user-supplied confidence -> weight): High = 1 (Serpine1, Smad7, Tgfbi, Pmepa1), Medium = 0.5 (Skil, Ccn2/Ctgf,
Fn1), Low = 0.25 (Smurf2). Context (pathway input, not part of the score): Tgfb1-3, Tgfbr1-2, Ltbp1.

1. Cell origin: dura-centric single-cell reference (CP10K per class) + bulk correlation with composition signatures.
2. Coherence: do the readouts co-vary across animals like a pathway (pairwise r; item-rest r)?
3. Panel scores (weighted mean z of blind VST): all-8 weighted, High-only (4), and a "residual" score in which each
   gene is first regressed on Composition_PC1/PC3 + fibroblast + endothelial abundance (composition-independent part).
4. Effects (SD units of the 12 NTC animals for genotype; 12 5xFAD animals for SPP1):
   genotype (WT_NTC vs 5xFAD_NTC): ~sex ; +Composition_PC1 ; +PC1+PC3 ; +pct_mito ; + fibroblast abundance
   WT 616 position (technical-axis test, see step 06b)
   SPP1 (5xFAD_SPP1 vs 5xFAD_NTC): ~sex ; +PC1 ; +PC1+PC3 ; exact within-sex permutation (400)
Outputs: 11_tgfb_panel/, tables/tgfb_*.tsv
"""
import itertools
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

HERE = os.path.dirname(os.path.abspath(__file__))
CA = os.path.dirname(HERE)
ROOT = os.path.dirname(CA)
OUT = os.path.join(CA, "11_tgfb_panel")
TAB = os.path.join(CA, "tables")
exec(open(os.path.join(HERE, "common.py")).read())

PANEL = {"Serpine1": ("High", 1.0), "Smad7": ("High", 1.0), "Tgfbi": ("High", 1.0), "Pmepa1": ("High", 1.0),
         "Skil": ("Medium", 0.5), "Ccn2": ("Medium", 0.5), "Fn1": ("Medium", 0.5), "Smurf2": ("Low", 0.25)}
CONTEXT = ["Tgfb1", "Tgfb2", "Tgfb3", "Tgfbr1", "Tgfbr2", "Ltbp1"]

meta = pd.read_csv(os.path.join(TAB, "sample_metadata_with_composition.tsv"), sep="\t", index_col=0)
vst = pd.read_csv(os.path.join(CA, "01_QC", "vst_blind.tsv.gz"), sep="\t", index_col=0)
S = pd.read_csv(os.path.join(TAB, "composition_signature_scores_all.tsv"), sep="\t", index_col=0)
raw = pd.read_csv(os.path.join(ROOT, "data", "VT3YXB-expression-matrix.tsv"), sep="\t").set_index("gene_id")
libs = list(meta.index)
cpm = raw[[f"{l}_cpm" for l in libs]]
cpm.columns = libs
n2i = pd.Series(raw.index, index=raw.gene_name.values)
n2i = n2i[~n2i.index.duplicated()]
ALL = list(PANEL) + CONTEXT
ids = {g: n2i[g] for g in ALL}
de = pd.read_csv(os.path.join(TAB, "tgfb_panel_DE_all_models.tsv"), sep="\t")

# ---------------------------------------------------------------- 1. cell origin
ref = pd.read_csv(os.path.join(CA, "reference", "dura_reference_pseudobulk_cp10k.csv.gz"), index_col=0)
origin = ref.reindex([ids[g] for g in ALL])
origin.index = ALL
origin.to_csv(os.path.join(TAB, "tgfb_panel_reference_expression_cp10k.tsv"), sep="\t")
share = origin.div(origin.sum(axis=1), axis=0)
top_cls = share.apply(lambda r: ", ".join(f"{c} {v*100:.0f}%" for c, v in r.sort_values(ascending=False).head(3).items()), axis=1)
sigs = ["fibroblast_stromal", "blood_endothelial", "macrophage_BAM", "bone_osteoblast", "osteoclast", "olfactory_nasal", "CNS_neural",
        "skeletal_muscle", "erythroid_blood"]
bulk_r = pd.DataFrame({s: [np.corrcoef(vst.loc[ids[g]], S[s])[0, 1] for g in ALL] for s in sigs}, index=ALL)
bulk_r.to_csv(os.path.join(TAB, "tgfb_panel_bulk_r_with_composition.tsv"), sep="\t")

# ---------------------------------------------------------------- 2. coherence
zv = vst.sub(vst.mean(axis=1), axis=0).div(vst.std(axis=1), axis=0)
Zp = zv.loc[[ids[g] for g in PANEL]]
Zp.index = list(PANEL)
Cm = Zp.T.corr()
Cm.to_csv(os.path.join(TAB, "tgfb_panel_gene_correlations.tsv"), sep="\t")
w = pd.Series({g: PANEL[g][1] for g in PANEL})
item_rest = {g: np.corrcoef(Zp.loc[g], Zp.drop(g).mul(w.drop(g), axis=0).sum() / w.drop(g).sum())[0, 1] for g in PANEL}
mean_pair = Cm.values[np.triu_indices(len(Cm), 1)].mean()
# null: mean pairwise r of 8 random expressed genes with matched expression (1000 draws)
rng = np.random.default_rng(7)
expr = vst.mean(axis=1)
pool = vst.index[(expr > expr.loc[[ids[g] for g in PANEL]].min()) & (expr < expr.loc[[ids[g] for g in PANEL]].max())]
null = []
for _ in range(1000):
    s = rng.choice(pool, len(PANEL), replace=False)
    cc = zv.loc[s].T.corr().values
    null.append(cc[np.triu_indices(len(PANEL), 1)].mean())
coh_p = (np.array(null) >= mean_pair).mean()
print(f"panel mean pairwise r = {mean_pair:.2f}; random-gene null mean {np.mean(null):.2f}, p = {coh_p:.3f}")

# ---------------------------------------------------------------- 3. scores
hi = [g for g in PANEL if PANEL[g][0] == "High"]
score_w = Zp.mul(w, axis=0).sum() / w.sum()
score_hi = Zp.loc[hi].mean()
# composition-residualised genes
covs = meta[["Composition_PC1", "Composition_PC3"]].join(S[["fibroblast_stromal", "blood_endothelial"]])
X = np.column_stack([np.ones(len(libs))] + [((covs[c] - covs[c].mean()) / covs[c].std()).values for c in covs])
B, *_ = np.linalg.lstsq(X, Zp.T.values, rcond=None)
R = pd.DataFrame((Zp.T.values - X @ B).T, index=Zp.index, columns=libs)
R = R.div(R.std(axis=1), axis=0)
score_res = R.mul(w, axis=0).sum() / w.sum()
scores = pd.DataFrame({"TGFb_weighted_all8": score_w, "TGFb_high_conf_4": score_hi, "TGFb_weighted_composition_residual": score_res})
scores.to_csv(os.path.join(TAB, "tgfb_panel_scores.tsv"), sep="\t")

# ---------------------------------------------------------------- 4. effects
P = meta[meta.group.isin(["WT_NTC", "5xFAD_NTC"])].copy()
P["g5"] = (P.genotype == "5xFAD").astype(int)
P = P.join(S[["fibroblast_stromal", "blood_endothelial"]])
F5 = meta[meta.genotype == "5xFAD"].copy()
F5["sp"] = (F5.treatment == "SPP1").astype(int)
a616 = meta.index[meta.animal.astype(str) == "616"][0]
WTo = P.index[(P.genotype == "WT") & (P.index != a616)]
FX = P.index[P.genotype == "5xFAD"]
fem, mal = F5.index[F5.sex == "F"], F5.index[F5.sex == "M"]
perms = [set(a) | set(b) for a in itertools.combinations(fem, 3) for b in itertools.combinations(mal, 3)]


def spp1_diff(sc, spset):
    d = [sc[[l for l in (fem if sx == "F" else mal) if l in spset]].mean() - sc[[l for l in (fem if sx == "F" else mal) if l not in spset]].mean() for sx in "FM"]
    return np.mean(d)


rows = []
targets = {**{k: v for k, v in scores.items()}, **{g: vst.loc[ids[g]] for g in ALL}}
for name, sc in targets.items():
    y = (sc - sc[P.index].mean()) / sc[P.index].std()
    d = P.assign(y=y[P.index])
    r = {"readout": name, "type": "score" if name.startswith("TGFb") else ("panel gene" if name in PANEL else "context gene"),
         "confidence": PANEL.get(name, ("", 0))[0]}
    for lab, f in [("geno_base", "y ~ C(sex) + g5"), ("geno_compPC1", "y ~ C(sex) + Composition_PC1 + g5"),
                   ("geno_compPC1PC3", "y ~ C(sex) + Composition_PC1 + Composition_PC3 + g5"), ("geno_mito", "y ~ C(sex) + pct_mito + g5"),
                   ("geno_fibroblast", "y ~ C(sex) + fibroblast_stromal + g5")]:
        m = smf.ols(f, d).fit()
        r[f"{lab}_SD"], r[f"{lab}_p"] = m.params["g5"], m.pvalues["g5"]
    gap = sc[FX].mean() - sc[WTo].mean()
    r["position_616"] = (sc[a616] - sc[WTo].mean()) / gap if abs(gap) > 1e-9 else np.nan
    y5 = (sc - sc[F5.index].mean()) / sc[F5.index].std()
    d5 = F5.assign(y=y5[F5.index])
    for lab, f in [("spp1_base", "y ~ C(sex) + sp"), ("spp1_compPC1", "y ~ C(sex) + Composition_PC1 + sp"),
                   ("spp1_compPC1PC3", "y ~ C(sex) + Composition_PC1 + Composition_PC3 + sp")]:
        m = smf.ols(f, d5).fit()
        r[f"{lab}_SD"], r[f"{lab}_p"] = m.params["sp"], m.pvalues["sp"]
    obs = spp1_diff(y5, set(F5.index[F5.sp == 1]))
    nullp = np.array([spp1_diff(y5, p) for p in perms])
    r["spp1_perm_p"] = (np.abs(nullp) >= abs(obs) - 1e-12).mean()
    rows.append(r)
eff = pd.DataFrame(rows)
eff.to_csv(os.path.join(TAB, "tgfb_panel_effects.tsv"), sep="\t", index=False)

summ = pd.DataFrame({"confidence": [PANEL.get(g, ("context",))[0] for g in ALL], "median_CPM_WT": [cpm.loc[ids[g], P.index[P.genotype == 'WT']].median() for g in ALL],
                     "median_CPM_5xFAD_NTC": [cpm.loc[ids[g], FX].median() for g in ALL],
                     "median_CPM_5xFAD_SPP1": [cpm.loc[ids[g], F5.index[F5.sp == 1]].median() for g in ALL],
                     "top_reference_cell_types": top_cls.values,
                     "item_rest_r": [item_rest.get(g, np.nan) for g in ALL]}, index=ALL)
summ.to_csv(os.path.join(TAB, "tgfb_panel_gene_summary.tsv"), sep="\t")
pd.set_option("display.width", 250)
print(summ.round(2).to_string())
print(bulk_r.round(2).to_string())
print(Cm.round(2).to_string())
print(eff[["readout", "geno_base_SD", "geno_base_p", "geno_compPC1_SD", "geno_compPC1PC3_SD", "geno_mito_SD", "geno_mito_p", "geno_fibroblast_SD",
           "position_616", "spp1_base_SD", "spp1_base_p", "spp1_compPC1PC3_SD", "spp1_perm_p"]].round(3).to_string())

# ---------------------------------------------------------------- figures
# A: per-gene strip plots (VST), 616 circled
fig, axes = plt.subplots(2, 4, figsize=(16, 8), sharex=True)
for ax, g in zip(axes.flat, PANEL):
    strip(ax, vst.loc[ids[g]], meta, label_animals=False)
    ax.scatter([0.12], [vst.loc[ids[g], a616]], s=160, facecolor="none", edgecolor=INK, lw=1.2, zorder=4)
    e = eff.set_index("readout").loc[g]
    ax.set_title(f"{g} ({PANEL[g][0]})\n5xFAD-WT p={e.geno_base_p:.3f} | +%mito p={e.geno_mito_p:.2f} | SPP1 perm p={e.spp1_perm_p:.2f}", fontsize=8.5)
axes[0, 0].set_ylabel("VST expression")
axes[1, 0].set_ylabel("VST expression")
fig.suptitle("TGF-β target panel per animal (circle F, triangle M; ringed = WT 616, the WT library with the 5xFAD-like technical profile)", fontweight="bold")
fig.tight_layout()
save(fig, os.path.join(OUT, "TGFb_A_panel_genes_per_animal.png"))

# B: effect heatmap across models
cols = ["geno_base_SD", "geno_compPC1_SD", "geno_compPC1PC3_SD", "geno_mito_SD", "geno_fibroblast_SD", "spp1_base_SD", "spp1_compPC1PC3_SD"]
pcols = [c.replace("_SD", "_p") for c in cols]
E = eff.set_index("readout").loc[list(scores.columns) + list(PANEL) + CONTEXT]
fig, ax = plt.subplots(figsize=(11, 9))
im = ax.imshow(E[cols].values, cmap=DIV, vmin=-2.5, vmax=2.5, aspect="auto")
ax.set_xticks(range(len(cols)), ["5xFAD vs WT\nbaseline", "+CompPC1", "+CompPC1\n+PC3", "+%mito", "+fibroblast\nabund.", "SPP1 vs NTC\nbaseline", "SPP1\n+PC1+PC3"], fontsize=8)
ax.set_yticks(range(len(E)), [f"{i}  [{E.loc[i, 'confidence'] or E.loc[i, 'type']}]" for i in E.index], fontsize=8)
for i in range(len(E)):
    for j, pc in enumerate(pcols):
        v, p = E[cols[j]].iloc[i], E[pc].iloc[i]
        ax.text(j, i, f"{v:.1f}{'*' if p < 0.05 else ''}", ha="center", va="center", fontsize=7, color="white" if abs(v) > 1.6 else INK)
ax.axhline(len(scores.columns) - 0.5, color=INK, lw=1)
ax.axhline(len(scores.columns) + len(PANEL) - 0.5, color=INK, lw=1)
ax.axvline(4.5, color=INK, lw=1)
ax.grid(False)
fig.colorbar(im, ax=ax, shrink=0.5, label="effect (SD units); * p<0.05 (unadjusted)")
ax.set_title("TGF-β panel: genotype and SPP1 effects under each model")
save(fig, os.path.join(OUT, "TGFb_B_effects_by_model.png"))

# C: scores per animal
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
for ax, c in zip(axes, scores.columns):
    strip(ax, scores[c], meta, label_animals=True)
    ax.scatter([0.12], [scores.loc[a616, c]], s=160, facecolor="none", edgecolor=INK, lw=1.2, zorder=4)
    e = eff.set_index("readout").loc[c]
    ax.set_title(f"{c}\n5xFAD-WT {e.geno_base_SD:.2f} SD (p={e.geno_base_p:.3f}); +%mito p={e.geno_mito_p:.2f}\nSPP1-NTC {e.spp1_base_SD:.2f} SD, perm p={e.spp1_perm_p:.2f}", fontsize=8.5)
save(fig, os.path.join(OUT, "TGFb_C_panel_scores.png"))

# D: coherence + cell origin
fig, axes = plt.subplots(1, 2, figsize=(16, 6.5), gridspec_kw={"width_ratios": [1, 1.6]})
ax = axes[0]
im = ax.imshow(Cm.values, cmap=DIV, vmin=-1, vmax=1)
ax.set_xticks(range(len(Cm)), Cm.columns, rotation=45, ha="right", fontsize=8)
ax.set_yticks(range(len(Cm)), Cm.index, fontsize=8)
for i in range(len(Cm)):
    for j in range(len(Cm)):
        ax.text(j, i, f"{Cm.values[i, j]:.2f}", ha="center", va="center", fontsize=7, color="white" if abs(Cm.values[i, j]) > 0.6 else INK)
ax.grid(False)
ax.set_title(f"Do the readouts co-vary? mean r = {mean_pair:.2f}\n(random expression-matched 8-gene sets: {np.mean(null):.2f}; p = {coh_p:.2f})", fontsize=9.5)
ax = axes[1]
L = np.log10(origin.loc[ALL] + 0.01)
im = ax.imshow(L.values, cmap=SEQ, aspect="auto", vmin=-2)
ax.set_yticks(range(len(ALL)), ALL, fontsize=8)
ax.set_xticks(range(L.shape[1]), L.columns, rotation=90, fontsize=7.5)
ax.axhline(len(PANEL) - 0.5, color=INK, lw=1)
ax.grid(False)
fig.colorbar(im, ax=ax, shrink=0.6, label="log10 CP10K (dura-centric single-cell reference)")
ax.set_title("Which cells express each gene?")
fig.tight_layout()
save(fig, os.path.join(OUT, "TGFb_D_coherence_and_cell_origin.png"))
