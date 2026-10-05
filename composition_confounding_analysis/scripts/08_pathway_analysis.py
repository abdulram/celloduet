"""08 - Pathway / module robustness matrix (central output).

Module score = mean z (blind VST) of the set's genes. Animals: WT_NTC + 5xFAD_NTC (n=12) unless stated.
Genotype effect in SD units (score standardised across the 12 animals), from OLS:
  baseline        score ~ sex + genotype
  comp_adjusted   score ~ sex + Composition_PC1 + genotype
  comp13_adjusted score ~ sex + Composition_PC1 + Composition_PC3 + genotype
  tech_adjusted   score ~ sex + pct_mito + genotype
  matched         mean(5xFAD) - mean(WT) in the MatchIt cohort from step 06 (n = 2 vs 2; direction only)
  616 position    WT animal 616 placed between other WT (0) and 5xFAD_NTC (1) (technical-axis test, step 06b)
camera FDRs from steps 04/05 are attached. Lineage association and residual genotype-given-abundance from step 07.

Confidence tiers (pre-specified rules; pathology unavailable so 'associated with pathology' cannot be met):
  no baseline association : baseline p >= 0.05
  Tier 4 uninterpretable  : baseline p < 0.05 but effect follows the technical axis
                            (tech-adjusted retention < 0.5 OR 616 position > 0.5)
  Tier 3 composition-assoc: not Tier 4, comp-adjusted retention < 0.5 (or comp13 retention < 0.5)
  Tier 2 moderately robust: retention >= 0.5 in comp and tech models, comp-adjusted p < 0.05, matched direction
                            agrees, but matched/common-support evidence limited
  Tier 1 most robust      : Tier 2 + comp13 p < 0.05 + tech-adjusted p < 0.05 (+ pathology association - unavailable)
"""
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

HERE = os.path.dirname(os.path.abspath(__file__))
CA = os.path.dirname(HERE)
OUT = os.path.join(CA, "08_pathways")
TAB = os.path.join(CA, "tables")
exec(open(os.path.join(HERE, "common.py")).read())

meta = pd.read_csv(os.path.join(TAB, "sample_metadata_with_composition.tsv"), sep="\t", index_col=0)
vst = pd.read_csv(os.path.join(CA, "01_QC", "vst_blind.tsv.gz"), sep="\t", index_col=0)
gs = pd.read_csv(os.path.join(TAB, "gene_sets.tsv"), sep="\t")
coh = pd.read_csv(os.path.join(TAB, "common_support_cohorts_animals.tsv"), sep="\t")
state = pd.read_csv(os.path.join(TAB, "cell_state_vs_abundance.tsv"), sep="\t")
S = pd.read_csv(os.path.join(TAB, "composition_signature_scores_all.tsv"), sep="\t", index_col=0)
cam = {m: pd.read_csv(os.path.join(TAB, f"pathways_camera_{m}.tsv"), sep="\t").set_index("set") for m in ["baseline", "comp1", "comp13", "tech"]}

LINEAGE = {"COMPLEMENT": "macrophage_BAM", "IFN": "macrophage_BAM", "INTERFERON": "macrophage_BAM", "MHCII": "macrophage_BAM",
           "BAM": "macrophage_BAM", "PHAGOCYTOSIS": "macrophage_BAM", "LIPID": "macrophage_BAM", "DAM": "macrophage_BAM",
           "INNATE": "macrophage_BAM", "INFLAMMATORY_RESPONSE": "macrophage_BAM", "TNFA": "macrophage_BAM", "ALLOGRAFT": "macrophage_BAM",
           "MONOCYTE": "monocyte", "T_CELL": "T_cell", "B_PLASMA": "B_cell", "ECM": "fibroblast_stromal", "FIBROBLAST": "fibroblast_stromal",
           "EPITHELIAL_MESENCHYMAL": "fibroblast_stromal", "TGF_BETA": "fibroblast_stromal", "VASCULAR": "blood_endothelial",
           "ANGIOGENESIS": "blood_endothelial", "COAGULATION": "blood_endothelial", "LYMPHATIC": "lymphatic_endothelial",
           "TRAFFICKING": "blood_endothelial", "CHEMOKINES": "macrophage_BAM", "IL6": "macrophage_BAM",
           "MITO": "technical", "RIBOSOME": "technical", "OXIDATIVE": "technical", "HYPOXIA": "-", "APOPTOSIS": "-",
           "REACTIVE_OXYGEN": "-", "CHOLESTEROL": "macrophage_BAM"}


def lineage_of(s):
    for k, v in LINEAGE.items():
        if k in s:
            return v
    return "-"


P = meta[meta.group.isin(["WT_NTC", "5xFAD_NTC"])].copy()
P["g5"] = (P.genotype == "5xFAD").astype(int)
zv = vst.sub(vst.mean(axis=1), axis=0).div(vst.std(axis=1).replace(0, np.nan), axis=0)
matched = coh[coh.cohort == "matched"].animal.astype(str)
idx_by_animal = pd.Series(meta.index, index=meta.animal.astype(str))
M_WT = idx_by_animal[matched[coh[coh.cohort == "matched"].group.values == "WT_NTC"]].values
M_FX = idx_by_animal[matched[coh[coh.cohort == "matched"].group.values == "5xFAD_NTC"]].values
a616 = idx_by_animal["616"]
WT_OTHER = P.index[(P.genotype == "WT") & (P.index != a616)]
FX = P.index[P.genotype == "5xFAD"]

rows, scores = [], {}
for s, g in gs.groupby("set"):
    ids = [i for i in g.gene_id if i in vst.index]
    if len(ids) < 5:
        continue
    sc = zv.loc[ids].mean()
    scores[s] = sc
    d = P.assign(y=(sc[P.index] - sc[P.index].mean()) / sc[P.index].std())
    fits = {"baseline": smf.ols("y ~ C(sex) + g5", d).fit(),
            "comp_adjusted": smf.ols("y ~ C(sex) + Composition_PC1 + g5", d).fit(),
            "comp13_adjusted": smf.ols("y ~ C(sex) + Composition_PC1 + Composition_PC3 + g5", d).fit(),
            "tech_adjusted": smf.ols("y ~ C(sex) + pct_mito + g5", d).fit()}
    r = {"pathway": s, "n_genes": len(ids), "associated_lineage": lineage_of(s)}
    for k, f in fits.items():
        r[f"{k}_effect_SD"], r[f"{k}_p"] = f.params["g5"], f.pvalues["g5"]
    b = r["baseline_effect_SD"]
    for k in ["comp_adjusted", "comp13_adjusted", "tech_adjusted"]:
        r[f"{k}_retention"] = r[f"{k}_effect_SD"] / b if b else np.nan
    ys = (sc - sc[P.index].mean()) / sc[P.index].std()
    r["matched_effect_SD"] = ys[M_FX].mean() - ys[M_WT].mean()
    r["matched_n"] = f"{len(M_WT)} vs {len(M_FX)}"
    gap = sc[FX].mean() - sc[WT_OTHER].mean()
    r["position_616"] = (sc[a616] - sc[WT_OTHER].mean()) / gap if gap else np.nan
    lin = r["associated_lineage"]
    r["r_with_lineage_abundance"] = np.corrcoef(sc[P.index], S.loc[P.index, lin])[0, 1] if lin in S.columns else np.nan
    st = state[state.state_module == s]
    r["residual_genotype_given_abundance_p"] = st.p_given_abundance.iloc[0] if len(st) and "p_given_abundance" in st else np.nan
    for m, c in cam.items():
        r[f"camera_FDR_{m}"] = c.FDR.get(s, np.nan)
    r["pathology_association"] = "not available"
    rows.append(r)
R = pd.DataFrame(rows)


def tier(r):
    if r.associated_lineage == "technical":
        return "technical sentinel (negative control)"
    if r.baseline_p >= 0.05:
        return "no baseline genotype association"
    if (r.tech_adjusted_retention < 0.5) or (r.position_616 > 0.5):
        return "Tier 4 - uninterpretable (follows technical/processing axis)"
    if (r.comp_adjusted_retention < 0.5) or (r.comp13_adjusted_retention < 0.5):
        return "Tier 3 - composition-associated"
    md = np.sign(r.matched_effect_SD) == np.sign(r.baseline_effect_SD)
    if r.comp_adjusted_p < 0.05 and md:
        if r.comp13_adjusted_p < 0.05 and r.tech_adjusted_p < 0.05:
            return "Tier 1 - most robust (pathology unavailable)"
        return "Tier 2 - moderately robust"
    return "Tier 3 - composition-associated"


R["confidence_tier"] = R.apply(tier, axis=1)
R = R.sort_values(["confidence_tier", "baseline_p"])
R.to_csv(os.path.join(TAB, "pathway_robustness_matrix.tsv"), sep="\t", index=False)
pd.DataFrame(scores).to_csv(os.path.join(TAB, "pathway_module_scores.tsv"), sep="\t")
pd.set_option("display.width", 250)
print(R[["pathway", "baseline_effect_SD", "baseline_p", "comp_adjusted_retention", "comp13_adjusted_retention", "tech_adjusted_retention",
         "matched_effect_SD", "position_616", "confidence_tier"]].round(2).to_string())

# Figure 8: effects across models
order = R.sort_values("baseline_effect_SD").pathway
fig, ax = plt.subplots(figsize=(11, 12))
models = [("baseline", "#0b0b0b", "o"), ("comp_adjusted", "#2a78d6", "s"), ("comp13_adjusted", "#1baf7a", "D"), ("tech_adjusted", "#eb6834", "^")]
y = np.arange(len(order))
Ri = R.set_index("pathway").loc[order]
for k, (m, c, mk) in enumerate(models):
    ax.scatter(Ri[f"{m}_effect_SD"], y + (k - 1.5) * 0.17, c=c, marker=mk, s=28, label=m, zorder=3)
ax.scatter(Ri.matched_effect_SD, y + 0.42, c="#9a9994", marker="|", s=60, label="matched (2 vs 2)", zorder=3)
for i, p in enumerate(order):
    if Ri.loc[p, "baseline_p"] < 0.05:
        ax.text(ax.get_xlim()[0] if False else -3.4, i, "*", fontsize=11, va="center", color="#d03b3b")
ax.axvline(0, color=INK2, lw=0.8)
ax.set_yticks(y, [p.replace("MSIGDB_HALLMARK_", "H: ").replace("CURATED_", "").replace("TECHNICAL_", "TECH: ") for p in order], fontsize=7.5)
ax.set_xlim(-3.5, 3.5)
ax.set_xlabel("5xFAD vs WT effect on module score (SD units)")
ax.set_title("Pathway/module genotype effects before and after adjustment (* baseline p<0.05)")
ax.legend(fontsize=8, loc="lower right")
ax.grid(axis="y", visible=False)
save(fig, os.path.join(OUT, "Fig8_pathway_effects_by_model.png"))

# Figure 12: robustness summary heatmap
cols = ["baseline_effect_SD", "comp_adjusted_effect_SD", "comp13_adjusted_effect_SD", "tech_adjusted_effect_SD", "matched_effect_SD"]
H = Ri[cols]
fig, axes = plt.subplots(1, 2, figsize=(14, 12), gridspec_kw={"width_ratios": [1.2, 1]})
ax = axes[0]
im = ax.imshow(H.values, cmap=DIV, vmin=-2.5, vmax=2.5, aspect="auto")
ax.set_yticks(y, [p.replace("MSIGDB_HALLMARK_", "H: ").replace("CURATED_", "").replace("TECHNICAL_", "TECH: ") for p in order], fontsize=7.5)
ax.set_xticks(range(len(cols)), ["baseline", "+CompPC1", "+CompPC1+PC3", "+%mito", "matched"], rotation=45, ha="right", fontsize=8)
for i in range(H.shape[0]):
    for j, m in enumerate(["baseline", "comp_adjusted", "comp13_adjusted", "tech_adjusted"]):
        if Ri[f"{m}_p"].iloc[i] < 0.05:
            ax.text(j, i, "*", ha="center", va="center", fontsize=9, color="white" if abs(H.values[i, j]) > 1.5 else INK)
ax.grid(False)
fig.colorbar(im, ax=ax, shrink=0.4, label="effect (SD); * p<0.05")
ax2 = axes[1]
ax2.axis("off")
tc = {"Tier 1": "#0ca30c", "Tier 2": "#2a78d6", "Tier 3": "#eda100", "Tier 4": "#d03b3b", "no baseline": "#9a9994", "technical": "#4a3aa7"}
for i, p in enumerate(order):
    t = Ri.loc[p, "confidence_tier"]
    c = next((v for k, v in tc.items() if t.startswith(k) or k in t), INK)
    ax2.text(0, i, t, va="center", fontsize=7.5, color=c, transform=ax2.get_yaxis_transform())
ax2.set_ylim(len(order) - 0.5, -0.5)
ax2.set_title("Confidence tier", fontsize=10)
fig.suptitle("Pathway robustness summary (WT_NTC vs 5xFAD_NTC)", fontweight="bold")
save(fig, os.path.join(OUT, "Fig12_pathway_robustness_summary.png"))
print(R.confidence_tier.value_counts())
