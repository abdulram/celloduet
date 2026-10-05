"""07 - Cell-STATE vs cell-ABUNDANCE (inferred; not single-cell evidence).

For each lineage with a VALIDATED abundance signature (step 02), functional state modules (gene_sets.tsv) are
scored after REMOVING any gene that is also an abundance marker of that lineage (no circularity).
Model (WT_NTC + 5xFAD_NTC, n = 12):
    state_score ~ lineage_abundance + sex + genotype                 (primary)
    state_score ~ lineage_abundance + sex + pct_mito + genotype      (technical sensitivity)
The genotype coefficient = genotype difference in state at equal lineage abundance.
p-values: OLS and exact permutation (genotype relabelled within sex, 20 x 20 = 400 relabellings).
Lineages whose abundance signature failed validation (T cell, lymphatic endothelial, plasma cell, mural) cannot be
analysed this way and are reported as such.
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
OUT = os.path.join(CA, "07_cell_state_analysis")
TAB = os.path.join(CA, "tables")
exec(open(os.path.join(HERE, "common.py")).read())

meta = pd.read_csv(os.path.join(TAB, "sample_metadata_with_composition.tsv"), sep="\t", index_col=0)
vst = pd.read_csv(os.path.join(CA, "01_QC", "vst_blind.tsv.gz"), sep="\t", index_col=0)
S = pd.read_csv(os.path.join(TAB, "composition_signature_scores_all.tsv"), sep="\t", index_col=0)
val = pd.read_csv(os.path.join(TAB, "signature_validation.tsv"), sep="\t").set_index("signature")
gs = pd.read_csv(os.path.join(TAB, "gene_sets.tsv"), sep="\t")
items = pd.read_csv(os.path.join(TAB, "signature_marker_item_rest.tsv"), sep="\t")

LINEAGE_STATES = {
    "macrophage_BAM": ["CURATED_INNATE_IMMUNE_ACTIVATION", "CURATED_TYPE_I_IFN_ISG", "CURATED_MHCII_ANTIGEN_PRESENTATION",
                       "CURATED_COMPLEMENT_CLASSICAL", "CURATED_PHAGOCYTOSIS_LYSOSOME", "CURATED_LIPID_HANDLING", "CURATED_DAM_KerenShaul2017"],
    "fibroblast_stromal": ["CURATED_ECM_REMODELING", "CURATED_FIBROBLAST_ACTIVATION", "CURATED_INFLAMMATORY_FIBROBLAST"],
    "blood_endothelial": ["CURATED_VASCULAR_INFLAMMATION_ADHESION", "CURATED_ANGIOGENESIS_ENDOTHELIAL", "CURATED_MENINGEAL_IMMUNE_TRAFFICKING"],
    "T_cell": ["CURATED_T_CELL_ACTIVATION", "CURATED_T_CELL_CYTOTOXICITY"],
    "lymphatic_endothelial": ["CURATED_LYMPHATIC_ENDOTHELIAL"],
    "B_cell": ["CURATED_B_PLASMA"],
    "monocyte": ["CURATED_MONOCYTE_RECRUITMENT"],
}
P = meta[meta.group.isin(["WT_NTC", "5xFAD_NTC"])].copy()
P["g5"] = (P.genotype == "5xFAD").astype(int)
zv = vst.sub(vst.mean(axis=1), axis=0).div(vst.std(axis=1).replace(0, np.nan), axis=0)
name2id = pd.Series(gs.gene_id.values, index=gs.gene_name.values)
fem, mal = P.index[P.sex == "F"], P.index[P.sex == "M"]
perms = [set(a) | set(b) for a in itertools.combinations(fem, 3) for b in itertools.combinations(mal, 3)]


def perm_p(df, formula, obs):
    null = []
    for s in perms:
        d = df.copy()
        d["g5"] = d.index.isin(list(s)).astype(int)
        null.append(smf.ols(formula, d).fit().params["g5"])
    null = np.array(null)
    return (np.abs(null) >= abs(obs) - 1e-12).mean()


rows = []
state_scores = {}
for lin, states in LINEAGE_STATES.items():
    ok = bool(val.loc[lin, "validated"])
    abund_markers = set(items[items.signature == lin].gene)
    for st in states:
        ids = gs[(gs.set == st) & (~gs.gene_name.isin(abund_markers))].gene_id
        ids = [i for i in ids if i in vst.index]
        if len(ids) < 4:
            continue
        sc = zv.loc[ids].mean()
        state_scores[f"{lin}|{st}"] = sc
        row = {"lineage": lin, "abundance_signature_validated": ok, "state_module": st, "n_state_genes": len(ids)}
        d = P.assign(state=sc.loc[P.index], ab=S.loc[P.index, lin])
        row["r_state_vs_abundance"] = np.corrcoef(d.state, d.ab)[0, 1]
        f0 = smf.ols("state ~ C(sex) + g5", d).fit()
        row["genotype_beta_unadjusted"], row["p_unadjusted"] = f0.params["g5"], f0.pvalues["g5"]
        if ok:
            f1 = smf.ols("state ~ ab + C(sex) + g5", d).fit()
            f2 = smf.ols("state ~ ab + C(sex) + pct_mito + g5", d).fit()
            row["genotype_beta_given_abundance"], row["p_given_abundance"] = f1.params["g5"], f1.pvalues["g5"]
            row["perm_p_given_abundance"] = perm_p(d, "state ~ ab + C(sex) + g5", f1.params["g5"])
            row["genotype_beta_given_abundance_and_mito"], row["p_given_abundance_and_mito"] = f2.params["g5"], f2.pvalues["g5"]
            row["abundance_beta"] = f1.params["ab"]
        else:
            row["note"] = "abundance signature not validated - state cannot be separated from abundance"
        rows.append(row)
res = pd.DataFrame(rows)
res.to_csv(os.path.join(TAB, "cell_state_vs_abundance.tsv"), sep="\t", index=False)
pd.DataFrame(state_scores).to_csv(os.path.join(TAB, "cell_state_scores.tsv"), sep="\t")
print(res.drop(columns=[c for c in ["note"] if c in res]).round(3).to_string())

# Figure 9: state vs abundance scatter per validated lineage x state
vr = res[res.abundance_signature_validated]
n = len(vr)
ncol = 4
fig, axes = plt.subplots(int(np.ceil(n / ncol)), ncol, figsize=(16, 3.6 * int(np.ceil(n / ncol))))
for ax, (_, r) in zip(axes.flat, vr.iterrows()):
    sc = state_scores[f"{r.lineage}|{r.state_module}"]
    for g in ["WT_NTC", "5xFAD_NTC", "5xFAD_SPP1"]:
        idx = meta.index[meta.group == g]
        ax.scatter(S.loc[idx, r.lineage], sc[idx], c=GROUP_COLORS[g], s=34, edgecolor=SURFACE, lw=1, alpha=1 if g != "5xFAD_SPP1" else 0.35)
    ax.set_xlabel(f"{r.lineage} abundance", fontsize=8)
    ax.set_ylabel("state score", fontsize=8)
    ax.set_title(f"{r.state_module.replace('CURATED_', '')}\nβ_geno|abund={r.genotype_beta_given_abundance:.2f}, perm p={r.perm_p_given_abundance:.2f}", fontsize=8)
for ax in list(axes.flat)[n:]:
    ax.axis("off")
fig.suptitle("Inferred cell state vs lineage abundance (WT blue, 5xFAD-NTC orange; SPP1 faded, not modelled)", fontweight="bold")
fig.tight_layout()
save(fig, os.path.join(OUT, "Fig9_state_vs_abundance.png"))
