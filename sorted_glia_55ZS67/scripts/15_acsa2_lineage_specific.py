"""15 - ACSA-2+ fraction: astrocyte-only and oligodendrocyte-only genes, SPP1 vs NTC, and whether astrocyte /
oligodendrocyte changes track each mouse's own microglial Spp1 knockdown (fractions are paired by animal; confirmed).

Lineage-only genes: for each gene, the share of its ACSA-2 signal expected from each cell type =
  (mean ACSA-2 RNA share of that type, step 08b deconvolution) x (brain reference expression, cp10k), normalised.
  astrocyte-only = >= 90% expected from astrocytes; oligodendrocyte-only = >= 90% from oligodendrocytes
  (80% used as a sensitivity check). Genes need mean ACSA-2 CPM >= 5 and must be protein-coding with a real symbol:
  pseudogenes / Gm / Rik loci track the arm-confounded library-quality axis (step 11) and are excluded.
Per-lineage DE (PyDESeq2, ~ treatment, n = 3 v 3): size factors fitted on that lineage's own genes (control_genes),
  so a change in how much astrocyte vs oligodendrocyte RNA a sort captured does not appear as differential
  expression. Consequence: a uniform up/down shift of the whole lineage is absorbed by normalisation; only changes
  relative to the lineage's typical gene are detectable. The lineage-content shift itself is reported separately.
Microglial link: per mouse, log2 CPM of Spp1 in that mouse's CD11b+ sample. Because treatment lowers Spp1, a
  knockdown-driven gene should (i) change with treatment and (ii) also track Spp1 within each arm. Test: across
  lineage genes, Spearman correlation between the treatment log2FC and the within-arm slope on Spp1; null from all
  36 within-arm permutations of the Spp1 values.
Outputs: 15_acsa2_lineage/ (tables, figures).
"""
import itertools
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "15_acsa2_lineage")
os.makedirs(OUT, exist_ok=True)
np.random.seed(1)
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 1.2, "savefig.dpi": 300, "axes.titleweight": "bold", "legend.frameon": False})
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
COL = {"NTC": "#eb6834", "SPP1": "#1baf7a"}
LINCOL = {"astrocyte": "#7a3fbf", "oligodendrocyte": "#2a78d6"}

# ---------------------------------------------------------------- data
m = pd.read_csv(os.path.join(ROOT, "data", "55ZS67-expression-matrix.tsv.gz"), sep="\t").set_index("gene_id")
m = m[~m.index.str.startswith("ERCC")]
sym = m.gene_name
cnt = m[[c for c in m.columns if c.endswith("_count")]].round().astype(int)
cnt.columns = [c[:-6] for c in cnt.columns]
ss = pd.read_csv(os.path.join(ROOT, "tables", "sample_sheet.tsv"), sep="\t", index_col=0)
ss["arm"] = np.where(ss.treatment == "Spp1", "SPP1", "NTC")
AC = ss.index[ss.fraction == "ACSA2"]
CD = ss.index[ss.fraction == "CD11b"]
cpm = cnt / cnt.sum() * 1e6
# mouse k: CD11b sample with same arm + replicate
pair = {a: ss.index[(ss.fraction == "CD11b") & (ss.arm == ss.loc[a, "arm"]) & (ss.replicate == ss.loc[a, "replicate"])][0] for a in AC}
arm = ss.loc[AC, "arm"].values
spp = (arm == "SPP1").astype(int)
mg_spp1 = np.log2(np.array([cpm.loc[sym[sym == "Spp1"].index[0], pair[a]] for a in AC]) + 1)

# ---------------------------------------------------------------- 1. lineage-only genes
ref = pd.read_csv(os.path.join(ROOT, "data", "reference", "brain_pseudobulk_cp10k.csv.gz"), index_col=0)
dec = pd.read_csv(os.path.join(ROOT, "02_purity_celltyping", "deconvolution_rna_share_pct_median.tsv"), sep="\t", index_col=0)
w = dec.loc[AC].mean() / 100
w = w[w > 0]
contrib = ref[w.index].mul(w, axis=1)
origin = contrib.div(contrib.sum(1) + 1e-12, axis=0)
expressed = cpm.loc[:, AC].mean(1) >= 5
gl = pd.DataFrame({"gene": sym.reindex(origin.index), "astro_origin": origin.astrocyte, "oligo_origin": origin.oligodendrocyte,
                   "ependymal_origin": origin.ependymal, "ACSA2_mean_cpm": cpm.loc[:, AC].mean(1).reindex(origin.index)})
biot = m.gene_biotype.reindex(gl.index)
real = (biot == "protein_coding") & ~gl.gene.fillna("").str.match(r"^(Gm\d|.*Rik$|.*-ps\d*$)")
gl = gl[expressed.reindex(gl.index).fillna(False) & real]
LIN = {}
for lin, colname in [("astrocyte", "astro_origin"), ("oligodendrocyte", "oligo_origin")]:
    for thr in (0.9, 0.8):
        LIN[(lin, thr)] = list(gl.index[gl[colname] >= thr])
gl["astrocyte_only_90"] = gl.index.isin(LIN[("astrocyte", 0.9)])
gl["oligo_only_90"] = gl.index.isin(LIN[("oligodendrocyte", 0.9)])
gl.to_csv(os.path.join(OUT, "ACSA2_gene_origin.tsv"), sep="\t", float_format="%.4g")
print({f"{k[0]}>={k[1]}": len(v) for k, v in LIN.items()})


# ---------------------------------------------------------------- 2. lineage content per sample
def sf_on(genes, samples):
    x = cnt.loc[genes, samples]
    x = x[(x > 0).all(1)]
    lg = np.log(x)
    return np.exp((lg.sub(lg.mean(1), axis=0)).median())


lib_sf = sf_on(list(cnt.index[(cnt[AC] > 0).all(1)]), AC)
content = {}
for lin in ("astrocyte", "oligodendrocyte"):
    content[lin] = np.log2(sf_on(LIN[(lin, 0.9)], AC) / lib_sf)  # lineage RNA relative to whole-library median
ep = list(gl.index[gl.ependymal_origin >= 0.9])
content["ependymal"] = np.log2(sf_on(ep, AC) / lib_sf)
crow = []
for lin, v in content.items():
    v = np.asarray(v)
    crow.append({"lineage": lin, "log2FC_SPP1_vs_NTC": v[spp == 1].mean() - v[spp == 0].mean(),
                 "welch_p": stats.ttest_ind(v[spp == 1], v[spp == 0], equal_var=False).pvalue,
                 **{f"{a}": round(float(x), 3) for a, x in zip(AC, v)}})
content_tab = pd.DataFrame(crow)
content_tab.to_csv(os.path.join(OUT, "lineage_content_per_sample.tsv"), sep="\t", index=False, float_format="%.4g")
print(content_tab.round(3).to_string())


# ---------------------------------------------------------------- 3. per-lineage DE
def lineage_de(genes, label):
    keep = cnt.loc[:, AC]
    keep = keep[(keep >= 10).sum(1) >= 3]
    meta = pd.DataFrame({"treatment": np.where(spp == 1, "Spp1", "NTC")}, index=AC)
    ctrl = [g for g in genes if g in keep.index]
    dds = DeseqDataSet(counts=keep.T, metadata=meta, design="~treatment", control_genes=ctrl, quiet=True, n_cpus=4)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["treatment", "Spp1", "NTC"], quiet=True)
    st.summary()
    r = st.results_df.loc[ctrl].copy()
    r["q_BH_within_lineage"] = multipletests(r.pvalue.fillna(1), method="fdr_bh")[1]
    r.insert(0, "gene", sym.reindex(r.index).values)
    r["normalisation"] = label
    sf = dds.obsm["size_factors"] if "size_factors" in dds.obsm else dds.obs["size_factors"]
    return r.sort_values("pvalue"), pd.Series(np.asarray(sf).ravel(), index=AC)


res, sfs = {}, {}
for lin in ("astrocyte", "oligodendrocyte"):
    for thr in (0.9, 0.8):
        r, sf = lineage_de(LIN[(lin, thr)], f"per-{lin}")
        res[(lin, thr)], sfs[(lin, thr)] = r, sf
        r.to_csv(os.path.join(OUT, f"DE_{lin}_only_genes_per_{lin}_{int(thr * 100)}.tsv"), sep="\t", float_format="%.4g")
        print(f"\n== {lin}-only (>= {thr}): {len(r)} genes; p<0.01: {(r.pvalue < 0.01).sum()}, q<0.1: {(r.q_BH_within_lineage < 0.1).sum()}, "
              f"q<0.05: {(r.q_BH_within_lineage < 0.05).sum()}")
        if thr == 0.9:
            print(r.head(15)[["gene", "baseMean", "log2FoldChange", "pvalue", "q_BH_within_lineage"]].round(4).to_string())

# library-normalised (standard) log2FC for the same genes, to show what composition does
std_de = pd.read_csv(os.path.join(ROOT, "06_DE_ACSA2", "DE_Spp1_vs_NTC_ACSA2.tsv"), sep="\t").set_index("gene_id")

# ---------------------------------------------------------------- 4. programs per host cell (lineage-attributable genes)
gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(os.path.join(ROOT, "data", "modules.gmt"))}
OLIGO_PROGRAMS = {
    "DOL_disease_oligo": ["Serpina3n", "C4b", "H2-D1", "H2-K1", "B2m", "Il33", "Klk6", "Cdkn1a", "Psmb8", "Stat1"],
    "Myelin_structural": ["Plp1", "Mbp", "Mog", "Mag", "Mobp", "Cldn11", "Opalin", "Ermn", "Mal", "Cnp", "Pllp", "Tspan2"],
    "Oligo_cholesterol": ["Hmgcr", "Hmgcs1", "Fdft1", "Sqle", "Cyp51", "Dhcr24", "Dhcr7", "Msmo1", "Idi1", "Ldlr", "Insig1", "Srebf2"],
}
prog_rows = []
name2id = {v: k for k, v in sym.items() if k in cnt.index}
for lin, sets, origin_col in [("astrocyte", {k: v for k, v in gmt.items() if k.startswith("A_")}, "astro_origin"),
                              ("oligodendrocyte", OLIGO_PROGRAMS, "oligo_origin")]:
    sf = sfs[(lin, 0.9)]
    for s, genes in sets.items():
        ids = [name2id[g] for g in genes if g in name2id and name2id[g] in gl.index and gl.loc[name2id[g], origin_col] >= 0.5]
        if len(ids) < 3:
            continue
        y = np.log2(cnt.loc[ids, AC].div(sf, axis=1) + 1).mean()
        d = y[spp == 1].mean() - y[spp == 0].mean()
        prog_rows.append({"lineage": lin, "program": s, "n_lineage_attributable_genes": len(ids), "genes": ";".join(sym[ids]),
                          "pct_change_per_host_cell": (2 ** d - 1) * 100,
                          "welch_p": stats.ttest_ind(y[spp == 1], y[spp == 0], equal_var=False).pvalue,
                          "r_with_microglial_Spp1_within_arm": stats.pearsonr(y - pd.Series(y).groupby(arm).transform("mean").values,
                                                                              mg_spp1 - pd.Series(mg_spp1).groupby(arm).transform("mean").values)[0]})
prog = pd.DataFrame(prog_rows)
prog["q_BH"] = multipletests(prog.welch_p, method="fdr_bh")[1]
prog.to_csv(os.path.join(OUT, "programs_per_host_cell.tsv"), sep="\t", index=False, float_format="%.4g")
print("\n", prog.drop(columns="genes").round(3).to_string())


# ---------------------------------------------------------------- 5. link to each mouse's microglial knockdown
def within(v):
    return v - pd.Series(v).groupby(arm).transform("mean").values


perms = []
idx = {a: np.where(arm == a)[0] for a in ("NTC", "SPP1")}
for p1 in itertools.permutations(idx["NTC"]):
    for p2 in itertools.permutations(idx["SPP1"]):
        o = np.arange(len(arm))
        o[idx["NTC"]] = list(p1)
        o[idx["SPP1"]] = list(p2)
        perms.append(o)
link_rows = []
mg_w = within(mg_spp1)
for lin in ("astrocyte", "oligodendrocyte"):
    r = res[(lin, 0.9)]
    sf = sfs[(lin, 0.9)]
    Y = np.log2(cnt.loc[r.index, AC].div(sf, axis=1) + 1).values
    Yw = Y - np.column_stack([Y[:, arm == a].mean(1) if False else None for a in []]) if False else np.array([within(y) for y in Y])

    def concord(xw):
        slope = (Yw @ xw) / (xw @ xw)  # pooled within-arm slope per gene
        return stats.spearmanr(r.log2FoldChange.values, slope)[0], slope

    obs, slope = concord(mg_w)
    null = np.array([concord(mg_w[o])[0] for o in perms])
    p_lo = float(np.mean(null <= obs + 1e-12))
    link_rows.append({"lineage": lin, "n_genes": len(r), "spearman_treatmentFC_vs_withinarm_slope": obs,
                      "expected_sign_if_knockdown_driven": "negative (lower Spp1 -> same direction as SPP1 effect)",
                      "perm_p_one_sided_negative": p_lo, "n_within_arm_permutations": len(perms)})
    r["within_arm_slope_per_log2_microglial_Spp1"] = slope
    r["r_total_with_microglial_Spp1"] = [stats.pearsonr(y, mg_spp1)[0] for y in Y]
    r.to_csv(os.path.join(OUT, f"DE_{lin}_only_genes_per_{lin}_90.tsv"), sep="\t", float_format="%.4g")
link = pd.DataFrame(link_rows)
link.to_csv(os.path.join(OUT, "microglial_knockdown_link.tsv"), sep="\t", index=False, float_format="%.4g")
print("\n", link.to_string())
print("microglial Spp1 log2CPM per mouse:", dict(zip([f"{a}:{x}" for a, x in zip(AC, arm)], mg_spp1.round(2))))


# ---------------------------------------------------------------- 6. figures
def fmt(p):
    return f"{p:.1e}" if p < 1e-3 else f"{p:.3f}" if p < 0.1 else f"{p:.2f}"


# 6a. lineage content
fig, axs = plt.subplots(1, 3, figsize=(9.5, 3.8))
for ax, (_, rr) in zip(axs, content_tab.iterrows()):
    v = np.array([rr[a] for a in AC])
    for i, a in enumerate(("NTC", "SPP1")):
        vv = v[arm == a]
        ax.scatter(np.full(3, i) + np.linspace(-0.12, 0.12, 3), vv, s=40, color=COL[a], edgecolor="white", zorder=3)
        ax.errorbar(i, vv.mean(), yerr=vv.std(ddof=1) / np.sqrt(3), color=INK, capsize=5, lw=1.2)
        ax.plot([i - 0.25, i + 0.25], [vv.mean()] * 2, color=INK, lw=1.6)
    ax.set_xticks([0, 1], ["NTC", "SPP1"]); ax.set_xlim(-0.6, 1.6)
    ax.set_title(f"{rr.lineage} RNA content\nWelch p={fmt(rr.welch_p)}; FC={2 ** rr.log2FC_SPP1_vs_NTC:.2f}", fontsize=10)
    ax.set_ylabel("log2 lineage content\n(rel. to library)" if ax is axs[0] else "")
fig.suptitle("ACSA-2+ fraction: how much of each lineage each sort captured", fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "15a_lineage_content.png"), bbox_inches="tight", facecolor="white"); plt.close(fig)

# 6b. volcanoes per lineage
fig, axs = plt.subplots(1, 2, figsize=(11, 4.8))
for ax, lin in zip(axs, ("astrocyte", "oligodendrocyte")):
    r = res[(lin, 0.9)]
    sig = r.q_BH_within_lineage < 0.1
    ax.scatter(r.log2FoldChange, -np.log10(r.pvalue), s=12, color=np.where(sig, LINCOL[lin], "#c9c8c3"), edgecolor="none")
    for _, rr in r.head(10).iterrows():
        ax.annotate(rr.gene, (rr.log2FoldChange, -np.log10(rr.pvalue)), fontsize=7.5, style="italic", xytext=(3, 2), textcoords="offset points")
    ax.axvline(0, color=INK2, lw=0.8)
    ax.set_xlabel("log2FC SPP1 vs NTC (per-" + lin + " normalised)")
    ax.set_ylabel("−log10 p (DESeq2)")
    ax.set_title(f"{lin}-only genes (n = {len(r)}); q<0.1: {int(sig.sum())}", fontsize=10)
fig.suptitle("Lineage-only genes in the ACSA-2+ fraction, SPP1 vs NTC (n = 3 v 3)", fontweight="bold")
fig.text(0.01, -0.02, "Each lineage normalised to its own genes (DESeq2 control genes). Coloured = BH q < 0.1 within that lineage's genes. "
         "A uniform shift of the whole lineage is absorbed by this normalisation (see 15a).", fontsize=7.5, color=INK2)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "15b_lineage_only_volcano.png"), bbox_inches="tight", facecolor="white"); plt.close(fig)

# 6c. microglial Spp1 vs lineage program score per mouse (top program per lineage + top gene)
fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
for ax, lin in zip(axs, ("astrocyte", "oligodendrocyte")):
    r = res[(lin, 0.9)]
    sf = sfs[(lin, 0.9)]
    up = r[(r.log2FoldChange > 0) & (r.pvalue < 0.05)].index
    dn = r[(r.log2FoldChange < 0) & (r.pvalue < 0.05)].index
    Y = np.log2(cnt.loc[r.index, AC].div(sf, axis=1) + 1)
    Z = Y.sub(Y.mean(1), axis=0).div(Y.std(1) + 1e-9, axis=0)
    score = (Z.loc[up].mean() if len(up) else 0) - (Z.loc[dn].mean() if len(dn) else 0)
    for a in ("NTC", "SPP1"):
        ax.scatter(mg_spp1[arm == a], np.asarray(score)[arm == a], s=60, color=COL[a], edgecolor="white", label=a, zorder=3)
    for i, a in enumerate(AC):
        ax.annotate(f"mouse {ss.loc[a, 'replicate']}", (mg_spp1[i], np.asarray(score)[i]), fontsize=7, xytext=(4, 3), textcoords="offset points", color=INK2)
    lk = link[link.lineage == lin].iloc[0]
    ax.set_xlabel("microglial Spp1 in the same mouse (log2 CPM, CD11b+ sort)")
    ax.set_ylabel("SPP1-response score\n(up − down genes, mean z)")
    ax.set_title(f"{lin}-only genes: {len(up)} up / {len(dn)} down (p<0.05)\n"
                 f"gene-wide concordance ρ={lk.spearman_treatmentFC_vs_withinarm_slope:.2f}, perm p={fmt(lk.perm_p_one_sided_negative)}", fontsize=9.5)
    ax.legend(fontsize=8, loc="best")
fig.tight_layout()
fig.text(0.01, -0.04, "Score built from the same 3 v 3 comparison, so the between-arm separation is circular; the informative part is whether mice "
         "WITHIN an arm line up with their own microglial Spp1. Concordance ρ: across lineage-only genes, treatment log2FC vs within-arm slope on "
         "microglial Spp1 (negative = knockdown-consistent); p from all 36 within-arm permutations.", fontsize=7.3, color=INK2, wrap=True)
fig.savefig(os.path.join(OUT, "15c_microglial_knockdown_link.png"), bbox_inches="tight", facecolor="white"); plt.close(fig)
print("figures written")


# ---------------------------------------------------------------- 7. sensitivity: leave one mouse out, and drop stress genes
STRESS = {"Sgk1", "Fkbp5", "Tsc22d3", "Klf9", "Per1", "Mt2", "Cdkn1a", "Zbtb16", "Klf15", "Errfi1", "Sult1a1", "Pdk4", "Arrdc3", "Arrdc2"}
loo = []
for lin in ("astrocyte", "oligodendrocyte"):
    r = res[(lin, 0.9)]
    sf = sfs[(lin, 0.9)]
    Yall = np.log2(cnt.loc[r.index, AC].div(sf, axis=1) + 1)
    for drop in [None] + list(AC) + ["stress_genes"]:
        keep = [a for a in AC if a != drop]
        Y = Yall.loc[[g for g in Yall.index if not (drop == "stress_genes" and sym[g] in STRESS)], keep].values
        ar = ss.loc[keep, "arm"].values
        x = mg_spp1[[list(AC).index(a) for a in keep]]
        fc = Y[:, ar == "SPP1"].mean(1) - Y[:, ar == "NTC"].mean(1)
        cen = lambda v: v - pd.Series(v).groupby(ar).transform("mean").values
        Yw = np.array([cen(y) for y in Y])
        def rho(xv):
            xw = cen(xv)
            return stats.spearmanr(fc, (Yw @ xw) / (xw @ xw))[0]
        obs = rho(x)
        ia = {a: np.where(ar == a)[0] for a in ("NTC", "SPP1")}
        null = []
        for p1 in itertools.permutations(ia["NTC"]):
            for p2 in itertools.permutations(ia["SPP1"]):
                o = np.arange(len(ar)); o[ia["NTC"]] = list(p1); o[ia["SPP1"]] = list(p2)
                null.append(rho(x[o]))
        lab = "none" if drop is None else drop if drop == "stress_genes" else f"{drop} ({ss.loc[drop, 'arm']} mouse {ss.loc[drop, 'replicate']})"
        loo.append({"lineage": lin, "dropped": lab, "rho": obs, "perm_p_one_sided": float(np.mean(np.array(null) <= obs + 1e-12)), "n_perm": len(null)})
loo = pd.DataFrame(loo)
loo.to_csv(os.path.join(OUT, "microglial_knockdown_link_leave_one_out.tsv"), sep="\t", index=False, float_format="%.4g")
print("\n", loo.round(3).to_string())


# ---------------------------------------------------------------- 8. is within-arm Spp1 knockdown, or disease severity?
ki = pd.read_csv(os.path.join(ROOT, "09_pathology", "spp1_knockdown_index.tsv"), sep="\t", index_col=0)
qc = pd.read_csv(os.path.join(ROOT, "01_QC", "qc_metrics.tsv"), sep="\t", index_col=0)
per_mouse = pd.DataFrame({"arm": arm, "replicate": ss.loc[AC, "replicate"].values,
                          "mg_Spp1_log2cpm": mg_spp1,
                          "mg_DAM2_log2": [ki.loc[pair[a], "DAM2_mean_log2"] for a in AC],
                          "mg_Spp1_minus_DAM2": [ki.loc[pair[a], "Spp1_minus_DAM2"] for a in AC],
                          "ACSA2_astro_content": np.asarray(content["astrocyte"]),
                          "ACSA2_pct_unmapped": qc.loc[AC, "pct_unmapped"].values,
                          "CD11b_pct_unmapped": [qc.loc[pair[a], "pct_unmapped"] for a in AC]}, index=AC)
per_mouse.to_csv(os.path.join(OUT, "per_mouse_covariates.tsv"), sep="\t", float_format="%.4g")
wc = {c: stats.pearsonr(within(per_mouse.mg_Spp1_log2cpm.values), within(per_mouse[c].values))[0]
      for c in ["mg_DAM2_log2", "mg_Spp1_minus_DAM2", "ACSA2_astro_content", "ACSA2_pct_unmapped", "CD11b_pct_unmapped"]}
print("\nwithin-arm r of microglial Spp1 with:", {k: round(v, 2) for k, v in wc.items()})
alt = []
for lin in ("astrocyte", "oligodendrocyte"):
    r = res[(lin, 0.9)]
    sf = sfs[(lin, 0.9)]
    Y = np.log2(cnt.loc[r.index, AC].div(sf, axis=1) + 1).values
    Yw = np.array([within(y) for y in Y])
    fc = Y[:, arm == "SPP1"].mean(1) - Y[:, arm == "NTC"].mean(1)
    for xname in ["mg_Spp1_log2cpm", "mg_Spp1_minus_DAM2", "mg_DAM2_log2", "ACSA2_pct_unmapped"]:
        x = per_mouse[xname].values
        def rho(xv):
            xw = within(xv)
            return stats.spearmanr(fc, (Yw @ xw) / (xw @ xw))[0]
        obs = rho(x)
        null = np.array([rho(x[o]) for o in perms])
        alt.append({"lineage": lin, "per_mouse_x": xname, "rho": obs, "perm_p_one_sided_negative": float(np.mean(null <= obs + 1e-12))})
alt = pd.DataFrame(alt)
alt.to_csv(os.path.join(OUT, "microglial_link_knockdown_vs_severity.tsv"), sep="\t", index=False, float_format="%.4g")
print(per_mouse.round(3).to_string())
print(alt.round(3).to_string())
