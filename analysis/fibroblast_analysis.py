"""Fibroblast-focused re-analysis to minimise contamination from other cell types.

Strategy
  1. Reference-based deconvolution (NNLS on class marker genes) of each sample against a CELLxGENE Census
     pseudobulk reference (analysis/reference/build_reference.py) to estimate cell-type composition.
  2. Fibroblast signal share per gene: using each sample's estimated composition, the fraction of the gene's
     predicted bulk signal that comes from fibroblast classes. Strict: median share >= 0.8 and >= 0.6 in every
     sample; relaxed: median >= 0.6 and >= 0.4 in every sample. Bulk sanity filter: drop genes whose expression
     still tracks a contaminant module (|r| >= 0.5), which catches cell types missing from the reference.
  3. (see 2)
  4. Fibroblast-normalised DE: DESeq2 on the fibroblast-specific genes only, so size factors come from
     fibroblast genes and effects are per unit of fibroblast content. Two models: ~sex + group and
     ~sex + pct_mito + group (technical axis).
  5. Rank-based enrichment within the fibroblast gene universe.

Run after run_analysis.py and reference/build_reference.py. Outputs: results/fibroblast/.
"""
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.optimize import nnls
from scipy.stats import mannwhitneyu, spearmanr
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import gseapy

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "results", "fibroblast")
os.makedirs(OUT, exist_ok=True)

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
GROUP_COLORS = {"WT_NTC": "#2a78d6", "5xFAD_NTC": "#eb6834", "5xFAD_SPP1": "#1baf7a"}
GROUP_ORDER = list(GROUP_COLORS)
SEX_MARKER = {"F": "o", "M": "^"}
DIV = LinearSegmentedColormap.from_list("bgr", ["#184f95", "#3987e5", "#86b6ef", "#f0efec", "#f0a3a2", "#e34948", "#a8322f"])
SEQ = LinearSegmentedColormap.from_list("blue", ["#f0efec", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
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


# ---------------------------------------------------------------- load bulk
qc = pd.read_csv(os.path.join(ROOT, "results", "qc_and_composition.csv"), index_col=0)
meta = qc[["sample", "genotype", "sex", "guide", "group"]].copy()
meta["animal"] = meta["sample"].str.split("_").str[0]
libs = list(qc.index)

raw = pd.read_csv(os.path.join(DATA, "VT3YXB-expression-matrix.tsv"), sep="\t")
raw = raw[~raw.gene_id.str.startswith("ERCC")]
raw["gene_name"] = raw["gene_name"].fillna(raw["gene_id"])
counts = raw.set_index("gene_id")[[f"{l}_count" for l in libs]].astype(int)
counts.columns = libs
cpm = raw.set_index("gene_id")[[f"{l}_cpm" for l in libs]]
cpm.columns = libs
g2name = raw.set_index("gene_id")["gene_name"]
g2type = raw.set_index("gene_id")["gene_biotype"]
expressed = ((counts >= 10).sum(axis=1) >= 6) & g2type.reindex(counts.index).isin(["protein_coding", "lncRNA"])
logc = np.log2(cpm + 1)

CONTAMINANTS = ["mod_osteoblast_bone", "mod_osteoclast", "mod_choroid_plexus", "mod_olfactory_mucosa",
                "mod_erythroid_blood", "mod_neutrophil_marrow", "mod_skeletal_muscle", "mod_brain_parenchyma"]

# ---------------------------------------------------------------- 1. reference specificity
ref = pd.read_csv(os.path.join(DATA, "reference", "census_pseudobulk_cp10k.csv.gz"), index_col=0)
ref_summary = pd.read_csv(os.path.join(DATA, "reference", "census_class_summary.csv"), index_col=0)
fib_cls = [c for c in ref.columns if c.startswith("FIB_")]
oth_cls = [c for c in ref.columns if not c.startswith("FIB_")]
ref = ref.reindex(counts.index).dropna(how="all").fillna(0)

EPS = 0.02
fib_max = ref[fib_cls].max(axis=1)
fib_best = ref[fib_cls].idxmax(axis=1)
oth_max = ref[oth_cls].max(axis=1)
oth_best = ref[oth_cls].idxmax(axis=1)
ratio = (fib_max + EPS) / (oth_max + EPS)
spec = pd.DataFrame({"gene_name": g2name.reindex(ref.index), "fib_max_cp10k": fib_max, "fib_class": fib_best,
                     "other_max_cp10k": oth_max, "other_class": oth_best, "specificity_ratio": ratio,
                     "bulk_expressed": expressed.reindex(ref.index).fillna(False)})
for c in fib_cls:
    spec[c] = ref[c]

# bulk contaminant check (Pearson on log2 CPM across the 18 samples)
def rowcorr(Y, v):
    Yc = Y.sub(Y.mean(1), axis=0)
    vc = v - v.mean()
    return (Yc @ vc) / (np.sqrt((Yc**2).sum(1)) * np.sqrt((vc**2).sum()) + 1e-12)


# ---------------------------------------------------------------- 2. deconvolution (NNLS on marker genes)
# markers per class: top 25 genes by (class / max other class) ratio, expressed >= 1 CP10K in that class
ref_e = ref.loc[ref.index.intersection(expressed.index[expressed])]
markers = {}
for c in ref.columns:
    others = ref_e.drop(columns=c).max(axis=1)
    r = (ref_e[c] + EPS) / (others + EPS)
    m = r[(ref_e[c] >= 1) & (r >= 3)].sort_values(ascending=False).head(25).index
    if len(m) >= 5:
        markers[c] = m
mk = pd.Index(sorted(set().union(*markers.values())))
S = ref.loc[mk, list(markers)] * 100  # CP10K -> CPM scale
B = cpm.loc[mk]
w = 1 / (S.mean(axis=1) + 1)  # down-weight very highly expressed markers
props = {}
for l in libs:
    coef, _ = nnls((S.values * w.values[:, None]), B[l].values * w.values)
    props[l] = coef
raw_coef = {l: v.copy() for l, v in props.items()}
props = pd.DataFrame(props, index=list(markers)).T
props = props.div(props.sum(axis=1), axis=0)
props.to_csv(os.path.join(OUT, "deconvolution_proportions.csv"))
fib_frac = props[[c for c in props if c.startswith("FIB_")]].sum(axis=1)
print("Deconvolved fibroblast fraction by group:\n", fib_frac.groupby(meta.group).median().round(3))

# ---------------------------------------------------------------- fibroblast signal share per gene
# Expected contribution of each reference class to each gene in each sample = NNLS coefficient x reference profile.
# fib_share = fraction of a gene's predicted bulk signal that comes from fibroblast classes. This weights
# every contaminant by how abundant it actually is in these samples (unlike a max-over-all-cell-types ratio).
coef = pd.DataFrame(raw_coef, index=list(markers)).T  # samples x classes (CPM scale)
R = ref[list(markers)] * 100
pred_total = R.values @ coef.T.values  # genes x samples
fcols = [i for i, c in enumerate(markers) if c.startswith("FIB_")]
pred_fib = R.values[:, fcols] @ coef.T.values[fcols]
share = pd.DataFrame(pred_fib / np.maximum(pred_total, 1e-9), index=R.index, columns=libs)
spec["fib_share_median"] = share.median(axis=1)
spec["fib_share_min"] = share.min(axis=1)
spec["top_contaminant_class"] = pd.DataFrame(R.values[:, [i for i in range(len(markers)) if i not in fcols]] * coef.T.values[[i for i in range(len(markers)) if i not in fcols]].mean(axis=1),
                                             index=R.index, columns=[c for c in markers if not c.startswith("FIB_")]).idxmax(axis=1)

cand = spec.index[(spec.fib_max_cp10k >= 0.3) & (spec.fib_share_median >= 0.6) & spec.bulk_expressed]
cont_r = pd.DataFrame({m: rowcorr(logc.loc[cand], qc[m]) for m in CONTAMINANTS})
spec["max_abs_r_contaminant_bulk"] = cont_r.abs().max(axis=1).reindex(spec.index)
spec["worst_contaminant_bulk"] = cont_r.abs().idxmax(axis=1).reindex(spec.index)
spec["tier"] = "not_specific"
ok = spec.index.isin(cand) & (spec.max_abs_r_contaminant_bulk < 0.5)
spec.loc[ok & (spec.fib_share_median >= 0.6) & (spec.fib_share_min >= 0.4), "tier"] = "relaxed"
spec.loc[ok & (spec.fib_share_median >= 0.8) & (spec.fib_share_min >= 0.6), "tier"] = "strict"
spec.loc[spec.index.isin(cand) & (spec.max_abs_r_contaminant_bulk >= 0.5), "tier"] = "dropped_bulk_contaminant"
spec.sort_values("fib_share_median", ascending=False).to_csv(os.path.join(OUT, "fibroblast_specificity_all_genes.csv.gz"))
fib_genes = spec.index[spec.tier.isin(["strict", "relaxed"])]
strict = spec.index[spec.tier == "strict"]
fs = spec.loc[fib_genes].sort_values("fib_share_median", ascending=False)
fs.to_csv(os.path.join(OUT, "fibroblast_specific_genes.csv"))
print(spec.tier.value_counts())
print("Top strict genes:", ", ".join(fs[fs.tier == "strict"].sort_values("fib_max_cp10k", ascending=False).gene_name.head(80).astype(str)))

# enrichment sanity check of the gene set
hall = {}
with open(os.path.join(DATA, "genesets", "h.all.v7.0.symbols.gmt")) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        hall[f[0]] = set(f[2:])
from scipy.stats import hypergeom
universe = set(g2name.reindex(expressed.index[expressed]).astype(str).str.upper())
fsym = set(fs.gene_name.astype(str).str.upper())
ora = []
for k, s in hall.items():
    s = s & universe
    x = len(s & fsym)
    ora.append({"set": k, "overlap": x, "set_size": len(s),
                "p": hypergeom.sf(x - 1, len(universe), len(s), len(fsym)) if x else 1.0})
ora = pd.DataFrame(ora).sort_values("p")
ora["fdr"] = np.minimum(1, ora.p * len(ora) / (np.arange(len(ora)) + 1))
ora.to_csv(os.path.join(OUT, "fibroblast_genes_hallmark_overlap.csv"), index=False)
print(ora.head(5).to_string())

# fibroblast content from specific genes (bulk)
strict_score = logc.loc[strict].sub(logc.loc[strict].mean(1), axis=0).mean()
qc["fibroblast_score"] = strict_score
qc["fibroblast_fraction_deconv"] = fib_frac

# ---------------------------------------------------------------- 3. fibroblast-normalised DE
fc = counts.loc[fib_genes]
md = meta[["sex", "group"]].copy()
md["group"] = pd.Categorical(md.group, categories=GROUP_ORDER)
md["pct_mito"] = ((qc.pct_mito - qc.pct_mito.mean()) / qc.pct_mito.std()).values
CONTRASTS = {"5xFAD_NTC_vs_WT_NTC": ("5xFAD_NTC", "WT_NTC"), "5xFAD_SPP1_vs_5xFAD_NTC": ("5xFAD_SPP1", "5xFAD_NTC")}
MODELS = {"fibnorm": "~sex + group", "fibnorm_mitoadj": "~sex + pct_mito + group"}
de, summary = {}, []
for mname, design in MODELS.items():
    dds = DeseqDataSet(counts=fc.T, metadata=md, design=design, refit_cooks=True, quiet=True)
    dds.deseq2()
    for cname, (a, b) in CONTRASTS.items():
        st = DeseqStats(dds, contrast=["group", a, b], quiet=True)
        st.summary()
        r = st.results_df.copy()
        r.insert(0, "gene_name", g2name.reindex(r.index).values)
        r["tier"] = spec.tier.reindex(r.index).values
        r["fib_class"] = spec.fib_class.reindex(r.index).values
        r["specificity_ratio"] = spec.specificity_ratio.reindex(r.index).values
        r = r.sort_values("pvalue")
        r.to_csv(os.path.join(OUT, f"DE_{cname}_{mname}.csv"))
        de[(cname, mname)] = r
        s = r[r.padj < 0.05]
        summary.append({"model": mname, "design": design, "contrast": cname, "genes_tested": len(r),
                        "padj<0.05": len(s), "up": int((s.log2FoldChange > 0).sum()),
                        "down": int((s.log2FoldChange < 0).sum()), "padj<0.1": int((r.padj < 0.1).sum())})
summary = pd.DataFrame(summary)
summary.to_csv(os.path.join(OUT, "DE_summary.csv"), index=False)
print(summary.to_string())

# whole-tissue vs fibroblast-normalised effect sizes for the same genes
for cname in CONTRASTS:
    base = pd.read_csv(os.path.join(ROOT, "results", f"DE_{cname}_base.csv"), index_col=0)
    j = de[(cname, "fibnorm")][["gene_name", "log2FoldChange", "padj"]].join(
        base[["log2FoldChange", "padj"]], rsuffix="_wholetissue", how="inner")
    j.to_csv(os.path.join(OUT, f"compare_fibnorm_vs_wholetissue_{cname}.csv"))

# ---------------------------------------------------------------- 4. enrichment within the fibroblast universe
sets = {k: sorted(v) for k, v in hall.items()}
custom = {
    "CUSTOM_INTERFERON_ISG": ["Ifit1", "Ifit3", "Isg15", "Irf7", "Oasl2", "Rsad2", "Usp18", "Stat1", "Bst2", "Ifi27l2a", "Mx1", "Oas1a", "Rtp4", "Ifitm3", "Xaf1"],
    "CUSTOM_FIBROBLAST_ACTIVATION": ["Postn", "Col1a1", "Col3a1", "Fn1", "Tnc", "Timp1", "Lox", "Loxl2", "Ctgf", "Serpine1", "Thbs1", "Acta2", "Cthrc1", "Sparc", "Col5a1", "Col12a1"],
    "CUSTOM_ARACHNOID_RA_BARRIER": ["Aldh1a2", "Crabp2", "Rbp1", "Cyp26b1", "Rdh10", "Slc47a1", "Cldn11", "Tjp1", "Slc22a6", "Dpp4", "Abcg2", "Slc6a13"],
}
sets.update({k: [g.upper() for g in v] for k, v in custom.items()})
gsea_out = {}
for key, r in de.items():
    rr = r.dropna(subset=["stat"])
    rnk = rr.assign(sym=rr.gene_name.astype(str).str.upper()).groupby("sym")["stat"].apply(
        lambda s: s.loc[s.abs().idxmax()]).sort_values(ascending=False)
    try:
        pre = gseapy.prerank(rnk=rnk, gene_sets=sets, min_size=5, max_size=500, permutation_num=2000,
                             seed=7, threads=4, outdir=None, verbose=False)
    except LookupError:
        print("GSEA skipped (no gene sets with >= 5 genes in universe):", key)
        gsea_out[key] = pd.DataFrame(columns=["Term", "NES", "FDR q-val"])
        continue
    g = pre.res2d.copy()
    g["NES"] = g.NES.astype(float)
    g["FDR q-val"] = g["FDR q-val"].astype(float)
    g = g.sort_values("NES", ascending=False)
    g.to_csv(os.path.join(OUT, f"GSEA_{key[0]}_{key[1]}.csv"), index=False)
    gsea_out[key] = g

# ---------------------------------------------------------------- figures
def scatter_groups(ax, x, y):
    for grp in GROUP_ORDER:
        for sx, mk_ in SEX_MARKER.items():
            idx = meta.index[(meta.group == grp) & (meta.sex == sx)]
            ax.scatter(x[idx], y[idx], c=GROUP_COLORS[grp], marker=mk_, s=70, edgecolor=SURFACE, lw=1.5, zorder=3)
    for l in libs:
        ax.annotate(meta.animal[l], (x[l], y[l]), xytext=(5, 4), textcoords="offset points", fontsize=7, color=INK2)


def legend(ax, loc="best"):
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="s", ls="", color=GROUP_COLORS[g], ms=9, label=g.replace("_", " ")) for g in GROUP_ORDER]
    h += [Line2D([], [], marker=m, ls="", color=INK2, ms=8, label={"F": "female", "M": "male"}[s]) for s, m in SEX_MARKER.items()]
    ax.legend(handles=h, loc=loc, fontsize=8)


# F1: reference heatmap of the top fibroblast-specific genes across all reference classes
top = fs[fs.tier == "strict"].sort_values("fib_max_cp10k", ascending=False).head(40).index
col_order = fib_cls + [c for c in oth_cls]
mat = np.log10(ref.loc[top, col_order] + 0.01)
fig, ax = plt.subplots(figsize=(14, 10))
im = ax.imshow(mat.values, cmap=SEQ, aspect="auto", vmin=-2, vmax=mat.values.max())
ax.set_yticks(range(len(top)), g2name.reindex(top).values, fontsize=7.5)
ax.set_xticks(range(len(col_order)), [f"{c} (n={ref_summary.cells.get(c, 0)})" for c in col_order], rotation=90, fontsize=7.5)
for t, c in zip(ax.get_xticklabels(), col_order):
    if c.startswith("FIB_"):
        t.set_fontweight("bold")
ax.axvline(len(fib_cls) - 0.5, color=INK, lw=1.5)
ax.grid(False)
fig.colorbar(im, ax=ax, shrink=0.5, label="log10 mean CP10K (single-cell pseudobulk)")
ax.set_title("40 most highly expressed strict fibroblast genes across the CELLxGENE mouse reference")
save(fig, "F1_reference_specificity_heatmap.png")

# F2: deconvolution
order = meta.sort_values(["group", "sex"]).index
pr = props.loc[order]
keep_cols = pr.columns[pr.max() >= 0.03]
pr_plot = pr[keep_cols].copy()
pr_plot["other (<3% max)"] = pr.drop(columns=keep_cols).sum(axis=1)
fig, axes = plt.subplots(1, 2, figsize=(17, 6), gridspec_kw={"width_ratios": [2.2, 1]})
ax = axes[0]
im = ax.imshow(pr_plot.T.values, cmap=SEQ, aspect="auto", vmin=0, vmax=max(0.5, pr_plot.values.max()))
ax.set_yticks(range(pr_plot.shape[1]), pr_plot.columns, fontsize=8)
ax.set_xticks(range(len(order)), [f"{meta.animal[l]} {meta.group[l].replace('_', ' ')} {meta.sex[l]}" for l in order], rotation=90, fontsize=7)
for t, l in zip(ax.get_xticklabels(), order):
    t.set_color(GROUP_COLORS[meta.group[l]])
for i in range(pr_plot.shape[1]):
    for j in range(pr_plot.shape[0]):
        v = pr_plot.values[j, i]
        if v >= 0.05:
            ax.text(j, i, f"{v*100:.0f}", ha="center", va="center", fontsize=6.5, color="white" if v > 0.3 else INK)
ax.grid(False)
fig.colorbar(im, ax=ax, shrink=0.6, label="estimated fraction of marker-gene signal")
ax.set_title("Reference-based deconvolution (NNLS on class marker genes), % shown where ≥5")
ax = axes[1]
for i, grp in enumerate(GROUP_ORDER):
    idx = meta.index[meta.group == grp]
    ax.scatter(np.full(len(idx), i) + np.linspace(-0.12, 0.12, len(idx)), fib_frac[idx] * 100, c=GROUP_COLORS[grp],
               s=60, edgecolor=SURFACE, lw=1.5, zorder=3)
    ax.hlines(fib_frac[idx].median() * 100, i - 0.25, i + 0.25, color=INK, lw=2)
p1 = mannwhitneyu(fib_frac[meta.group == "5xFAD_NTC"], fib_frac[meta.group == "WT_NTC"]).pvalue
p2 = mannwhitneyu(fib_frac[meta.group == "5xFAD_SPP1"], fib_frac[meta.group == "5xFAD_NTC"]).pvalue
ax.set_xticks(range(3), [g.replace("_", " ") for g in GROUP_ORDER])
ax.set_ylabel("Estimated fibroblast fraction (%)")
ax.set_title(f"Fibroblast fraction by group\n5xFAD vs WT p={p1:.3f}; SPP1 vs NTC p={p2:.2f}", fontsize=9.5)
ax.grid(axis="x", visible=False)
fig.tight_layout()
save(fig, "F2_deconvolution.png")

# F3: fibroblast content vs technical axis
fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
scatter_groups(axes[0], qc.pct_mito, strict_score)
axes[0].set_xlabel("% mitochondrial reads (technical axis)")
axes[0].set_ylabel("Fibroblast score (strict genes, centred log2 CPM)")
axes[0].set_title(f"Fibroblast signal vs technical axis (r = {np.corrcoef(qc.pct_mito, strict_score)[0, 1]:.2f})")
legend(axes[0], loc="lower left")
scatter_groups(axes[1], fib_frac * 100, strict_score)
axes[1].set_xlabel("Deconvolved fibroblast fraction (%)")
axes[1].set_ylabel("Fibroblast score (strict genes)")
axes[1].set_title(f"Two independent fibroblast estimates agree (ρ = {spearmanr(fib_frac, strict_score)[0]:.2f})")
fig.tight_layout()
save(fig, "F3_fibroblast_content_vs_technical.png")

# F4: volcanoes
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
for ax, (key, r) in zip(axes.flat, de.items()):
    r = r.dropna(subset=["padj"])
    y = -np.log10(r.pvalue.clip(lower=1e-300))
    sig = r.padj < 0.05
    ax.scatter(r.log2FoldChange[~sig], y[~sig], s=8, c="#c9c8c3", lw=0)
    for m, c, lab in [((r.log2FoldChange > 0) & sig, "#e34948", "up"), ((r.log2FoldChange < 0) & sig, "#2a78d6", "down")]:
        ax.scatter(r.log2FoldChange[m], y[m], s=16, c=c, lw=0, label=f"{lab} (n={m.sum()})")
    for gid, row in r.head(15).iterrows():
        ax.annotate(row.gene_name, (row.log2FoldChange, -np.log10(max(row.pvalue, 1e-300))), xytext=(3, 2),
                    textcoords="offset points", fontsize=7)
    ax.axvline(0, color=INK2, lw=0.6)
    ax.set_xlabel("log2 fold change (per unit fibroblast content)")
    ax.set_ylabel("-log10 p-value")
    ax.set_title(f"{key[0].replace('_vs_', ' vs ').replace('_', ' ')}\nmodel: {MODELS[key[1]]}", fontsize=9.5)
    ax.legend(loc="upper left", fontsize=8)
fig.suptitle(f"Fibroblast-normalised DE on {len(fib_genes)} fibroblast-specific genes; red/blue = padj < 0.05", fontweight="bold")
fig.tight_layout()
save(fig, "F4_volcano_fibroblast.png")

# F5: GSEA
fig, axes = plt.subplots(1, 4, figsize=(22, 7))
for ax, (key, g) in zip(axes, gsea_out.items()):
    d = g[g["FDR q-val"] < 0.25]
    d = pd.concat([d.head(10), d.tail(10)]).drop_duplicates("Term").sort_values("NES")
    if d.empty:
        ax.text(0.5, 0.5, "No gene sets at FDR < 0.25", ha="center", transform=ax.transAxes)
    else:
        ax.barh(range(len(d)), d.NES, color=["#e34948" if v > 0 else "#2a78d6" for v in d.NES], height=0.75,
                edgecolor=SURFACE, lw=1.5)
        ax.set_yticks(range(len(d)), [t.replace("HALLMARK_", "").replace("CUSTOM_", "*")[:36] for t in d.Term], fontsize=7.5)
        for i, (nes, q) in enumerate(zip(d.NES, d["FDR q-val"])):
            ax.text(0.05 if nes < 0 else nes + 0.05, i, f"q={q:.2g}", va="center", ha="left", fontsize=6.5, color=INK2)
    ax.axvline(0, color=INK2, lw=0.6)
    ax.set_xlabel("NES")
    ax.set_title(f"{key[0].replace('_vs_', ' vs ').replace('_', ' ')}\n({MODELS[key[1]]})", fontsize=9)
    ax.grid(axis="y", visible=False)
fig.suptitle("GSEA within fibroblast-specific genes (FDR < 0.25); red = higher in first group", fontweight="bold")
fig.tight_layout()
save(fig, "F5_gsea_fibroblast.png")

# F6: heatmap of top fibroblast genes for each contrast (fibroblast-normalised VST)
dds_v = DeseqDataSet(counts=fc.T, metadata=md, design="~sex + group", quiet=True)
dds_v.fit_size_factors()
dds_v.vst(use_design=False)
vst = pd.DataFrame(dds_v.layers["vst_counts"], index=libs, columns=fc.index).T
vst.to_csv(os.path.join(OUT, "fibroblast_normalised_vst.csv"))
fig, axes = plt.subplots(1, 2, figsize=(16, 9))
for ax, cname in zip(axes, CONTRASTS):
    r = de[(cname, "fibnorm")].dropna(subset=["padj"])
    top = r.head(30).index
    a, b = CONTRASTS[cname]
    order = meta[meta.group.isin([b, a])].sort_values(["group", "sex"], key=lambda s: s.map({b: 0, a: 1}) if s.name == "group" else s).index
    sub = vst.loc[top, order]
    zz = sub.sub(sub.mean(1), axis=0).div(sub.std(1), axis=0)
    im = ax.imshow(zz.values, cmap=DIV, vmin=-2.5, vmax=2.5, aspect="auto")
    ax.set_yticks(range(len(top)), [f"{g2name[g]}  ({r.padj[g]:.2g})" for g in top], fontsize=7)
    ax.set_xticks(range(len(order)), [f"{meta.animal[l]} {meta.group[l].replace('_', ' ')} {meta.sex[l]}" for l in order], rotation=90, fontsize=7)
    for t, l in zip(ax.get_xticklabels(), order):
        t.set_color(GROUP_COLORS[meta.group[l]])
    ax.grid(False)
    ax.set_title(f"Top 30 fibroblast genes: {a.replace('_', ' ')} vs {b.replace('_', ' ')}\n(label = padj, ~sex + group)", fontsize=9.5)
fig.colorbar(im, ax=axes, shrink=0.5, label="row z-score (fibroblast-normalised VST)")
save(fig, "F6_heatmaps_fibroblast.png")

pd.concat([qc[["sample", "group", "sex", "pct_mito", "fibroblast_score", "fibroblast_fraction_deconv"]], props.add_prefix("deconv_")], axis=1).to_csv(
    os.path.join(OUT, "fibroblast_scores_per_sample.csv"))
print("Done:", OUT)

# ---------------------------------------------------------------- F7: SPP1 fibroblast ECM program vs Spp1 level
gs = gsea_out[("5xFAD_SPP1_vs_5xFAD_NTC", "fibnorm")]
emt = gs[gs.Term == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]
if len(emt):
    lead = set(emt.Lead_genes.iloc[0].split(";"))
    r = de[("5xFAD_SPP1_vs_5xFAD_NTC", "fibnorm")]
    ids = r.index[r.gene_name.astype(str).str.upper().isin(lead)]
    ecm = vst.loc[ids].sub(vst.loc[ids].mean(1), axis=0).mean()
    qc["fibroblast_ECM_leading_edge_score"] = ecm
    five = meta.index[meta.genotype == "5xFAD"]
    rho, p = spearmanr(ecm[five], np.log2(qc.Spp1_cpm[five] + 1))
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    ax = axes[0]
    for i, grp in enumerate(GROUP_ORDER):
        for j, sx in enumerate(["F", "M"]):
            idx = meta.index[(meta.group == grp) & (meta.sex == sx)]
            xpos = i + (-0.15 if sx == "F" else 0.15)
            ax.scatter(np.full(len(idx), xpos), ecm[idx], c=GROUP_COLORS[grp], marker=SEX_MARKER[sx], s=60,
                       edgecolor=SURFACE, lw=1.5, zorder=3)
            ax.hlines(ecm[idx].median(), xpos - 0.1, xpos + 0.1, color=INK, lw=1.5)
    ax.set_xticks(range(3), [g.replace("_", " ") for g in GROUP_ORDER])
    ax.set_ylabel("Fibroblast ECM score (EMT leading-edge genes,\nfibroblast-normalised VST, centred)")
    ax.set_title(f"Fibroblast ECM program is lower with SPP1 in both sexes\n({len(ids)} genes: {', '.join(sorted(g2name[ids].astype(str))[:9])} ...)", fontsize=9)
    ax.grid(axis="x", visible=False)
    legend(ax, loc="upper left")
    ax = axes[1]
    for grp in ["5xFAD_NTC", "5xFAD_SPP1"]:
        for sx, mk_ in SEX_MARKER.items():
            idx = meta.index[(meta.group == grp) & (meta.sex == sx)]
            ax.scatter(np.log2(qc.Spp1_cpm[idx] + 1), ecm[idx], c=GROUP_COLORS[grp], marker=mk_, s=70, edgecolor=SURFACE, lw=1.5, zorder=3)
    for l in five:
        ax.annotate(meta.animal[l], (np.log2(qc.Spp1_cpm[l] + 1), ecm[l]), xytext=(5, 4), textcoords="offset points", fontsize=7, color=INK2)
    ax.set_xlabel("Bulk Spp1 log2(CPM+1)")
    ax.set_ylabel("Fibroblast ECM score")
    ax.set_title(f"Within 5xFAD, fibroblast ECM score tracks Spp1 (Spearman ρ = {rho:.2f}, p = {p:.3f})", fontsize=9.5)
    fig.tight_layout()
    save(fig, "F7_spp1_fibroblast_ecm.png")
    pd.Series({"n_genes": len(ids), "genes": ";".join(g2name[ids].astype(str)), "spearman_rho_vs_Spp1_5xFAD": rho, "p": p}).to_csv(
        os.path.join(OUT, "ecm_score_vs_spp1.csv"))
    qc[["sample", "group", "sex", "Spp1_cpm", "fibroblast_ECM_leading_edge_score"]].to_csv(os.path.join(OUT, "ecm_score_per_sample.csv"))

# ---------------------------------------------------------------- ECM score robustness to contamination (5xFAD only)
if len(emt):
    import statsmodels.formula.api as smf
    dfm = pd.DataFrame({"ecm": ecm, "spp1": (meta.guide == "SPP1").astype(int), "sex": meta.sex,
                        "olfactory": qc.mod_olfactory_mucosa, "osteoblast": props["osteoblast"] if "osteoblast" in props else 0.0,
                        "osteoclast": props["osteoclast"] if "osteoclast" in props else 0.0}).loc[five]
    rows = []
    for label, formula, data in [
        ("SPP1 + sex", "ecm ~ spp1 + sex", dfm),
        ("+ olfactory contamination", "ecm ~ spp1 + sex + olfactory", dfm),
        ("+ deconvolved osteoblast", "ecm ~ spp1 + sex + osteoblast", dfm),
        ("+ deconvolved osteoclast", "ecm ~ spp1 + sex + osteoclast", dfm),
        ("+ olfactory + osteoblast", "ecm ~ spp1 + sex + olfactory + osteoblast", dfm),
        ("SPP1 + sex, excluding outlier 620", "ecm ~ spp1 + sex", dfm.drop(meta.index[meta.animal == "620"], errors="ignore")),
    ]:
        fit = smf.ols(formula, data).fit()
        rows.append({"model": label, "SPP1_effect": fit.params["spp1"], "p": fit.pvalues["spp1"], "n": int(fit.nobs)})
    rob = pd.DataFrame(rows)
    rob.to_csv(os.path.join(OUT, "ecm_score_spp1_robustness.csv"), index=False)
    print(rob.round(3).to_string())
