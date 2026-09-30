"""Bulk RNA-seq analysis of run VT3YXB (5xFAD vs WT, Spp1-targeted vs non-targeting control).

Usage: python3 analysis/run_analysis.py
Inputs are read from data/, all tables and figures are written to results/.
"""
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import spearmanr, mannwhitneyu
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import gseapy

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "results")
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)

# ---------------------------------------------------------------- styling
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
GROUP_COLORS = {"WT_NTC": "#2a78d6", "5xFAD_NTC": "#eb6834", "5xFAD_SPP1": "#1baf7a"}
GROUP_ORDER = list(GROUP_COLORS)
SEX_MARKER = {"F": "o", "M": "^"}
DIVERGING = LinearSegmentedColormap.from_list(
    "bgr", ["#184f95", "#3987e5", "#86b6ef", "#f0efec", "#f0a3a2", "#e34948", "#a8322f"])
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold", "legend.frameon": False,
})


def savefig(fig, name):
    fig.savefig(os.path.join(FIG, name), dpi=160, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- load
raw = pd.read_csv(os.path.join(DATA, "VT3YXB-expression-matrix.tsv"), sep="\t")
raw["gene_name"] = raw["gene_name"].fillna(raw["gene_id"])
raw = raw[~raw.gene_id.str.startswith("ERCC")]  # spike-ins carry no reads here

# Sample sheet: library index -> animal/label (from the facility's sample names)
labels = {
    1: "0453_5xFAD_M_SPP1", 2: "612_5xFAD_M_SPP1", 3: "613_5xFAD_M_SPP1",
    4: "621_5xFAD_F_SPP1", 5: "0463_5xFAD_F_SPP1", 6: "0455_5xFAD_F_SPP1",
    7: "616_WT_M_NTC", 8: "611_WT_M_NTC", 9: "626_WT_M_NTC",
    10: "624_WT_F_NTC", 11: "623_WT_F_NTC", 12: "0466_WT_F_NTC",
    13: "617_5xFAD_M_NTC", 14: "0452_5xFAD_M_NTC", 15: "0451_5xFAD_M_NTC",
    16: "629_5xFAD_F_NTC", 17: "0469_5xFAD_F_NTC", 18: "620_5xFAD_F_NTC",
}
libs = [f"VT3YXB_{i}" for i in labels]
meta = pd.DataFrame({"library": libs, "sample": [labels[i] for i in labels]}).set_index("library")
parts = meta["sample"].str.split("_", expand=True)
meta["animal"], meta["genotype"], meta["sex"], meta["guide"] = parts[0], parts[1], parts[2], parts[3]
meta["group"] = meta.genotype + "_" + meta.guide
meta["cohort"] = np.where(meta.animal.str.startswith("0"), "0xxx", "6xx")

counts = raw.set_index("gene_id")[[f"{l}_count" for l in libs]].astype(int)
counts.columns = libs
cpm_raw = raw.set_index("gene_id")[[f"{l}_cpm" for l in libs]]
cpm_raw.columns = libs
g2name = raw.set_index("gene_id")["gene_name"]
g2type = raw.set_index("gene_id")["gene_biotype"]
cpm_sym = cpm_raw.groupby(g2name).sum()
log_sym = np.log2(cpm_sym + 1)

# ---------------------------------------------------------------- QC
mapping = pd.read_csv(os.path.join(DATA, "VT3YXB-mapping-stats.csv"), index_col=0)
qc = meta.copy()
m = mapping.reindex(meta["sample"].values)
qc["total_reads"] = m.sum(axis=1).values
qc["uniquely_mapped"] = m["Uniquely Mapped"].values
qc["pct_unique"] = 100 * m["Uniquely Mapped"].values / qc["total_reads"]
qc["pct_multi"] = 100 * m["Multi-mapped"].values / qc["total_reads"]
qc["assigned_counts"] = counts.sum().values
qc["pct_assigned_of_unique"] = 100 * qc.assigned_counts / qc.uniquely_mapped
qc["pct_mito"] = 100 * cpm_sym[cpm_sym.index.str.startswith("mt-")].sum() / 1e6
qc["pct_hemoglobin"] = 100 * cpm_sym.loc[["Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2"]].sum() / 1e6
qc["genes_detected"] = (counts >= 5).sum().values
qc["Xist_cpm"] = cpm_sym.loc["Xist"]
qc["chrY_cpm"] = cpm_sym.loc[["Ddx3y", "Uty", "Kdm5d", "Eif2s3y"]].sum()
qc["sex_check"] = np.where(qc.Xist_cpm > qc.chrY_cpm, "F", "M")
qc["sex_ok"] = qc.sex_check == qc.sex

# Tissue composition modules: mean log2 CPM of marker genes, then z-scored across samples
MODULES = {
    "osteoclast": ["Ctsk", "Acp5", "Mmp9", "Atp6v0d2", "Dcstamp", "Ocstamp", "Calcr", "Oscar"],
    "osteoblast_bone": ["Bglap", "Bglap2", "Ibsp", "Sp7", "Dmp1", "Phex", "Mepe", "Sost"],
    "olfactory_mucosa": ["Omp", "Cyp2g1", "Ugt2a1", "Cyp2a5", "Gnal", "Scgb1c1"],
    "choroid_plexus": ["Ttr", "Kl", "Folr1", "Kcnj13", "Clic6", "Sostdc1", "F5"],
    "meningeal_fibroblast": ["Ptgds", "Slc47a1", "Crabp2", "Dcn", "Col1a1", "Col1a2"],
    "erythroid_blood": ["Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Alas2"],
    "neutrophil_marrow": ["S100a8", "S100a9", "Ngp", "Camp", "Ltf", "Chil3"],
    "brain_parenchyma": ["Snap25", "Syt1", "Stmn2", "Plp1", "Mobp"],
    "skeletal_muscle": ["Acta1", "Ckm", "Myh1", "Tnnt3", "Mylpf"],
    "border_macrophage": ["Mrc1", "Lyve1", "Cd163", "Pf4", "F13a1"],
    "microglia_homeostatic": ["P2ry12", "Tmem119", "Sall1", "Fcrls", "Gpr34"],
    "DAM_microglia": ["Cst7", "Itgax", "Clec7a", "Lpl", "Spp1", "Gpnmb", "Cd9", "Apoe"],
    "B_plasma_cell": ["Cd79a", "Cd79b", "Ms4a1", "Jchain", "Igkc", "Mzb1"],
}
modscore = pd.DataFrame({k: log_sym.reindex(v).dropna().mean() for k, v in MODULES.items()})
qc = qc.join(modscore.add_prefix("mod_"))
qc["Spp1_cpm"] = cpm_sym.loc["Spp1"]
qc.to_csv(os.path.join(OUT, "qc_and_composition.csv"))

# ---------------------------------------------------------------- gene filtering + VST
keep = (counts >= 10).sum(axis=1) >= 3
keep &= ~g2type.reindex(counts.index).isin(["rRNA", "Mt_rRNA", "Mt_tRNA", "misc_RNA", "snRNA", "snoRNA", "scaRNA", "miRNA", "ribozyme"])
cf = counts[keep]
print(f"Genes kept for modelling: {keep.sum()} / {len(keep)}")

design_meta = meta[["sex", "group", "cohort"]].copy()
design_meta["group"] = pd.Categorical(design_meta.group, categories=GROUP_ORDER)
dds0 = DeseqDataSet(counts=cf.T, metadata=design_meta, design="~sex + group", quiet=True)
dds0.fit_size_factors()
dds0.vst(use_design=False)
vst = pd.DataFrame(dds0.layers["vst_counts"], index=libs, columns=cf.index).T
vst_sym = vst.copy()
vst_sym.index = g2name.reindex(vst.index).values
qc["size_factor"] = dds0.obs["size_factors"].values if "size_factors" in dds0.obs else dds0.obsm.get("size_factors")

# ---------------------------------------------------------------- PCA
top = vst.var(axis=1).sort_values(ascending=False).index[:2000]
X = vst.loc[top].T
X = X - X.mean()
U, S, Vt = np.linalg.svd(X.values, full_matrices=False)
var_exp = S**2 / (S**2).sum()
pcs = pd.DataFrame(U * S, index=libs, columns=[f"PC{i+1}" for i in range(len(S))]).iloc[:, :6]
loadings = pd.DataFrame(Vt[:6].T, index=g2name.reindex(top).values, columns=pcs.columns)
loadings.to_csv(os.path.join(OUT, "pca_loadings_top2000.csv"))
pcs.join(meta).to_csv(os.path.join(OUT, "pca_scores.csv"))

covars = ["pct_mito", "pct_hemoglobin", "pct_unique", "assigned_counts"] + [f"mod_{k}" for k in MODULES] + ["Spp1_cpm"]
rows = {}
for pc in pcs.columns[:5]:
    r = {c: spearmanr(pcs[pc], qc[c])[0] for c in covars}
    r["genotype(5xFAD)"] = spearmanr(pcs[pc], (meta.genotype == "5xFAD").astype(int))[0]
    r["guide(SPP1)"] = spearmanr(pcs[pc], (meta.guide == "SPP1").astype(int))[0]
    r["sex(M)"] = spearmanr(pcs[pc], (meta.sex == "M").astype(int))[0]
    r["cohort(0xxx)"] = spearmanr(pcs[pc], (meta.cohort == "0xxx").astype(int))[0]
    rows[pc] = r
pc_cor = pd.DataFrame(rows)
pc_cor.to_csv(os.path.join(OUT, "pc_covariate_spearman.csv"))


def scatter_groups(ax, x, y, annotate=True):
    for g in GROUP_ORDER:
        for sx, mk in SEX_MARKER.items():
            idx = meta.index[(meta.group == g) & (meta.sex == sx)]
            ax.scatter(x[idx], y[idx], c=GROUP_COLORS[g], marker=mk, s=70,
                       edgecolor=SURFACE, linewidth=1.5, zorder=3)
    if annotate:
        for l in libs:
            ax.annotate(meta.animal[l], (x[l], y[l]), xytext=(5, 4), textcoords="offset points",
                        fontsize=7, color=INK2)


def group_legend(ax, loc="best"):
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="s", ls="", color=GROUP_COLORS[g], ms=9, label=g.replace("_", " ")) for g in GROUP_ORDER]
    h += [Line2D([], [], marker=mk, ls="", color=INK2, ms=8, label={"F": "female", "M": "male"}[s]) for s, mk in SEX_MARKER.items()]
    ax.legend(handles=h, loc=loc, fontsize=8)


fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
for ax, (a, b) in zip(axes, [("PC1", "PC2"), ("PC3", "PC4")]):
    scatter_groups(ax, pcs[a], pcs[b])
    ia, ib = int(a[2:]) - 1, int(b[2:]) - 1
    ax.set_xlabel(f"{a} ({var_exp[ia]*100:.1f}% variance)")
    ax.set_ylabel(f"{b} ({var_exp[ib]*100:.1f}% variance)")
    top_a = loadings[a].abs().sort_values(ascending=False).index[:6]
    ax.set_title(f"{a} vs {b}\n{a} top genes: {', '.join(top_a)}", fontsize=9)
group_legend(axes[0])
fig.suptitle("PCA of 2,000 most variable genes (VST)", fontweight="bold")
savefig(fig, "02_pca.png")

# PC-covariate correlation heatmap
fig, ax = plt.subplots(figsize=(6.5, 8))
im = ax.imshow(pc_cor.values, cmap=DIVERGING, vmin=-1, vmax=1, aspect="auto")
ax.set_xticks(range(pc_cor.shape[1]), [f"{c}\n{var_exp[i]*100:.0f}%" for i, c in enumerate(pc_cor.columns)])
ax.set_yticks(range(pc_cor.shape[0]), pc_cor.index)
ax.grid(False)
for i in range(pc_cor.shape[0]):
    for j in range(pc_cor.shape[1]):
        v = pc_cor.values[i, j]
        ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7,
                color="white" if abs(v) > 0.6 else INK)
fig.colorbar(im, ax=ax, shrink=0.5, label="Spearman ρ")
ax.set_title("What drives each principal component?")
savefig(fig, "03_pc_covariate_correlation.png")

# ---------------------------------------------------------------- QC figure
fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
order_idx = meta.sort_values(["group", "sex"]).index
xt = [meta.animal[l] for l in order_idx]
colors = [GROUP_COLORS[meta.group[l]] for l in order_idx]
for ax, (col, title) in zip(axes.flat, [
    ("total_reads", "Total reads (M)"), ("pct_unique", "% uniquely mapped"),
    ("pct_mito", "% reads mitochondrial"), ("pct_hemoglobin", "% reads hemoglobin"),
    ("genes_detected", "Genes with ≥5 reads"), ("Spp1_cpm", "Spp1 CPM")]):
    vals = qc.loc[order_idx, col] / (1e6 if col == "total_reads" else 1)
    ax.bar(range(len(vals)), vals, color=colors, width=0.8, edgecolor=SURFACE, linewidth=2)
    ax.set_xticks(range(len(vals)), xt, rotation=90, fontsize=7)
    ax.set_title(title)
    ax.grid(axis="x", visible=False)
group_legend(axes[0, 0], loc="lower right")
fig.suptitle("Library QC (bars coloured by group, sorted by group then sex)", fontweight="bold")
fig.tight_layout()
savefig(fig, "01_qc_overview.png")

# Composition heatmap
comp = qc.loc[order_idx, [f"mod_{k}" for k in MODULES]]
compz = (comp - comp.mean()) / comp.std()
fig, ax = plt.subplots(figsize=(11, 5.5))
im = ax.imshow(compz.T.values, cmap=DIVERGING, vmin=-2.5, vmax=2.5, aspect="auto")
ax.set_yticks(range(len(MODULES)), [k.replace("_", " ") for k in MODULES])
ax.set_xticks(range(len(order_idx)), [f"{meta.animal[l]} {meta.group[l].replace('_', ' ')} {meta.sex[l]}" for l in order_idx], rotation=90, fontsize=7)
for t, l in zip(ax.get_xticklabels(), order_idx):
    t.set_color(GROUP_COLORS[meta.group[l]])
ax.grid(False)
fig.colorbar(im, ax=ax, shrink=0.7, label="z-score (module mean log2 CPM)")
ax.set_title("Tissue / cell-type composition scores per sample")
savefig(fig, "04_composition_modules.png")

# Spp1 vs osteoclast module
fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
ax = axes[0]
scatter_groups(ax, qc.mod_osteoclast, np.log2(qc.Spp1_cpm + 1))
r, p = spearmanr(qc.mod_osteoclast, qc.Spp1_cpm)
ax.set_xlabel("Osteoclast module score (mean log2 CPM)")
ax.set_ylabel("Spp1 log2(CPM+1)")
ax.set_title(f"Spp1 tracks osteoclast content (Spearman ρ={r:.2f}, p={p:.3f})")
group_legend(ax, loc="upper left")
ax = axes[1]
for i, g in enumerate(GROUP_ORDER):
    idx = meta.index[meta.group == g]
    y = np.log2(qc.Spp1_cpm[idx] + 1)
    ax.scatter(np.full(len(idx), i) + np.linspace(-0.12, 0.12, len(idx)), y, c=GROUP_COLORS[g],
               s=60, edgecolor=SURFACE, linewidth=1.5, zorder=3)
    ax.hlines(np.median(y), i - 0.25, i + 0.25, color=INK, lw=2)
ax.set_xticks(range(3), [g.replace("_", " ") for g in GROUP_ORDER])
ax.set_ylabel("Spp1 log2(CPM+1)")
ax.set_title("Spp1 by group (line = median)")
ax.grid(axis="x", visible=False)
savefig(fig, "05_spp1_vs_osteoclast.png")

# ---------------------------------------------------------------- differential expression
def run_deseq(design, extra_cols, contrasts, tag):
    md = design_meta.copy()
    for c in extra_cols:
        v = qc[c]
        md[c] = ((v - v.mean()) / v.std()).values
    dds = DeseqDataSet(counts=cf.T, metadata=md, design=design, refit_cooks=True, quiet=True)
    dds.deseq2()
    out = {}
    for name, (a, b) in contrasts.items():
        st = DeseqStats(dds, contrast=["group", a, b], quiet=True)
        st.summary()
        res = st.results_df.copy()
        res.insert(0, "gene_name", g2name.reindex(res.index).values)
        res.insert(1, "biotype", g2type.reindex(res.index).values)
        res = res.sort_values("pvalue")
        res.to_csv(os.path.join(OUT, f"DE_{name}_{tag}.csv"))
        out[name] = res
    return out


CONTRASTS = {
    "5xFAD_NTC_vs_WT_NTC": ("5xFAD_NTC", "WT_NTC"),
    "5xFAD_SPP1_vs_5xFAD_NTC": ("5xFAD_SPP1", "5xFAD_NTC"),
    "5xFAD_SPP1_vs_WT_NTC": ("5xFAD_SPP1", "WT_NTC"),
}
print("Running DESeq2: base model ~ sex + group")
de_base = run_deseq("~sex + group", [], CONTRASTS, "base")
cov_cols = ["mod_osteoclast", "mod_choroid_plexus", "mod_olfactory_mucosa", "pct_mito"]
print("Running DESeq2: adjusted model ~ sex + composition covariates + group")
de_adj = run_deseq("~sex + mod_osteoclast + mod_choroid_plexus + mod_olfactory_mucosa + pct_mito + group",
                   cov_cols, CONTRASTS, "adjusted")

summary_rows = []
for tag, de in [("base", de_base), ("adjusted", de_adj)]:
    for name, res in de.items():
        sig = res[res.padj < 0.05]
        summary_rows.append({
            "model": tag, "contrast": name,
            "padj<0.05": len(sig), "up": (sig.log2FoldChange > 0).sum(), "down": (sig.log2FoldChange < 0).sum(),
            "padj<0.05 & |LFC|>1": (sig.log2FoldChange.abs() > 1).sum(),
            "padj<0.1": (res.padj < 0.1).sum(),
            "Spp1_LFC": res.loc[res.gene_name == "Spp1", "log2FoldChange"].squeeze(),
            "Spp1_padj": res.loc[res.gene_name == "Spp1", "padj"].squeeze(),
        })
de_summary = pd.DataFrame(summary_rows)
de_summary.to_csv(os.path.join(OUT, "DE_summary.csv"), index=False)
print(de_summary.to_string())


def volcano(ax, res, title, n_label=14):
    r = res.dropna(subset=["padj"])
    y = -np.log10(r.pvalue.clip(lower=1e-300))
    sig = r.padj < 0.05
    ax.scatter(r.log2FoldChange[~sig], y[~sig], s=6, c="#c9c8c3", lw=0, rasterized=True)
    up, dn = sig & (r.log2FoldChange > 0), sig & (r.log2FoldChange < 0)
    ax.scatter(r.log2FoldChange[up], y[up], s=12, c="#e34948", lw=0, label=f"up (n={up.sum()})")
    ax.scatter(r.log2FoldChange[dn], y[dn], s=12, c="#2a78d6", lw=0, label=f"down (n={dn.sum()})")
    lab = r[sig].head(n_label)
    for gid, row in lab.iterrows():
        ax.annotate(row.gene_name, (row.log2FoldChange, -np.log10(max(row.pvalue, 1e-300))),
                    xytext=(3, 2), textcoords="offset points", fontsize=7, color=INK)
    ax.axvline(0, color=INK2, lw=0.6)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 p-value")
    ax.set_title(title, fontsize=10)
    ax.legend(loc="upper left", fontsize=8)
    lim = np.nanpercentile(np.abs(r.log2FoldChange), 99.9)
    ax.set_xlim(-max(lim, 2), max(lim, 2))


for tag, de in [("base", de_base), ("adjusted", de_adj)]:
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
    for ax, (name, res) in zip(axes, de.items()):
        volcano(ax, res, name.replace("_vs_", " vs ").replace("_", " "))
    fig.suptitle(f"Differential expression ({'~ sex + group' if tag == 'base' else '~ sex + osteoclast + choroid plexus + olfactory + %mito + group'}), red/blue = padj<0.05",
                 fontweight="bold")
    fig.tight_layout()
    savefig(fig, f"06_volcano_{tag}.png")

# Heatmap of top genes in the genotype contrast
res = de_base["5xFAD_NTC_vs_WT_NTC"]
top_ids = res[res.padj < 0.05].head(50).index
if len(top_ids):
    sub = vst.loc[top_ids, order_idx]
    z = sub.sub(sub.mean(1), axis=0).div(sub.std(1), axis=0)
    fig, ax = plt.subplots(figsize=(10, max(4, len(top_ids) * 0.2)))
    im = ax.imshow(z.values, cmap=DIVERGING, vmin=-2.5, vmax=2.5, aspect="auto")
    ax.set_yticks(range(len(top_ids)), g2name.reindex(top_ids).values, fontsize=7)
    ax.set_xticks(range(len(order_idx)), [f"{meta.animal[l]} {meta.group[l].replace('_', ' ')} {meta.sex[l]}" for l in order_idx], rotation=90, fontsize=7)
    for t, l in zip(ax.get_xticklabels(), order_idx):
        t.set_color(GROUP_COLORS[meta.group[l]])
    ax.grid(False)
    fig.colorbar(im, ax=ax, shrink=0.5, label="row z-score (VST)")
    ax.set_title("Top 50 DE genes, 5xFAD NTC vs WT NTC (base model), all samples shown")
    savefig(fig, "07_heatmap_top_genotype_genes.png")

# ---------------------------------------------------------------- GSEA (preranked on Wald statistic)
def read_gmt(path):
    sets = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sets[f[0]] = [g for g in f[2:] if g]
    return sets


gene_sets = read_gmt(os.path.join(DATA, "genesets", "h.all.v7.0.symbols.gmt"))  # human symbols
# Custom mouse sets (upper-cased to match the human-symbol space used for ranking)
custom = {
    "CUSTOM_DAM_MICROGLIA_KerenShaul2017": ["Cst7", "Itgax", "Clec7a", "Lpl", "Spp1", "Gpnmb", "Cd9", "Apoe", "Trem2", "Tyrobp", "Axl", "Ank", "Csf1", "Ctsb", "Ctsd", "Ctsl", "Cd63", "Lgals3", "Igf1", "Ch25h", "Fabp5", "Lilrb4a", "Cd68"],
    "CUSTOM_HOMEOSTATIC_MICROGLIA": ["P2ry12", "P2ry13", "Tmem119", "Cx3cr1", "Sall1", "Fcrls", "Gpr34", "Olfml3", "Siglech", "Hexb", "Csf1r", "Selplg", "Serinc3"],
    "CUSTOM_OSTEOCLAST": MODULES["osteoclast"] + ["Nfatc1", "Car2", "Ctsk", "Itgb3", "Tnfrsf11a", "Atp6v1c1", "Snx10", "Clcn7", "Ostm1"],
    "CUSTOM_OSTEOBLAST": MODULES["osteoblast_bone"] + ["Runx2", "Alpl", "Col1a1", "Spp1", "Satb2", "Ifitm5"],
    "CUSTOM_CHOROID_PLEXUS": MODULES["choroid_plexus"] + ["Otx2", "Aqp1", "Slc4a5", "Enpp2", "Igfbp2", "Prlr"],
    "CUSTOM_INTERFERON_ISG": ["Ifit1", "Ifit3", "Isg15", "Irf7", "Oasl2", "Rsad2", "Usp18", "Stat1", "Bst2", "Ifi27l2a", "Mx1", "Oas1a", "Rtp4", "Ifitm3", "Xaf1"],
    "CUSTOM_MHCII_ANTIGEN": ["H2-Aa", "H2-Ab1", "H2-Eb1", "Cd74", "H2-DMa", "H2-DMb1", "Ciita"],
    "CUSTOM_PLASMA_B_CELL": MODULES["B_plasma_cell"] + ["Igha", "Ighm", "Xbp1", "Prdm1"],
}
gene_sets.update({k: [g.upper() for g in v] for k, v in custom.items()})


def gsea(res, name, tag):
    r = res.dropna(subset=["stat"]).copy()
    r["sym"] = r.gene_name.astype(str).str.upper()
    rnk = r.groupby("sym")["stat"].apply(lambda s: s.loc[s.abs().idxmax()]).sort_values(ascending=False)
    pre = gseapy.prerank(rnk=rnk, gene_sets=gene_sets, min_size=10, max_size=500,
                         permutation_num=2000, seed=7, threads=4, outdir=None, verbose=False)
    df = pre.res2d.copy()
    df["NES"] = df["NES"].astype(float)
    df["FDR q-val"] = df["FDR q-val"].astype(float)
    df = df.sort_values("NES", ascending=False)
    df.to_csv(os.path.join(OUT, f"GSEA_{name}_{tag}.csv"), index=False)
    return df


gsea_res = {}
for tag, de in [("base", de_base), ("adjusted", de_adj)]:
    for name in ["5xFAD_NTC_vs_WT_NTC", "5xFAD_SPP1_vs_5xFAD_NTC"]:
        print(f"GSEA {name} {tag}")
        gsea_res[(name, tag)] = gsea(de[name], name, tag)

fig, axes = plt.subplots(1, 4, figsize=(22, 7.5))
for ax, ((name, tag), df) in zip(axes, gsea_res.items()):
    d = df[df["FDR q-val"] < 0.25]
    d = pd.concat([d.head(12), d.tail(12)]).drop_duplicates("Term").sort_values("NES")
    if d.empty:
        ax.text(0.5, 0.5, "No gene sets at FDR < 0.25", ha="center", transform=ax.transAxes)
    else:
        cols = ["#e34948" if v > 0 else "#2a78d6" for v in d.NES]
        ax.barh(range(len(d)), d.NES, color=cols, height=0.75, edgecolor=SURFACE, lw=1.5)
        ax.set_yticks(range(len(d)), [t.replace("HALLMARK_", "").replace("CUSTOM_", "*")[:38] for t in d.Term], fontsize=7.5)
        for i, (nes, q) in enumerate(zip(d.NES, d["FDR q-val"])):
            ax.text(0.05 if nes < 0 else nes + 0.05, i, f"q={q:.2g}", va="center",
                    ha="left", fontsize=6.5, color=INK2)
    ax.axvline(0, color=INK2, lw=0.6)
    ax.set_xlabel("NES")
    ax.set_title(f"{name.replace('_vs_', ' vs ').replace('_', ' ')}\n({tag} model)", fontsize=9.5)
    ax.grid(axis="y", visible=False)
fig.suptitle("GSEA (Hallmark + *custom sets), FDR < 0.25; red = higher in first group", fontweight="bold")
fig.tight_layout()
savefig(fig, "08_gsea.png")

# Marker panel for biology of interest
panel = ["Spp1", "Cst7", "Itgax", "Clec7a", "Lpl", "Gpnmb", "Trem2", "Apoe", "Cd68", "Lgals3",
         "P2ry12", "Tmem119", "Cx3cr1", "Mrc1", "Lyve1", "Cd74", "H2-Aa", "Isg15", "Ifit3", "Gfap",
         "Ctsk", "Acp5", "Ibsp", "Ttr", "Snap25"]
panel = [g for g in panel if g in log_sym.index]
fig, axes = plt.subplots(5, 5, figsize=(16, 14), sharex=True)
for ax, g in zip(axes.flat, panel):
    for i, grp in enumerate(GROUP_ORDER):
        for sx, mk in SEX_MARKER.items():
            idx = meta.index[(meta.group == grp) & (meta.sex == sx)]
            ax.scatter(np.full(len(idx), i) + (-0.12 if sx == "F" else 0.12), log_sym.loc[g, idx],
                       c=GROUP_COLORS[grp], marker=mk, s=40, edgecolor=SURFACE, lw=1, zorder=3)
        ax.hlines(log_sym.loc[g, meta.index[meta.group == grp]].median(), i - 0.3, i + 0.3, color=INK, lw=1.5)
    ax.set_title(g, fontsize=10)
    ax.grid(axis="x", visible=False)
for ax in axes[-1]:
    ax.set_xticks(range(3), ["WT\nNTC", "5xFAD\nNTC", "5xFAD\nSPP1"], fontsize=8)
axes[0, 0].set_ylabel("log2(CPM+1)")
group_legend(axes[0, 0], loc="lower left")
fig.suptitle("Marker genes by group (circle = female, triangle = male, line = median)", fontweight="bold")
fig.tight_layout()
savefig(fig, "09_marker_panel.png")

# Is the Spp1 reduction in SPP1 animals explained by osteoclast content?
five = meta.index[meta.genotype == "5xFAD"]
ratio = np.log2(qc.Spp1_cpm + 1) - qc.mod_osteoclast
stats = {
    "Spp1 log2CPM, SPP1 vs NTC (5xFAD) MWU p": mannwhitneyu(np.log2(qc.Spp1_cpm[five][meta.guide[five] == "SPP1"] + 1),
                                                          np.log2(qc.Spp1_cpm[five][meta.guide[five] == "NTC"] + 1)).pvalue,
    "Osteoclast module, SPP1 vs NTC (5xFAD) MWU p": mannwhitneyu(qc.mod_osteoclast[five][meta.guide[five] == "SPP1"],
                                                               qc.mod_osteoclast[five][meta.guide[five] == "NTC"]).pvalue,
    "Spp1 minus osteoclast score, SPP1 vs NTC (5xFAD) MWU p": mannwhitneyu(ratio[five][meta.guide[five] == "SPP1"],
                                                                         ratio[five][meta.guide[five] == "NTC"]).pvalue,
}
pd.Series(stats).to_csv(os.path.join(OUT, "spp1_stats.csv"))
print(pd.Series(stats))
qc.to_csv(os.path.join(OUT, "qc_and_composition.csv"))
print("Done. Outputs in", OUT)
