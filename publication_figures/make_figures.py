"""Publication-style (GraphPad-like) figures for the SPP1 vs NTC claims, dura (run VT3YXB) and sorted CD11b+/ACSA-2+
brain cells (run 55ZS67). One PNG per claim, plus a stats CSV per figure, written to publication_figures/png and
publication_figures/stats.

Statistics (stated on every figure and in FIGURE_LEGENDS.md):
- Gene panels: DESeq2 Wald test (dura: ~ sex + treatment within 5xFAD, n = 6 v 6; CD11b: ~ treatment, n = 3 v 3).
  Asterisks = Benjamini-Hochberg q across the genes in that figure (the family). Genome-wide DESeq2 padj (BH across
  all ~19-20k tested genes) and an exact label-permutation p are given in the stats CSV.
- Signature / program scores: linear model (dura: + sex) with exact within-sex permutation p; BH across every score
  of that kind tested, not only the ones plotted.
- Gene sets: limma camera with inter-gene correlation estimated per set (conservative); FDR = BH across all sets.
- WT is shown as a reference only. WT vs 5xFAD is not tested: genotype is confounded with library processing.
Asterisks: ns q >= 0.05, * < 0.05, ** < 0.01, *** < 0.001, **** < 0.0001.
"""
import itertools
import os
import re

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DURA = os.path.join(REPO, "composition_confounding_analysis")
SORT = os.path.join(REPO, "sorted_glia_55ZS67")
PNG, STAT = os.path.join(HERE, "png"), os.path.join(HERE, "stats")
for d in (PNG, STAT):
    os.makedirs(d, exist_ok=True)

# ---------------------------------------------------------------- style (GraphPad-like)
plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 10, "axes.linewidth": 1.2, "axes.edgecolor": "#0b0b0b",
    "axes.labelcolor": "#0b0b0b", "xtick.color": "#0b0b0b", "ytick.color": "#0b0b0b",
    "xtick.major.width": 1.2, "ytick.major.width": 1.2, "xtick.major.size": 4, "ytick.major.size": 4,
    "xtick.direction": "out", "ytick.direction": "out", "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.dpi": 300, "axes.titleweight": "bold",
    "axes.titlesize": 11, "legend.frameon": False,
})
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
COL = {"WT_NTC": "#2a78d6", "5xFAD_NTC": "#eb6834", "5xFAD_SPP1": "#1baf7a", "NTC": "#eb6834", "SPP1": "#1baf7a"}
LAB = {"WT_NTC": "WT\nNTC", "5xFAD_NTC": "5xFAD\nNTC", "5xFAD_SPP1": "5xFAD\nSPP1", "NTC": "NTC", "SPP1": "SPP1"}
UP, DOWN = "#e34948", "#2a78d6"


def stars(q):
    if not np.isfinite(q):
        return "n/a"
    return "****" if q < 1e-4 else "***" if q < 1e-3 else "**" if q < 0.01 else "*" if q < 0.05 else "ns"


def fmt_p(p):
    if not np.isfinite(p):
        return "n/a"
    return f"{p:.1e}" if p < 1e-3 else f"{p:.3f}" if p < 0.1 else f"{p:.2f}"


def bh(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def save(fig, name, note=None):
    if note:
        fig.text(0.01, 0.005, note, fontsize=7.2, color=INK2, ha="left", va="bottom", wrap=True)
    fig.savefig(os.path.join(PNG, name), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", name)


def dot_panel(ax, values, groups, order, sex=None, ylabel=None, title=None, tested=("NTC", "SPP1"), q=None, ptxt=None,
              ylim0=True):
    """Scatter of individual animals + mean +/- SEM, GraphPad style; bracket between the tested pair."""
    values, groups = np.asarray(values, float), np.asarray(groups)
    for i, g in enumerate(order):
        v = values[groups == g]
        sx = None if sex is None else np.asarray(sex)[groups == g]
        n = len(v)
        jit = np.linspace(-0.16, 0.16, n) if n > 1 else np.zeros(1)
        jit = jit[np.argsort(np.argsort(v))]  # spread by rank so points do not overlap
        for k in range(n):
            mk = "^" if (sx is not None and sx[k] == "M") else "o"
            ax.scatter(i + jit[k], v[k], s=34, marker=mk, color=COL[g], edgecolor="white", linewidth=0.8, zorder=3)
        m, se = v.mean(), v.std(ddof=1) / np.sqrt(n)
        ax.plot([i - 0.28, i + 0.28], [m, m], color=INK, lw=1.6, zorder=4)
        ax.errorbar(i, m, yerr=se, color=INK, lw=1.2, capsize=5, capthick=1.2, zorder=4)
    ax.set_xticks(range(len(order)), [LAB[g] for g in order])
    ax.set_xlim(-0.6, len(order) - 0.4)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, style="italic" if re.match(r"^[A-Z][a-z0-9]+$", title) else "normal")
    top = np.nanmax(values)
    lo = 0 if ylim0 and np.nanmin(values) >= 0 else np.nanmin(values) - 0.1 * (top - np.nanmin(values))
    span = top - lo
    if q is not None and tested is not None:
        a, b = [order.index(t) for t in tested]
        y = top + 0.08 * span
        ax.plot([a, a, b, b], [y - 0.025 * span, y, y, y - 0.025 * span], color=INK, lw=1.0)
        ax.text((a + b) / 2, y + 0.01 * span, stars(q), ha="center", va="bottom", fontsize=10 if stars(q) != "ns" else 8.5)
        if ptxt:
            ax.text((a + b) / 2, y + 0.115 * span, ptxt, ha="center", va="bottom", fontsize=6.8, color=INK2)
        ax.set_ylim(lo, top + 0.32 * span)
    else:
        ax.set_ylim(lo, top + 0.1 * span)


def perm_within(y, trt, strata, stat="diff"):
    """Exact two-sided label-permutation p (within strata). Returns p and number of relabelings."""
    y, trt, strata = np.asarray(y, float), np.asarray(trt, int), np.asarray(strata)
    blocks = [np.where(strata == s)[0] for s in np.unique(strata)]
    choices = [list(itertools.combinations(b, int(trt[b].sum()))) for b in blocks]
    obs = y[trt == 1].mean() - y[trt == 0].mean()
    null = []
    for combo in itertools.product(*choices):
        t = np.zeros(len(y), int)
        for c in combo:
            t[list(c)] = 1
        null.append(y[t == 1].mean() - y[t == 0].mean())
    null = np.array(null)
    return float(np.mean(np.abs(null) >= abs(obs) - 1e-12)), len(null)


def forest(ax, labels, est, lo, hi, colors=None, p=None, xlabel="log2 fold change SPP1 vs NTC (95% CI)"):
    yy = np.arange(len(labels))[::-1]
    ax.axvline(0, color=INK2, lw=1)
    for i, y in enumerate(yy):
        c = colors[i] if colors is not None else INK
        ax.plot([lo[i], hi[i]], [y, y], color=c, lw=1.6, solid_capstyle="round")
        ax.scatter(est[i], y, s=36, color=c, edgecolor="white", linewidth=0.8, zorder=3)
    ax.set_yticks(yy, labels)
    ax.set_xlabel(xlabel)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    return yy


# ---------------------------------------------------------------- load dura
dm = pd.read_csv(os.path.join(REPO, "data", "VT3YXB-expression-matrix.tsv"), sep="\t")
dm = dm.drop_duplicates("gene_name").set_index("gene_name")
dcpm = dm[[c for c in dm.columns if c.endswith("_cpm")]]
dcpm.columns = [c.replace("_cpm", "") for c in dcpm.columns]
meta = pd.read_csv(os.path.join(DURA, "tables", "sample_metadata_with_composition.tsv"), sep="\t", index_col=0)
dcpm = dcpm[meta.index]
GORD = ["WT_NTC", "5xFAD_NTC", "5xFAD_SPP1"]
F5 = meta.index[meta.genotype == "5xFAD"]
de_u = pd.read_csv(os.path.join(DURA, "tables", "DE_SPP1_vs_NTC_5xFADonly.tsv"), sep="\t").drop_duplicates("gene_name").set_index("gene_name")
de_f = pd.read_csv(os.path.join(DURA, "tables", "DE_SPP1_vs_NTC_fibroblast_adjusted.tsv"), sep="\t").drop_duplicates("gene_name").set_index("gene_name")
sig = pd.read_csv(os.path.join(DURA, "tables", "composition_signature_scores_all.tsv"), sep="\t", index_col=0)
trt5 = (meta.loc[F5, "treatment"] == "SPP1").astype(int).values
sex5 = meta.loc[F5, "sex"].values

# ---------------------------------------------------------------- load sorted
sm_ = pd.read_csv(os.path.join(SORT, "data", "55ZS67-expression-matrix.tsv.gz"), sep="\t")
sm_ = sm_[~sm_.gene_id.str.startswith("ERCC")].drop_duplicates("gene_name").set_index("gene_name")
scnt = sm_[[c for c in sm_.columns if c.endswith("_count")]]
scnt.columns = [c.replace("_count", "") for c in scnt.columns]
scpm = scnt / scnt.sum() * 1e6
ss = pd.read_csv(os.path.join(SORT, "tables", "sample_sheet.tsv"), sep="\t", index_col=0)
CD = ss.index[ss.fraction == "CD11b"]
cd_trt = np.where(ss.loc[CD, "treatment"] == "Spp1", "SPP1", "NTC")
de_c = pd.read_csv(os.path.join(SORT, "tables", "DE_CD11b_SPP1_vs_NTC_unadjusted.tsv"), sep="\t").drop_duplicates("gene_name").set_index("gene_name")
de_a = pd.read_csv(os.path.join(SORT, "06_DE_ACSA2", "DE_Spp1_vs_NTC_ACSA2.tsv"), sep="\t").drop_duplicates("gene").set_index("gene")
de_a = de_a.rename(columns={"log2FC": "log2FoldChange"})


def gene_grid(genes, fname, title, dataset, ncol=4, note_extra="", ylab="CPM"):
    """One dot panel per gene; BH across the genes in this figure."""
    rows = []
    for g in genes:
        if dataset == "dura":
            d = de_u.loc[g] if g in de_u.index else None
            pp, nrel = perm_within(np.log2(dcpm.loc[g, F5].values + 1), trt5, sex5)
        else:
            d = de_c.loc[g] if g in de_c.index else None
            pp, nrel = perm_within(np.log2(scpm.loc[g, CD].values + 1), (cd_trt == "SPP1").astype(int), np.zeros(len(CD)))
        rows.append({"gene": g, "log2FC": d.log2FoldChange if d is not None else np.nan,
                     "lfcSE": d.lfcSE if d is not None else np.nan,
                     "DESeq2_p": d.pvalue if d is not None else np.nan,
                     "genomewide_padj": d.padj if d is not None else np.nan,
                     "exact_perm_p": pp, "n_relabelings": nrel})
    st = pd.DataFrame(rows)
    st["q_BH_within_figure"] = bh(st.DESeq2_p)
    st.to_csv(os.path.join(STAT, fname.replace(".png", ".csv")), index=False)
    nrow = int(np.ceil(len(genes) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(2.35 * ncol, 3.0 * nrow), squeeze=False)
    for ax, (_, r) in zip(axs.flat, st.iterrows()):
        g = r.gene
        if dataset == "dura":
            dot_panel(ax, dcpm.loc[g].values, meta.group.values, GORD, sex=meta.sex.values, title=g,
                      tested=("5xFAD_NTC", "5xFAD_SPP1"), q=r.q_BH_within_figure,
                      ptxt=f"q={fmt_p(r.q_BH_within_figure)}; FC={2 ** r.log2FC:.2f}\ngenome-wide padj={fmt_p(r.genomewide_padj)}")
        else:
            dot_panel(ax, scpm.loc[g, CD].values, cd_trt, ["NTC", "SPP1"], title=g, q=r.q_BH_within_figure,
                      ptxt=f"q={fmt_p(r.q_BH_within_figure)}; FC={2 ** r.log2FC:.2f}\ngenome-wide padj={fmt_p(r.genomewide_padj)}")
        ax.set_ylabel(ylab if ax in axs[:, 0] else "")
    for ax in list(axs.flat)[len(genes):]:
        ax.axis("off")
    fig.suptitle(title, fontsize=12, fontweight="bold", y=1.0)
    if dataset == "dura":
        leg = [Line2D([], [], marker="o", ls="", color=INK2, label="female"), Line2D([], [], marker="^", ls="", color=INK2, label="male")]
        fig.legend(handles=leg, loc="upper right", ncol=2, fontsize=8.5, bbox_to_anchor=(1.0, 1.0))
        note = ("Dots = animals (n = 6/group); line = mean ± SEM. Test: DESeq2 Wald, 5xFAD SPP1 vs 5xFAD NTC, ~ sex + treatment. "
                "Asterisks = Benjamini-Hochberg q across the genes in this figure. WT shown for reference only (not tested: "
                "genotype confounded with library processing). Genome-wide padj and exact within-sex permutation p in stats CSV. " + note_extra)
    else:
        note = ("Sorted CD11b+ brain cells, 5xFAD males, n = 3/group; dots = animals; line = mean ± SEM. Test: DESeq2 Wald (~ treatment). "
                "Asterisks = Benjamini-Hochberg q across the genes in this figure. With n = 3 v 3 the smallest exact permutation p is 0.1. "
                "Genome-wide padj in stats CSV. " + note_extra)
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    save(fig, fname, note)
    return st


# ================================================================ 0. workflow schematics
def box(ax, x, y, w, h, text, fc="#f4f4f2", ec=INK, fs=8.6, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.015", fc=fc, ec=ec, lw=1.2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, fontweight="bold" if bold else "normal", wrap=True)


def arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=13, lw=1.3, color=INK))


def schematic_dura():
    fig, ax = plt.subplots(figsize=(11, 5.6))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.0, 0.97, "Dura bulk RNA-seq workflow (run VT3YXB)", fontsize=13, fontweight="bold", va="top")
    gx = [0.03, 0.36, 0.69]
    for x, g, t in zip(gx, GORD, ["WT + NTC", "5xFAD + NTC", "5xFAD + SPP1-targeting"]):
        box(ax, x, 0.72, 0.28, 0.13, f"{t}\nn = 6 (3 F / 3 M)", fc="white", ec=COL[g], fs=9.5, bold=True)
    for x in gx:
        arrow(ax, x + 0.14, 0.71, 0.5 if x == gx[1] else (0.3 if x < 0.3 else 0.7), 0.645)
    steps = [("Dura mater dissected\n(whole dura, one per animal)", 0.03),
             ("Bulk RNA extraction\n(tissue, no sorting)", 0.22),
             ("Library prep +\nsequencing (facility)\n17.7 M reads/sample\n(range 12-21 M)", 0.41),
             ("Alignment + gene counts\n(facility; mouse reference,\n78,317 genes)", 0.60),
             ("Count matrix +\nmapping statistics\n(18 libraries)", 0.79)]
    for t, x in steps:
        box(ax, x, 0.46, 0.155, 0.17, t, fs=8.0)
    for (_, x0), (_, x1) in zip(steps[:-1], steps[1:]):
        arrow(ax, x0 + 0.158, 0.545, x1 - 0.003, 0.545)
    box(ax, 0.03, 0.12, 0.93, 0.25,
        "Known limitations carried into every analysis\n"
        "• Genotype confounded with library processing (5xFAD libraries ~10% mito / 17% multimapped vs ~5% / 10% in 5 of 6 WT) → WT vs 5xFAD not interpreted\n"
        "• Dissection carry-over (bone/marrow, nasal mucosa, muscle, blood) varies by animal → composition modelled explicitly\n"
        "• SPP1 dura has less fibroblast RNA (signature −1.1 SD) → fibroblast-adjusted model run alongside\n"
        "• Mouse-only reference: human APP/PSEN1 transgene reads not countable → genotype taken from lab records; sex verified (Xist / chrY genes)",
        fc="#fbf3ec", ec="#eb6834", fs=8.2)
    save(fig, "00a_workflow_dura.png")


def schematic_sorted():
    fig, ax = plt.subplots(figsize=(11, 6.2))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.0, 0.97, "Sorted brain cell workflow (run 55ZS67)", fontsize=13, fontweight="bold", va="top")
    box(ax, 0.02, 0.76, 0.32, 0.12, "5xFAD + NTC L2-siRNA\nn = 3 males", fc="white", ec=COL["NTC"], fs=9.5, bold=True)
    box(ax, 0.66, 0.76, 0.32, 0.12, "5xFAD + Spp1 L2-siRNA\nn = 3 males", fc="white", ec=COL["SPP1"], fs=9.5, bold=True)
    ax.text(0.5, 0.82, "15 nmol, single ICV injection\nat 6 months", ha="center", va="center", fontsize=8.4)
    ax.text(0.5, 0.715, "2 months", ha="center", va="center", fontsize=8.4, fontweight="bold")
    arrow(ax, 0.18, 0.75, 0.40, 0.685); arrow(ax, 0.82, 0.75, 0.60, 0.685)
    box(ax, 0.30, 0.58, 0.40, 0.10, "Harvest at 8 months: whole brain → warm (enzymatic)\ndissociation → single-cell suspension", fs=8.4)
    arrow(ax, 0.5, 0.575, 0.5, 0.52)
    box(ax, 0.30, 0.42, 0.40, 0.10, "CD11b+ selection (MACS)", fs=9, bold=True)
    arrow(ax, 0.30, 0.47, 0.20, 0.38); arrow(ax, 0.70, 0.47, 0.80, 0.38)
    box(ax, 0.03, 0.24, 0.34, 0.14, "CD11b+ fraction\n(microglia; 96-99% microglial RNA)\n6 libraries", fc="#eef5fd", ec=COL["WT_NTC"], fs=8.6)
    box(ax, 0.60, 0.28, 0.37, 0.10, "CD11b-depleted flow-through\n→ ACSA-2+ selection (MACS)", fs=8.6)
    arrow(ax, 0.785, 0.275, 0.785, 0.235)
    box(ax, 0.60, 0.10, 0.37, 0.13, "ACSA-2+ fraction (\"astrocyte\")\n~30-40% astrocyte, ~55% oligodendrocyte,\n~6% ependymal RNA; 6 libraries", fc="#eef5fd", ec=COL["WT_NTC"], fs=8.2)
    box(ax, 0.03, 0.04, 0.50, 0.15,
        "Bulk RNA-seq of each fraction (libraries re-run after poor first RNA)\n"
        "~19.9 M reads/sample; facility alignment + counts\n"
        "Both fractions come from the same animal (replicate k = mouse k), 6 mice total\n"
        "Not checkable in RNA: 5xFAD transgenes (mouse-only reference)", fc="#fbf3ec", ec="#eb6834", fs=7.9)
    save(fig, "00b_workflow_sorted_cells.png")


def schematic_analysis():
    fig, ax = plt.subplots(figsize=(12, 6.4))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.0, 0.98, "Analysis workflow (both datasets)", fontsize=13, fontweight="bold", va="top")
    cols = [
        ("1  QC", ["Mapping, % mito,\n% multimapped,\ngenes detected", "Sex check (Xist / chrY)", "Test: do QC metrics\nseparate the groups?"]),
        ("2  Composition / purity", ["Validated marker signatures\n(item-rest r ≥ 0.4)", "Deconvolution vs CELLxGENE\nreference (trimmed NNLS)", "Dura: fibroblast, myeloid,\ncarry-over tissues\nSorts: purity + carry-over"]),
        ("3  Differential expression", ["DESeq2, SPP1 vs NTC\n(dura: ~ sex + treatment,\n5xFAD only)", "Adjusted models: dura +\nfibroblast score; sorts:\nmicroglia-intrinsic genes", "Exact permutation,\nleave-one-out"]),
        ("4  Pathways", ["limma camera on Hallmark,\ncurated + mechanistic sets", "Conservative variant\n(inter-gene correlation\nestimated) used for claims", "Link map: shared\nleading genes"]),
        ("5  Integration", ["Dura vs microglia\ngene + pathway comparison", "Composition checks on every\nhit (e.g. IFN vs monocytes)", "PANTHER / IPA upload kits\n(expressed-gene background)"]),
    ]
    w, gap = 0.172, 0.026
    for i, (head, items) in enumerate(cols):
        x = 0.01 + i * (w + gap)
        box(ax, x, 0.80, w, 0.09, head, fc="#0b0b0b", ec=INK, fs=9.5, bold=True)
        ax.texts[-1].set_color("white")
        for j, t in enumerate(items):
            box(ax, x, 0.60 - j * 0.20, w, 0.16, t, fs=8.2)
        if i < len(cols) - 1:
            arrow(ax, x + w + 0.002, 0.845, x + w + gap - 0.002, 0.845)
    ax.text(0.01, 0.02, "Multiple testing: genome-wide BH (DESeq2 padj); BH within each figure's gene family; BH across gene sets (camera FDR). "
            "Claims rest on effects that survive composition checks, not on single nominal p-values.", fontsize=8.2, color=INK2)
    save(fig, "00c_workflow_analysis.png")


schematic_dura(); schematic_sorted(); schematic_analysis()

# ================================================================ 1. Spp1 knockdown (both tissues)
st = []
fig, axs = plt.subplots(1, 3, figsize=(9.6, 3.9), gridspec_kw={"width_ratios": [1.5, 1, 1]})
d = de_u.loc["Spp1"]; c = de_c.loc["Spp1"]
pd_, _ = perm_within(np.log2(dcpm.loc["Spp1", F5].values + 1), trt5, sex5)
pc_, _ = perm_within(np.log2(scpm.loc["Spp1", CD].values + 1), (cd_trt == "SPP1").astype(int), np.zeros(len(CD)))
dot_panel(axs[0], dcpm.loc["Spp1"].values, meta.group.values, GORD, sex=meta.sex.values, ylabel="Spp1 (CPM)",
          title="Dura (bulk)", tested=("5xFAD_NTC", "5xFAD_SPP1"), q=d.pvalue, ptxt=f"p={fmt_p(d.pvalue)} (padj {fmt_p(d.padj)})\nFC={2 ** d.log2FoldChange:.2f}")
dot_panel(axs[1], scpm.loc["Spp1", CD].values, cd_trt, ["NTC", "SPP1"], ylabel="Spp1 (CPM)", title="CD11b+ microglia",
          q=c.padj, ptxt=f"padj={fmt_p(c.padj)}\nFC={2 ** c.log2FoldChange:.2f}")
ki = pd.read_csv(os.path.join(SORT, "09_pathology", "spp1_knockdown_index.tsv"), sep="\t", index_col=0)
kt = np.where(ki.treatment == "Spp1", "SPP1", "NTC")
wp = stats.ttest_ind(ki.Spp1_minus_DAM2[kt == "SPP1"], ki.Spp1_minus_DAM2[kt == "NTC"], equal_var=False).pvalue
dot_panel(axs[2], ki.Spp1_minus_DAM2.values, kt, ["NTC", "SPP1"], ylabel="log2 Spp1 − log2 DAM2 score", title="Spp1 per DAM response",
          q=wp, ptxt=f"Welch p={fmt_p(wp)}", ylim0=False)
fig.tight_layout(rect=(0, 0.1, 1, 1))
save(fig, "01_Spp1_knockdown_dura_and_microglia.png",
     "Dura: DESeq2 Wald within 5xFAD (~ sex + treatment), single pre-specified gene, so no within-figure correction; genome-wide padj shown "
     "(drop is male-driven). Microglia: asterisks use genome-wide DESeq2 padj (Spp1 = #1 gene of ~15,500). Right: Spp1 relative to the DAM2 "
     "program (Welch t-test), showing the drop is specific rather than a weaker plaque response. Exact permutation p: dura "
     f"{fmt_p(pd_)}, microglia {fmt_p(pc_)} (minimum possible 0.1 at n = 3 v 3).")
pd.DataFrame([{"panel": "dura", "log2FC": d.log2FoldChange, "lfcSE": d.lfcSE, "DESeq2_p": d.pvalue, "genomewide_padj": d.padj, "exact_perm_p": pd_},
              {"panel": "CD11b", "log2FC": c.log2FoldChange, "lfcSE": c.lfcSE, "DESeq2_p": c.pvalue, "genomewide_padj": c.padj, "exact_perm_p": pc_},
              {"panel": "Spp1_minus_DAM2", "welch_p": wp}]).to_csv(os.path.join(STAT, "01_Spp1_knockdown_dura_and_microglia.csv"), index=False)

# ================================================================ 2. dura ECM genes
ECM = ["Fn1", "Col1a1", "Col1a2", "Col3a1", "Col12a1", "Lox", "Loxl3", "Fbln1", "Timp1", "Cthrc1", "Prg4", "Matn4"]
ecm = gene_grid(ECM, "02_dura_ECM_genes.png", "Dura extracellular-matrix genes, SPP1 vs NTC (5xFAD)", "dura")

# ================================================================ 3. ECM attenuation by fibroblast adjustment
fig, ax = plt.subplots(figsize=(6.4, 6.2))
labs, est, lo, hi, cols = [], [], [], [], []
rows = []
for g in ECM:
    for tag, dd, colr in [("unadjusted", de_u, "#eb6834"), ("fibroblast-adjusted", de_f, "#52514e")]:
        r = dd.loc[g]
        rows.append({"gene": g, "model": tag, "log2FC": r.log2FoldChange, "lfcSE": r.lfcSE, "p": r.pvalue, "genomewide_padj": r.padj})
s3 = pd.DataFrame(rows)
s3["q_BH_within_figure_model"] = s3.groupby("model").p.transform(bh)
yy = np.arange(len(ECM))[::-1] * 1.0
ax.axvline(0, color=INK2, lw=1)
for k, (tag, colr, off) in enumerate([("unadjusted", "#eb6834", 0.17), ("fibroblast-adjusted", "#52514e", -0.17)]):
    sub = s3[s3.model == tag].set_index("gene").loc[ECM]
    for i, g in enumerate(ECM):
        r = sub.loc[g]; y = yy[i] + off
        ax.plot([r.log2FC - 1.96 * r.lfcSE, r.log2FC + 1.96 * r.lfcSE], [y, y], color=colr, lw=1.6)
        ax.scatter(r.log2FC, y, s=34, color=colr, edgecolor="white", lw=0.8, zorder=3, marker="o" if k == 0 else "s")
        ax.text(1.05, y, stars(r.q_BH_within_figure_model), fontsize=7.5, va="center", ha="left", color=INK2, transform=ax.get_yaxis_transform())
ax.set_yticks(yy, ECM, style="italic")
ax.set_xlabel("log2 fold change, 5xFAD SPP1 vs NTC (95% CI)")
ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
ax.legend(handles=[Line2D([], [], color="#eb6834", marker="o", label="~ sex + treatment"),
                   Line2D([], [], color="#52514e", marker="s", label="~ sex + fibroblast score + treatment")],
          loc="lower left", fontsize=8.5, bbox_to_anchor=(0, 1.0))
ax.set_title("ECM decreases shrink after adjusting for fibroblast content", pad=36)
fig.tight_layout(rect=(0, 0.08, 1, 1))
s3.to_csv(os.path.join(STAT, "03_dura_ECM_fibroblast_adjustment.csv"), index=False)
save(fig, "03_dura_ECM_fibroblast_adjustment.png",
     "DESeq2 log2FC ± 1.96 × SE, n = 6 v 6. Right-hand asterisks = BH q across these 12 genes within each model (upper = unadjusted, "
     "lower = adjusted). Caveat: Fbln1 is part of the fibroblast signature, so its adjusted value is partly circular.")

# ================================================================ 4. dura composition: fibroblast & monocyte signatures
rows = []
for s in sig.columns:
    y = sig.loc[F5, s].values
    X = sm.add_constant(np.column_stack([(sex5 == "M").astype(float), trt5]))
    fit = sm.OLS(y, X).fit()
    pp, _ = perm_within(y, trt5, sex5)
    rows.append({"signature": s, "SPP1_minus_NTC": fit.params[2], "ci_lo": fit.conf_int()[2][0], "ci_hi": fit.conf_int()[2][1], "ols_p": fit.pvalues[2], "exact_perm_p": pp})
s4 = pd.DataFrame(rows)
s4["q_BH_all_20_signatures"] = bh(s4.exact_perm_p)
s4.to_csv(os.path.join(STAT, "04_dura_composition_signatures.csv"), index=False)
fig, axs = plt.subplots(1, 3, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1, 1, 1.6]})
for ax, s, nm in [(axs[0], "fibroblast_stromal", "Fibroblast signature"), (axs[1], "monocyte", "Monocyte signature")]:
    r = s4.set_index("signature").loc[s]
    dot_panel(ax, sig.loc[meta.index, s].values, meta.group.values, GORD, sex=meta.sex.values, ylabel="score (mean z of markers)",
              title=nm, tested=("5xFAD_NTC", "5xFAD_SPP1"), q=r.q_BH_all_20_signatures,
              ptxt=f"perm p={fmt_p(r.exact_perm_p)}; q={fmt_p(r.q_BH_all_20_signatures)}", ylim0=False)
o = s4.sort_values("SPP1_minus_NTC")
ax = axs[2]
yy = np.arange(len(o))
ax.axvline(0, color=INK2, lw=1)
for i, (_, r) in enumerate(o.iterrows()):
    cc = DOWN if r.SPP1_minus_NTC < 0 else UP
    ax.plot([r.ci_lo, r.ci_hi], [i, i], color=cc, lw=1.4, alpha=0.9 if r.exact_perm_p < 0.05 else 0.45)
    ax.scatter(r.SPP1_minus_NTC, i, s=26, color=cc, edgecolor="white", lw=0.6, zorder=3, alpha=1 if r.exact_perm_p < 0.05 else 0.5)
ax.set_yticks(yy, [s.replace("_", " ") for s in o.signature], fontsize=7.5)
ax.set_xlabel("SPP1 − NTC (signature units, 95% CI)")
ax.set_title("All 20 composition signatures")
ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
fig.tight_layout(rect=(0, 0.1, 1, 1))
save(fig, "04_dura_composition_fibroblast_monocyte.png",
     "Signature score = mean z-score of validated markers (blind VST). SPP1 vs NTC within 5xFAD: OLS ~ sex + treatment for CI; exact within-sex "
     "permutation p (400 relabelings); BH across all 20 signatures tested (asterisks). Faded = perm p ≥ 0.05. The fibroblast decrease is nominal "
     "(perm p < 0.05) but does not survive correction across 20 signatures.")

# ================================================================ 5. dura robust single genes
rob = gene_grid(["Cldn11", "Prg4", "Rgs4", "Ifit3b", "Matn4"], "05_dura_robust_SPP1_genes.png",
                "Dura genes most robust to composition/technical models", "dura", ncol=5,
                note_extra="Cldn11, Prg4, Rgs4 are padj < 0.05 in 4-5 of 5 covariate models; Ifit3b and Matn4 strengthen after fibroblast adjustment.")

# ================================================================ 6. camera bar charts (dura, microglia)
ACR = {"uv": "UV", "dn": "down", "up": "up", "yap": "YAP", "taz": "TAZ", "ifn": "IFN", "isg": "ISG", "tgf": "TGF", "tgfb": "TGF-β",
       "smad": "SMAD", "dam": "DAM", "il6": "IL6", "jak": "JAK", "stat3": "STAT3", "mhcii": "MHC-II", "kras": "KRAS",
       "mtorc1": "mTORC1", "srebp": "SREBP", "ecm": "ECM", "arm": "ARM", "mgnd": "MGnD", "dna": "DNA", "i": "I", "a2": "A2",
       "kerenshaul2017": "(Keren-Shaul 2017)", "dam1": "DAM1", "dam2": "DAM2", "nospp1": "(no Spp1)", "a": "A", "m": "M"}


def nice(s):
    pre = "Astro: " if s.startswith("G_A_") else "Microglia: " if s.startswith("G_M_") else ""
    w = re.sub(r"^(H_|C_|M_|G_M_|G_A_|G_)", "", s).replace("_", " ").split()
    w = [ACR.get(t.lower(), t.lower()) for t in w]
    if w and w[0][0].islower():
        w[0] = w[0][0].upper() + w[0][1:]
    return pre + " ".join(w).replace("TGF beta", "TGF-β")


def camera_bars(cam, models, names, fname, title, n_top=14, note=""):
    fig, axs = plt.subplots(1, len(models), figsize=(6.2 * len(models), 6.2))
    axs = np.atleast_1d(axs)
    out = []
    for ax, m, nm in zip(axs, models, names):
        c = cam[(cam.model == m) & ~cam.set.str.startswith("TECHNICAL")].copy()
        c = c.sort_values("p_conservative").head(n_top).iloc[::-1]
        val = -np.log10(c.p_conservative) * np.where(c.direction == "Up", 1, -1)
        y = np.arange(len(c))
        ax.barh(y, val, height=0.62, color=[UP if v > 0 else DOWN for v in val], edgecolor="white", lw=0)
        ax.axvline(0, color=INK, lw=1)
        for v in (-np.log10(0.05), np.log10(0.05)):
            ax.axvline(v, color=MUTED, lw=0.9, ls=(0, (3, 2)))
        lab = [nice(s) for s in c.set]
        ax.set_yticks(y, lab, fontsize=8.2)
        ax.tick_params(axis="y", length=0)
        ax.set_xlabel("signed −log10 p (camera, conservative)\n← down with SPP1   |   up with SPP1 →")
        ax.set_title(nm)
        ax.spines["left"].set_visible(False)
        lim = max(2.2, np.abs(val).max() * 1.15)
        ax.set_xlim(-lim, lim)
        for yi, (_, r) in zip(y, c.iterrows()):
            ax.text(1.01, yi, f"FDR {r.FDR_conservative:.2f}", ha="left", va="center", fontsize=6.8, color=INK2,
                    transform=ax.get_yaxis_transform())
        out.append(c.assign(panel=nm))
    fig.suptitle(title, fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0.06, 0.97, 0.97), w_pad=6)
    pd.concat(out).to_csv(os.path.join(STAT, fname.replace(".png", ".csv")), index=False)
    save(fig, fname, "limma camera, competitive test, inter-gene correlation estimated per set (conservative). Dashed line = p 0.05. "
         "FDR = Benjamini-Hochberg across all sets tested in that model; no set reaches FDR < 0.05, so these are ranked hypotheses. " + note)


dcam = pd.read_csv(os.path.join(DURA, "tables", "spp1_pathway_camera_both_models.tsv"), sep="\t")
camera_bars(dcam, ["unadjusted", "fibro_adj"], ["unadjusted (~ sex + treatment)", "fibroblast-adjusted"],
            "06_dura_pathways_camera.png", "Dura gene-set changes, 5xFAD SPP1 vs NTC (n = 6 v 6)",
            note="Matrix / YAP-TAZ / TGF-β sets fade after fibroblast adjustment; interferon, eosinophil and IL6/JAK-STAT sets strengthen.")
ccam = pd.read_csv(os.path.join(SORT, "tables", "CD11b_pathway_camera_both_universes.tsv"), sep="\t")
camera_bars(ccam, ["all_genes", "mg_intrinsic"], ["all genes", "microglia-intrinsic genes only"],
            "10_microglia_pathways_camera.png", "CD11b+ microglia gene-set changes, SPP1 vs NTC (n = 3 v 3)",
            note="Up-regulated sets in 'all genes' (ECM, angiogenesis, EMT) are carried by glial/BAM/vascular carry-over and vanish when only "
                 "microglia-intrinsic genes are tested.")

# ================================================================ 7. dura IFN & eosinophil scores vs composition
sc = pd.read_csv(os.path.join(DURA, "13_pathway_links", "eosinophil_IFN_scores_per_animal.csv"), index_col=0)
rob_ = pd.read_csv(os.path.join(DURA, "13_pathway_links", "IFN_eosinophil_robustness.csv"))
std = lambda x: (x - x.mean()) / x.std()
male = (sc.sex == "M").astype(float).values
spp = (sc.group == "5xFAD_SPP1").astype(float).values
fig, axs = plt.subplots(1, 3, figsize=(11.5, 4.3), gridspec_kw={"width_ratios": [0.8, 0.8, 1.6]})
q7 = bh(rob_[rob_.covariates == "sex only"].perm_p.values)
for ax, col, nm, q in [(axs[0], "IFN_score", "Interferon score", q7[0]), (axs[1], "eosinophil_score", "Eosinophil granule score", q7[1])]:
    pp = rob_[(rob_.outcome == col) & (rob_.covariates == "sex only")].perm_p.iloc[0]
    dot_panel(ax, sc[col].values, sc.group.values, ["5xFAD_NTC", "5xFAD_SPP1"], sex=sc.sex.values, ylabel="score (mean z)",
              title=nm, tested=("5xFAD_NTC", "5xFAD_SPP1"), q=q, ptxt=f"perm p={fmt_p(pp)}; q={fmt_p(q)}", ylim0=False)
ax = axs[2]
rows = []
for col, colr, off in [("IFN_score", "#7a3fbf", 0.15), ("eosinophil_score", "#c2410c", -0.15)]:
    y = std(sc[col]).values
    for i, covs in enumerate([[], ["fibroblast_stromal"], ["monocyte"], ["neutrophil"], ["monocyte", "neutrophil"], ["fibroblast_stromal", "monocyte"]]):
        X = sm.add_constant(np.column_stack([male, spp] + [std(sig.loc[sc.index, cv]).values for cv in covs]))
        f = sm.OLS(y, X).fit()
        lab = "sex only" if not covs else "+ " + " + ".join(covs)
        pp = rob_[(rob_.outcome == col) & (rob_.covariates == lab)].perm_p.iloc[0]
        rows.append({"outcome": col, "covariates": lab, "effect_SD": f.params[2], "ci_lo": f.conf_int()[2][0], "ci_hi": f.conf_int()[2][1], "ols_p": f.pvalues[2], "exact_perm_p": pp})
        yv = 5 - i + off
        ax.plot(f.conf_int()[2], [yv, yv], color=colr, lw=1.6)
        ax.scatter(f.params[2], yv, s=34, color=colr, edgecolor="white", lw=0.8, zorder=3)
ax.axvline(0, color=INK2, lw=1)
ax.set_yticks(range(6)[::-1], ["sex only", "+ fibroblast", "+ monocyte", "+ neutrophil", "+ monocyte + neutrophil", "+ fibroblast + monocyte"], fontsize=8.5)
ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
ax.set_xlabel("SPP1 effect (SD units, 95% CI)")
ax.set_title("Effect after composition covariates")
ax.legend(handles=[Line2D([], [], color="#7a3fbf", marker="o", label="interferon"), Line2D([], [], color="#c2410c", marker="o", label="eosinophil")],
          fontsize=8.5, loc="lower right")
pd.DataFrame(rows).to_csv(os.path.join(STAT, "07_dura_IFN_eosinophil_vs_composition.csv"), index=False)
fig.tight_layout(rect=(0, 0.1, 1, 1))
save(fig, "07_dura_IFN_eosinophil_vs_composition.png",
     "Scores = mean z of module genes (blind VST), 12 5xFAD animals. Left: exact within-sex permutation p (400 relabelings), BH across the two scores. "
     "Right: OLS effect of SPP1 on the standardised score with sex plus composition covariates (95% CI); exact permutation p per model in stats CSV. "
     "The interferon effect largely disappears with monocyte adjustment.")

# ================================================================ 8. microglia DAM arms
dam = gene_grid(["Gpnmb", "Apoe", "Ch25h", "Fabp5", "Igf1", "Cst7", "Itgax", "Trem2", "Tyrobp", "Lpl", "P2ry12", "Tmem119", "Cx3cr1", "Sall1", "Csf1r"],
                "08_microglia_lipid_arm_vs_core_DAM.png",
                "Microglia: lipid-associated genes fall (row 1); core DAM (row 2) and homeostatic genes (row 3) do not", "cd11b", ncol=5,
                note_extra="Igf1 is BAM-enriched in the reference (possible carry-over); all others are microglia-intrinsic.")

# ================================================================ 9. microglia program scores (power)
pt = pd.read_csv(os.path.join(SORT, "09_pathology", "pathology_programs_SPP1_vs_NTC.tsv"), sep="\t")
pt = pt[pt.fraction == "CD11b"].copy()
pt["q_BH"] = bh(pt.welch_p)
pt.to_csv(os.path.join(STAT, "09_microglia_program_scores.csv"), index=False)
pt = pt.sort_values("pct_change")
fig, ax = plt.subplots(figsize=(6.6, 5.2))
yy = np.arange(len(pt))
ax.axvspan(-10, 10, color=GRID, alpha=0.6, lw=0)
ax.axvline(0, color=INK2, lw=1)
for i, (_, r) in enumerate(pt.iterrows()):
    ax.plot([r.ci95_lo_pct, r.ci95_hi_pct], [i, i], color=INK, lw=1.6)
    ax.scatter(r["pct_change"], i, s=40, color=COL["SPP1"], edgecolor="white", lw=0.8, zorder=3)
    ax.text(1.02, i, f"q={fmt_p(r.q_BH)}", transform=ax.get_yaxis_transform(), fontsize=7.5, va="center", color=INK2)
ax.set_yticks(yy, [s.replace("_", " ") for s in pt.set], fontsize=8.6)
ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
ax.set_xlabel("% change in program score, SPP1 vs NTC (95% CI)")
ax.set_title("Microglial disease programs are unchanged")
fig.tight_layout(rect=(0, 0.1, 0.93, 1))
save(fig, "09_microglia_program_scores.png",
     "Program score = mean log2 CPM of published gene sets (Spp1 removed; dissociation genes excluded), n = 3 v 3. Welch t-test; "
     "BH q across the 13 programs. Grey band = ±10%. The CI width shows the smallest change this design could detect (e.g. DAM2 −15 to +10%).")

# ================================================================ 10/11. cross-tissue forests (IFN, cholesterol)
def cross_forest(genes, fname, title, include_acsa=False, note=""):
    sets = [("dura, unadjusted", de_u, "#eb6834"), ("dura, fibroblast-adjusted", de_f, "#52514e"), ("CD11b+ microglia", de_c, "#1baf7a")]
    if include_acsa:
        sets.append(("ACSA-2+ fraction", de_a, "#7a3fbf"))
    rows = []
    for nm, dd, _ in sets:
        for g in genes:
            if g in dd.index:
                r = dd.loc[g]
                rows.append({"gene": g, "dataset": nm, "log2FC": r.log2FoldChange, "lfcSE": r.lfcSE, "p": r.pvalue, "genomewide_padj": r.padj})
            else:
                rows.append({"gene": g, "dataset": nm, "log2FC": np.nan, "lfcSE": np.nan, "p": np.nan, "genomewide_padj": np.nan})
    s = pd.DataFrame(rows)
    s["q_BH_within_figure_dataset"] = s.groupby("dataset").p.transform(bh)
    s.to_csv(os.path.join(STAT, fname.replace(".png", ".csv")), index=False)
    k = len(sets)
    fig, ax = plt.subplots(figsize=(7.2, 0.62 * len(genes) * (k / 3) + 1.8))
    base = np.arange(len(genes))[::-1] * (k * 0.24 + 0.35)
    ax.axvline(0, color=INK2, lw=1)
    for j, (nm, dd, colr) in enumerate(sets):
        sub = s[s.dataset == nm].set_index("gene")
        for i, g in enumerate(genes):
            r = sub.loc[g]
            y = base[i] + (k - 1) * 0.12 - j * 0.24
            if not np.isfinite(r.log2FC):
                ax.text(0, y, "not expressed", fontsize=6.5, color=MUTED, va="center", ha="center")
                continue
            ax.plot([r.log2FC - 1.96 * r.lfcSE, r.log2FC + 1.96 * r.lfcSE], [y, y], color=colr, lw=1.5)
            ax.scatter(r.log2FC, y, s=26, color=colr, edgecolor="white", lw=0.6, zorder=3)
            ax.text(1.01, y, stars(r.q_BH_within_figure_dataset), transform=ax.get_yaxis_transform(), fontsize=7, va="center", color=INK2)
    ax.set_yticks(base, genes, style="italic")
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    ax.set_xlabel("log2 fold change SPP1 vs NTC (95% CI)")
    ax.legend(handles=[Line2D([], [], color=cc, marker="o", label=nm) for nm, _, cc in sets], fontsize=8, loc="lower left",
              bbox_to_anchor=(0, 1.0), ncol=2)
    ax.set_title(title, pad=8 + 14 * int(np.ceil(k / 2)))
    fig.tight_layout(rect=(0, 0.09, 0.95, 1))
    save(fig, fname, "DESeq2 log2FC ± 1.96 × SE. Dura n = 6 v 6 (5xFAD, ~ sex + treatment); CD11b/ACSA-2 n = 3 v 3 (~ treatment). Right-hand "
         "asterisks = BH q across the genes in this figure, within each dataset. " + note)
    return s


cross_forest(["Isg15", "Irf7", "Ifit1", "Ifit3", "Ifit3b", "Rsad2", "Isg20", "Stat1", "Irf9", "Cxcl10"],
             "11_interferon_dura_vs_microglia.png", "Interferon-stimulated genes: up in dura, flat in brain microglia",
             note="Dura interferon tracks monocyte content (Fig. 07).")
cross_forest(["Insig1", "Hmgcs1", "Hmgcr", "Sqle", "Idi1", "Msmo1", "Cyp51", "Dhcr24", "Ldlr", "Srebf2"],
             "12_cholesterol_SREBP_both_tissues.png", "SREBP2 / cholesterol-synthesis genes across tissues", include_acsa=True,
             note="Down in dura (after fibroblast adjustment) and microglia; the astrocyte/oligodendrocyte-rich ACSA-2 fraction does NOT show it "
                  "(sterol program trends up ~10% there), so this is not a body-wide effect.")

# ================================================================ 13. CD11b carry-over by arm
pur = pd.read_csv(os.path.join(SORT, "02_purity_celltyping", "purity_carryover_per_sample.tsv"), sep="\t", index_col=0)
dec = pd.read_csv(os.path.join(SORT, "02_purity_celltyping", "deconvolution_rna_share_pct_median.tsv"), sep="\t", index_col=0)
meas = [("astro_marker_pct_vs_ACSA2", pur, "Astrocyte marker share (%)"), ("BAM", dec, "BAM RNA share (%)"),
        ("endothelial", dec, "Endothelial RNA share (%)"), ("microglia", dec, "Microglial RNA share (%)")]
rows = []
for col, src, nm in meas:
    v = src.loc[CD, col].values
    rows.append({"measure": nm, "NTC_mean": v[cd_trt == "NTC"].mean(), "SPP1_mean": v[cd_trt == "SPP1"].mean(),
                 "welch_p": stats.ttest_ind(v[cd_trt == "SPP1"], v[cd_trt == "NTC"], equal_var=False).pvalue,
                 "exact_perm_p": perm_within(v, (cd_trt == "SPP1").astype(int), np.zeros(len(v)))[0]})
s13 = pd.DataFrame(rows)
s13["q_BH"] = bh(s13.welch_p)
s13.to_csv(os.path.join(STAT, "13_microglia_sort_carryover.csv"), index=False)
fig, axs = plt.subplots(1, 4, figsize=(10, 3.9))
for ax, (col, src, nm), (_, r) in zip(axs, meas, s13.iterrows()):
    dot_panel(ax, src.loc[CD, col].values, cd_trt, ["NTC", "SPP1"], ylabel=nm, q=r.q_BH, ptxt=f"q={fmt_p(r.q_BH)}", ylim0=col != "microglia")
fig.suptitle("Every SPP1 CD11b+ sort carries more non-microglial RNA than every NTC sort", fontsize=12, fontweight="bold")
fig.tight_layout(rect=(0, 0.1, 1, 0.95))
save(fig, "13_microglia_sort_carryover.png",
     "Per-sample RNA shares (marker-based astrocyte share; reference deconvolution for BAM, endothelial, microglia), n = 3 v 3. Welch t-test, "
     "BH q across the 4 measures (complete separation in all 4, but with n = 3 v 3 the exact permutation p cannot go below 0.1). Small but consistent: it is why genes from these cells (Col3a1, Ccn1/2, Pf4, Mrc1, Gpr37l1) rise in SPP1 sorts, "
     "and why microglial claims use microglia-intrinsic genes.")

# ================================================================ 14. cross-tissue pathway scatter
X = pd.read_csv(os.path.join(SORT, "14_pathway_links", "pathways_dura_vs_CD11b_signed_log10p.csv"), index_col=0)
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.8))
thr = -np.log10(0.1)
for ax, (xc, yc, ttl) in zip(axs, [("dura_unadj", "CD11b_all", "Unadjusted, both tissues"),
                                   ("dura_fibro_adj", "CD11b_intrinsic", "Dura fibroblast-adjusted vs microglia-intrinsic")]):
    rho, prho = stats.spearmanr(X[xc], X[yc])
    ax.axhline(0, color=INK2, lw=0.9); ax.axvline(0, color=INK2, lw=0.9)
    for v in (thr, -thr):
        ax.axhline(v, color=MUTED, lw=0.7, ls=(0, (3, 2))); ax.axvline(v, color=MUTED, lw=0.7, ls=(0, (3, 2)))
    ax.scatter(X[xc], X[yc], s=24, color="#c9c8c3", edgecolor="white", lw=0.6, zorder=2)
    keyset = ["M_CHOLESTEROL_SREBP", "H_CHOLESTEROL_HOMEOSTASIS", "H_MTORC1_SIGNALING", "C_DAM_KerenShaul2017", "C_LIPID_HANDLING",
              "H_INTERFERON_ALPHA_RESPONSE", "C_TYPE_I_IFN_ISG", "M_TYPE_I_IFN_EXTENDED", "C_ECM_REMODELING", "H_EPITHELIAL_MESENCHYMAL_TRANSITION",
              "M_YAP_TAZ_PATHWAY", "M_TGFB_SMAD_TARGETS", "M_EOSINOPHIL_GRANULE"]
    for k_, s in enumerate(keyset):
        if s in X.index:
            r = X.loc[s]
            same = np.sign(r[xc]) == np.sign(r[yc])
            cc = "#1baf7a" if (same and abs(r[xc]) > thr * 0.6 and abs(r[yc]) > thr * 0.6) else INK
            ax.scatter(r[xc], r[yc], s=40, color=cc, edgecolor="white", lw=0.8, zorder=3)
            ax.annotate(nice(s), (r[xc], r[yc]), xytext=(5, 4 if k_ % 2 == 0 else -10), textcoords="offset points", fontsize=7.2)
    ax.set_xlabel("dura: signed −log10 p")
    ax.set_ylabel("CD11b microglia: signed −log10 p")
    ax.set_title(f"{ttl}\nSpearman ρ = {rho:.2f} (p = {fmt_p(prho)}, {len(X)} sets)", fontsize=10)
fig.tight_layout(rect=(0, 0.07, 1, 1))
X.to_csv(os.path.join(STAT, "14_pathways_dura_vs_microglia.csv"))
save(fig, "14_pathways_dura_vs_microglia.png",
     "Each dot = one gene set; signed −log10 conservative camera p (+ = up with SPP1). Dashed = p 0.1. Green = same direction and p < ~0.25 in both "
     "tissues (cholesterol/SREBP, mTORC1). Spearman p treats sets as independent, which they are not (shared genes), so it is optimistic.")
print("done")
