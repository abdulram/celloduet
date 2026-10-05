# Shared plotting helpers for the Python steps (exec'd by each script).
import matplotlib.pyplot as _plt
from matplotlib.colors import LinearSegmentedColormap as _LSC
import numpy as _np

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
GROUP_COLORS = {"WT_NTC": "#2a78d6", "5xFAD_NTC": "#eb6834", "5xFAD_SPP1": "#1baf7a"}
GENO_COLORS = {"WT": "#2a78d6", "5xFAD": "#eb6834"}
GROUP_ORDER = list(GROUP_COLORS)
SEX_MARKER = {"F": "o", "M": "^"}
DIV = _LSC.from_list("bgr", ["#184f95", "#3987e5", "#86b6ef", "#f0efec", "#f0a3a2", "#e34948", "#a8322f"])
SEQ = _LSC.from_list("blue", ["#f0efec", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
_plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "font.size": 10, "axes.titlesize": 11, "axes.titleweight": "bold", "legend.frameon": False,
})


def save(fig, path):
    fig.savefig(path, dpi=150, bbox_inches="tight")
    _plt.close(fig)


def strip(ax, vals, meta, groups=None, label_animals=False):
    """Individual animals by group (circle = F, triangle = M), median line."""
    groups = groups or GROUP_ORDER
    for i, g in enumerate(groups):
        for sx, mk in SEX_MARKER.items():
            idx = meta.index[(meta.group == g) & (meta.sex == sx)]
            xs = _np.full(len(idx), i) + (-0.12 if sx == "F" else 0.12)
            ax.scatter(xs, vals[idx], c=GROUP_COLORS[g], marker=mk, s=40, edgecolor=SURFACE, lw=1.1, zorder=3)
            if label_animals:
                for x, l in zip(xs, idx):
                    ax.annotate(str(meta.animal[l]), (x, vals[l]), xytext=(4, 0), textcoords="offset points", fontsize=6.5, color=INK2)
        ax.hlines(_np.median(vals[meta.index[meta.group == g]]), i - 0.3, i + 0.3, color=INK, lw=1.4)
    ax.set_xticks(range(len(groups)), [g.replace("_", "\n") for g in groups], fontsize=8)
    ax.grid(axis="x", visible=False)
