"""09: per-animal CPM plots for Spp1, Spp1-co-regulated candidates, DAM controls, glucocorticoid and carry-over genes."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ann, cnt = load_matrix(); ss = sample_sheet()
cpm = cnt / cnt.sum() * 1e6; cpm.index = ann.gene_name; cpm = cpm.groupby(level=0).sum()
panels = {
 "fig05_spp1_and_candidates_CD11b": ("CD11b", ["Spp1", "Gpnmb", "Ch25h", "Fabp5", "Rgs16", "Lgi2", "Mcoln3", "Hpse", "Apoe", "Igf1"]),
 "fig06_DAM_controls_CD11b": ("CD11b", ["Cst7", "Itgax", "Clec7a", "Lpl", "Cd9", "Trem2", "Tyrobp", "Axl", "P2ry12", "Tmem119"]),
 "fig07_spp1_axis_ACSA2": ("ACSA2", ["Spp1", "Hexb", "Cx3cr1", "Cd44", "Itgav", "Itgb1", "Itgb5", "Gfap", "Serpina3n", "C4b"]),
 "fig08_glucocorticoid_both": (None, ["Fkbp5", "Tsc22d3", "Sgk1", "Klf9", "Per1", "Mt2", "Cdkn1a", "Dbp"]),
 "fig09_carryover_CD11b": ("CD11b", ["Plp1", "Mbp", "Slc1a2", "Atp1a2", "Gpr37l1", "Sparcl1", "S100a8", "Ngp", "Pf4", "Cd163"]),
}
for fname, (frac, genes) in panels.items():
    fr = [frac] if frac else ["CD11b", "ACSA2"]
    fig, axes = plt.subplots(len(fr), len(genes), figsize=(1.55 * len(genes), 2.6 * len(fr)), squeeze=False)
    for i, f in enumerate(fr):
        s = ss.index[ss.fraction == f]
        for j, g in enumerate(genes):
            ax = axes[i, j]
            for k, t in enumerate(["NTC", "Spp1"]):
                st = [x for x in s if ss.loc[x, "treatment"] == t]
                y = cpm.loc[g, st].values
                ax.scatter(np.full(3, k) + np.array([-0.08, 0, 0.08]), y, marker="^", c=COLORS[t], edgecolor="k", s=40, zorder=3)
                ax.hlines(y.mean(), k - 0.25, k + 0.25, color="k", lw=1)
                for xx, yy, nm in zip(np.full(3, k) + np.array([-0.08, 0, 0.08]), y, st):
                    ax.annotate(nm.split("_")[1], (xx, yy), fontsize=5, xytext=(3, 0), textcoords="offset points")
            ax.set_xticks([0, 1]); ax.set_xticklabels(["NTC", "SPP1"], fontsize=7); ax.set_xlim(-0.5, 1.5)
            ax.set_title(f"{g}" + (f"\n{f}" if frac is None else ""), fontsize=8, style="italic"); ax.tick_params(labelsize=6)
            ax.set_ylim(bottom=0)
            if j == 0: ax.set_ylabel("CPM", fontsize=7)
    fig.suptitle(f"{fname.split('_',1)[1].replace('_',' ')}  (genotype not stated; all male by Xist/chrY; triangle = one sample)", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(ROOT, "figures", fname + ".png"), dpi=150); plt.close(fig)
print("ok")
