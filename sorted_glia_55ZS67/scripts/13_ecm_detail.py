"""13: extracellular-matrix (matrisome) genes, SPP1 vs NTC, both fractions.
Families follow the Naba/MatrisomeDB divisions (core matrisome + ECM regulators/affiliated) plus CNS-specific
perineuronal-net/hyaluronan genes and ECM receptors. Every member expressed at >= 5 CPM in the fraction is reported.
Each gene is tagged with its dominant brain cell type in the Census reference so composition shifts (astrocyte share
up, ependymal/microglial carry-over down in SPP1 ACSA-2; astro/BAM carry-over up in SPP1 CD11b) can be recognised."""
import os, sys, itertools
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch

FAM = {
 "Collagens": "Col1a1 Col1a2 Col3a1 Col4a1 Col4a2 Col4a5 Col5a1 Col5a2 Col5a3 Col6a1 Col6a2 Col6a3 Col8a1 Col9a3 Col11a1 Col11a2 Col12a1 Col13a1 Col15a1 Col16a1 Col18a1 Col19a1 Col20a1 Col23a1 Col25a1 Col27a1",
 "Basement membrane/glycoproteins": "Lama1 Lama2 Lama4 Lama5 Lamb1 Lamb2 Lamc1 Lamc3 Nid1 Nid2 Hspg2 Fn1 Vtn Fbln1 Fbln2 Fbln5 Efemp1 Emilin1 Fbn1 Eln Mfap4 Tinagl1 Smoc1 Smoc2 Sned1 Matn2 Matn4 Ecm1 Igfbp7",
 "Matricellular": "Spp1 Sparc Sparcl1 Thbs1 Thbs2 Thbs3 Thbs4 Tnc Tnr Tnn Postn Ccn1 Ccn2 Ccn3 Lgals1 Lgals3 Lgals9 Cthrc1 Mgp",
 "Proteoglycans (CNS/PNN)": "Acan Bcan Ncan Vcan Hapln1 Hapln2 Hapln4 Ptprz1 Cspg4 Cspg5 Dcn Bgn Lum Ogn Prelp Fmod Aspn Prg4 Srgn Gpc4 Gpc5 Gpc6 Sdc1 Sdc2 Sdc3 Sdc4 Cd44 Tnr",
 "GAG/hyaluronan synthesis & sulfation": "Has1 Has2 Has3 Hyal1 Hyal2 Cemip Cemip2 Tmem2 Chst1 Chst2 Chst3 Chst7 Chst11 Chst12 Chst14 Chst15 Ndst1 Ndst2 Hs3st1 Hs6st1 Hs6st2 Glce Xylt1 Xylt2 Csgalnact1 Chsy1 Ext1 Ext2 Ust Hs2st1 Sulf1 Sulf2",
 "ECM degradation": "Hpse Mmp2 Mmp12 Mmp14 Adamts1 Adamts4 Adamts5 Adamts9 Adamts15 Adamts18 Ctsk Ctsb Ctss Ctsl Plau Plat Mrc2 Arsb Gusb Naglu Hexa Hexb Galns Idua Ids",
 "ECM crosslinking/collagen processing": "Lox Loxl1 Loxl2 Loxl3 Loxl4 P4ha1 P4ha2 P3h1 Plod1 Plod2 Plod3 Serpinh1 Tgm2 Pcolce Pcolce2 Bmp1 Adamts2 Pxdn",
 "ECM inhibitors/regulators": "Timp1 Timp2 Timp3 Timp4 Reck Serpine1 Serpine2 Serpinf1 Cst3 Ambp Itih5",
 "ECM receptors": "Itga1 Itga4 Itga5 Itga6 Itga7 Itga9 Itgav Itgb1 Itgb3 Itgb4 Itgb5 Itgb8 Itgam Itgax Itgb2 Cd44 Ddr1 Ddr2 Dag1 Lrp1 Hmmr Lyve1 Stab1 Cd47 Rpsa",
 "ECM-related TGF-beta": "Tgfb1 Tgfb2 Tgfb3 Ltbp1 Ltbp3 Ltbp4 Tgfbr1 Tgfbr2 Tgfbr3 Thbs1",
}
ann, cnt = load_matrix(); ss = sample_sheet()
keep = ~ann.index.str.startswith("ERCC")
cpm = cnt[keep] / cnt[keep].sum() * 1e6
sym2id = ann.loc[keep].reset_index().dropna(subset=["gene_name"]).drop_duplicates("gene_name").set_index("gene_name")["gene_id"]
ref = pd.read_csv(os.path.join(DATA, "reference", "brain_pseudobulk_cp10k.csv.gz"), index_col=0)
ref.index = ref.index.str.split(".").str[0]; ref = ref.groupby(level=0).sum()
ROUND1 = {"Spp1", "Gpnmb", "Ch25h", "Fabp5", "Rgs16", "Chst1", "Mrc2", "Lat2", "Hpse", "Lmo2", "Mcoln3", "Lgi2", "Mir142hg", "Vtn", "Sparcl1", "Ptn"}

def dominant(gid):
    if gid not in ref.index or ref.loc[gid].max() <= 0: return ("NA", np.nan)
    r = ref.loc[gid].sort_values(ascending=False)
    return (r.index[0], r.iloc[0] / max(r.iloc[1], 1e-3))

allrows = []
for frac, d in [("CD11b", "05_DE_CD11b"), ("ACSA2", "06_DE_ACSA2")]:
    de = pd.read_csv(os.path.join(ROOT, d, f"DE_Spp1_vs_NTC_{frac}_annotated.tsv"), sep="\t").set_index("gene_id")
    s = ss.index[ss.fraction == frac]
    n = [x for x in s if ss.loc[x, "treatment"] == "NTC"]; p = [x for x in s if ss.loc[x, "treatment"] == "Spp1"]
    seen = set()
    for fam, genes in FAM.items():
        for g in genes.split():
            if g not in sym2id.index: continue
            gid = sym2id[g]
            if cpm.loc[gid, s].mean() < 5 or (g, fam) in seen: continue
            seen.add((g, fam))
            r = de.loc[gid] if gid in de.index else None
            dom, dratio = dominant(gid)
            sep = (cpm.loc[gid, p].max() < cpm.loc[gid, n].min()) or (cpm.loc[gid, p].min() > cpm.loc[gid, n].max())
            allrows.append(dict(fraction=frac, family=fam, gene=g, NTC_cpm=cpm.loc[gid, n].mean(), SPP1_cpm=cpm.loc[gid, p].mean(),
                pct_change=(cpm.loc[gid, p].mean() / cpm.loc[gid, n].mean() - 1) * 100,
                log2FC=r.log2FC if r is not None else np.nan, pvalue=r.pvalue if r is not None else np.nan,
                padj=r.padj if r is not None else np.nan, loo_max_p=r.loo_max_pvalue if r is not None else np.nan,
                separated=sep, carryover_pred_log2FC=r.carryover_pred_log2FC if r is not None else np.nan,
                other_over_this=r.ratio_other_over_this if r is not None else np.nan,
                dominant_brain_celltype=dom, dominant_fold_over_next=dratio, flagged_in_round1=g in ROUND1))
t = pd.DataFrame(allrows)
VASC = ("VLMC_fibroblast", "pericyte", "endothelial", "smooth_muscle", "ependymal")
def call(r):
    if r.fraction == "CD11b":
        if r.other_over_this >= 4 or (np.sign(r.carryover_pred_log2FC) == np.sign(r.log2FC) and abs(r.carryover_pred_log2FC) >= 0.5 * abs(r.log2FC) and abs(r.log2FC) > 0.2):
            return "carry-over"
        # vascular/meningeal genes absent from both sorts' own programs: SPP1 CD11b sorts carry more vascular RNA (PC1)
        if r.dominant_brain_celltype in VASC and r.log2FC > 0 and r.NTC_cpm < 60:
            return "carry-over"
    if r.fraction == "ACSA2":
        if r.dominant_brain_celltype in ("ependymal", "microglia", "BAM", "monocyte") and r.log2FC < 0:
            return "composition (ependymal/microglial share down)"
        if r.dominant_brain_celltype in ("VLMC_fibroblast", "pericyte") and r.log2FC < 0:
            return "composition (fibroblast/VLMC share down)"
        if r.dominant_brain_celltype in ("astrocyte", "OPC") and 0 < r.log2FC < 0.45:
            return "composition (astrocyte share up ~+0.24 log2)"
    if r.padj < 0.1 and r.loo_max_p < 0.05: return "changed (padj<0.1, LOO-robust)"
    if r.pvalue < 0.05 and r.separated: return "nominal (p<0.05, groups separate)"
    if r.pvalue < 0.05: return "nominal (p<0.05)"
    return "no change"
t["call"] = t.apply(call, axis=1)
# microglia-intrinsic flag for CD11b: expressed at least 2x higher in CD11b than ACSA2 (reference labels reflect healthy brain)
t["microglia_intrinsic"] = (t.fraction == "CD11b") & (t.other_over_this < 0.5)
# ACSA2: change per astrocyte (astrocyte lineage score +0.243 log2 in SPP1; scripts/10)
t["pct_change_per_astrocyte"] = np.where(t.fraction == "ACSA2", (2 ** (t.log2FC - 0.243) - 1) * 100, np.nan)
t.to_csv(os.path.join(ROOT, "07_modules", "ECM_gene_detail_both_fractions.tsv"), sep="\t", index=False, float_format="%.4g")
pd.set_option("display.width", 260); pd.set_option("display.max_rows", 600)
cols = ["gene", "NTC_cpm", "SPP1_cpm", "pct_change", "pct_change_per_astrocyte", "padj", "pvalue", "dominant_brain_celltype", "microglia_intrinsic", "call", "flagged_in_round1"]
for frac in ["CD11b", "ACSA2"]:
    tf = t[t.fraction == frac]
    print(f"\n######## {frac}: {tf.gene.nunique()} ECM genes expressed; nominal p<0.05: {int((tf.drop_duplicates('gene').pvalue < 0.05).sum())}")
    print(tf.drop_duplicates("gene")["call"].value_counts().to_string())
    for fam, g in tf.groupby("family", sort=False):
        g = g[(g.pvalue < 0.05) | (g["pct_change"].abs() > 25) | g.flagged_in_round1]
        if len(g): print(f"--- {fam}"); print(g[cols].round(3).to_string(index=False))

# Set-level: microglia-expressed ECM genes in CD11b (dominant cell type microglia/BAM/monocyte, not carry-over)
cd = ss.index[ss.fraction == "CD11b"]; lab = (ss.loc[cd, "treatment"] == "Spp1").values
tm = t[(t.fraction == "CD11b") & t.microglia_intrinsic & (t.call != "carry-over")].drop_duplicates("gene")
ids = [sym2id[g] for g in tm.gene if g != "Spp1"]
lc = np.log2(cpm.loc[ids, cd] + 1); sc = lc.sub(lc.loc[:, ~lab].mean(1), axis=0).mean()
obs = sc[lab].mean() - sc[~lab].mean()
null = [sc.values[list(c)].mean() - np.delete(sc.values, list(c)).mean() for c in itertools.combinations(range(6), 3)]
print(tm[tm.pvalue < 0.05][["gene", "family", "pct_change", "padj", "flagged_in_round1"]].round(3).to_string(index=False))
print(f"\nMicroglia-intrinsic ECM genes in CD11b (n={len(ids)}, Spp1 excluded): set {(2**obs-1)*100:+.1f}%, perm p={np.mean(np.abs(null) >= abs(obs)-1e-12):.2f}; "
      f"down {int((tm[tm.gene!='Spp1'].log2FC < 0).sum())}/{len(ids)}")

# Figure: both fractions
fig, axes = plt.subplots(1, 2, figsize=(11, 0.13 * t.groupby("fraction").size().max() + 1.5))
col = {"changed (padj<0.1, LOO-robust)": "#c0392b", "nominal (p<0.05, groups separate)": "#e67e22", "nominal (p<0.05)": "#f5b041",
       "no change": "#9aa0a6", "carry-over": "#7a5cd6", "composition (ependymal/microglial share down)": "#5dade2",
       "composition (astrocyte share up ~+0.24 log2)": "#48c9b0", "composition (fibroblast/VLMC share down)": "#a3e4d7"}
for ax, frac in zip(axes, ["CD11b", "ACSA2"]):
    tt = t[t.fraction == frac].iloc[::-1].reset_index(drop=True)
    ax.barh(range(len(tt)), tt["pct_change"].clip(-80, 150), color=[col[c] for c in tt["call"]])
    ax.set_yticks(range(len(tt))); ax.set_yticklabels([f"{g} [{f[:14]}]" for g, f in zip(tt.gene, tt.family)], fontsize=4.5)
    ax.axvline(0, color="k", lw=0.6); ax.set_title(f"{frac}+ ECM genes, SPP1 vs NTC", fontsize=9); ax.set_xlabel("% change in CPM (clipped -80..150)")
axes[1].legend(handles=[Patch(color=v, label=k) for k, v in col.items()], fontsize=6, loc="lower right")
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "figures", "fig12_ECM_both_fractions.png"), dpi=150)
