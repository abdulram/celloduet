"""12: gene-level detail in CD11b+ microglia for MMP/ECM-protease, phagocytosis and lipid-handling families.
Families are defined a priori (functional categories); every member expressed at >=5 CPM is reported, not just hits.
Columns: CPM per arm, DESeq2 log2FC/padj, leave-one-out max p, separation, carry-over prediction, enrichment in
the other fraction (>=4x means likely contaminant-derived)."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

FAM = {
 "MMP": "Mmp2 Mmp3 Mmp8 Mmp9 Mmp12 Mmp13 Mmp14 Mmp15 Mmp16 Mmp17 Mmp19 Mmp24 Mmp25 Mmp28",
 "ADAM/ADAMTS": "Adam8 Adam9 Adam10 Adam15 Adam17 Adam19 Adamts1 Adamts4 Adamts9 Adamts10",
 "TIMP/other protease regulators": "Timp1 Timp2 Timp3 Timp4 Reck Serpine1 Serpine2 Plau Plaur Plat",
 "Cathepsins (ECM/lysosomal)": "Ctsb Ctsd Ctsk Ctsl Ctss Ctsz Ctsh Ctsc",
 "Phago: recognition receptors": "Trem2 Tyrobp Axl Mertk Gas6 Pros1 Cd36 Msr1 Marco Scarb1 Lrp1 Cd14 Itgam Itgb2 Fcgr1 Fcgr2b Fcgr3 Clec7a Cd93 Siglech Cd300a Lilrb4a P2ry6 Megf10",
 "Phago: eat-me/bridging & complement": "C1qa C1qb C1qc C3 C3ar1 C5ar1 Mfge8 Lgals3 Cd47 Sirpa Itgav Itgb5",
 "Phago: engulfment/phagosome": "Elmo1 Dock2 Rac1 Rac2 Cdc42 Cd68 Lamp1 Lamp2 Rab7 Rubcn Ncf1 Ncf2 Ncf4 Cybb Cyba Atp6v0d2 Atp6v1a Tfeb",
 "Lipid: uptake": "Cd36 Msr1 Lrp1 Ldlr Vldlr Lpl Scarb1 Trem2 Olr1",
 "Lipid: apolipoprotein/efflux": "Apoe Apoc1 Apoc2 Clu Abca1 Abca7 Abcg1 Pltp Nr1h3 Nr1h2",
 "Lipid: storage/droplets": "Plin2 Plin3 Soat1 Dgat1 Dgat2 Lipa Pnpla2 Lpcat2 Cidea",
 "Lipid: cholesterol synthesis/oxysterol": "Hmgcr Hmgcs1 Sqle Fdft1 Cyp51 Dhcr7 Dhcr24 Srebf2 Insig1 Ldlr Ch25h Cyp27a1 Cyp7b1 Cyp46a1 Npc1 Npc2",
 "Lipid: fatty acid binding/oxidation/synthesis": "Fabp3 Fabp4 Fabp5 Fabp7 Cpt1a Acadl Acadm Ppara Ppard Pparg Acsl1 Fasn Scd1 Scd2 Srebf1 Elovl1 Elovl6",
 "Lipid: phospholipid/sphingolipid & lysosomal lipid": "Gba Gla Asah1 Smpd1 Psap Grn Pla2g7 Pla2g15 Lpcat2 Plcg2 Sgpl1 Sphk1 Cerk",
 "Lipid/lysosome stress markers": "Gpnmb Lgals3 Cd9 Itgax Igf1 Spp1 Cst7",
}
ROUND1 = {"Spp1", "Gpnmb", "Ch25h", "Fabp5", "Rgs16", "Chst1", "Mrc2", "Lat2", "Hpse", "Lmo2", "Mcoln3", "Lgi2", "Mir142hg"}
ann, cnt = load_matrix(); ss = sample_sheet()
de = pd.read_csv(os.path.join(ROOT, "05_DE_CD11b", "DE_Spp1_vs_NTC_CD11b_annotated.tsv"), sep="\t").dropna(subset=["gene"]).drop_duplicates("gene").set_index("gene")
keep = ~ann.index.str.startswith("ERCC")
cpm = cnt[keep] / cnt[keep].sum() * 1e6; cpm.index = ann.gene_name[keep]; cpm = cpm.groupby(level=0).sum()
cd = ss.index[ss.fraction == "CD11b"]; n = [s for s in cd if ss.loc[s, "treatment"] == "NTC"]; p = [s for s in cd if ss.loc[s, "treatment"] == "Spp1"]
rows = []
for fam, genes in FAM.items():
    for g in genes.split():
        if g not in cpm.index or cpm.loc[g, cd].mean() < 5: continue
        d = de.loc[g] if g in de.index else None
        sep = (cpm.loc[g, p].max() < cpm.loc[g, n].min()) or (cpm.loc[g, p].min() > cpm.loc[g, n].max())
        rows.append(dict(family=fam, gene=g, NTC_cpm=cpm.loc[g, n].mean(), SPP1_cpm=cpm.loc[g, p].mean(),
                         pct_change=(cpm.loc[g, p].mean() / cpm.loc[g, n].mean() - 1) * 100,
                         log2FC=d.log2FC if d is not None else np.nan, padj=d.padj if d is not None else np.nan,
                         pvalue=d.pvalue if d is not None else np.nan,
                         loo_max_p=d.loo_max_pvalue if d is not None else np.nan, separated=sep,
                         carryover_pred_log2FC=d.carryover_pred_log2FC if d is not None else np.nan,
                         ACSA2_over_CD11b=d.ratio_other_over_this if d is not None else np.nan,
                         flagged_in_round1=g in ROUND1))
t = pd.DataFrame(rows)
def call(r):
    if r.ACSA2_over_CD11b >= 4 or (abs(r.carryover_pred_log2FC) >= 0.5 * abs(r.log2FC) and np.sign(r.carryover_pred_log2FC) == np.sign(r.log2FC) and abs(r.log2FC) > 0.2):
        return "carry-over"
    if r.padj < 0.1 and r.loo_max_p < 0.05: return "changed (padj<0.1, LOO-robust)"
    if r.pvalue < 0.05 and r.separated: return "nominal (p<0.05, groups separate)"
    if r.pvalue < 0.05: return "nominal (p<0.05)"
    return "no change"
t["call"] = t.apply(call, axis=1)
t.to_csv(os.path.join(ROOT, "07_modules", "CD11b_MMP_phagocytosis_lipid_gene_detail.tsv"), sep="\t", index=False, float_format="%.4g")
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
for fam, g in t.groupby("family", sort=False):
    print(f"\n### {fam}")
    print(g[["gene", "NTC_cpm", "SPP1_cpm", "pct_change", "padj", "pvalue", "separated", "call", "flagged_in_round1"]].round(3).to_string(index=False))
# figure: per-gene % change, coloured by call
fig, ax = plt.subplots(figsize=(7, 0.16 * len(t) + 1))
col = {"changed (padj<0.1, LOO-robust)": "#c0392b", "nominal (p<0.05, groups separate)": "#e67e22", "nominal (p<0.05)": "#f5b041",
       "no change": "#9aa0a6", "carry-over": "#7a5cd6"}
tt = t.iloc[::-1].reset_index(drop=True)
ax.barh(range(len(tt)), tt["pct_change"], color=[col[c] for c in tt["call"]])
ax.set_yticks(range(len(tt))); ax.set_yticklabels([f"{g} [{f.split(':')[0][:12]}]" for g, f in zip(tt.gene, tt.family)], fontsize=5)
ax.axvline(0, color="k", lw=0.6); ax.set_xlabel("SPP1 vs NTC, % change in CPM (CD11b+, n=3 vs 3)")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=v, label=k) for k, v in col.items()], fontsize=6, loc="lower right"); fig.tight_layout()
fig.savefig(os.path.join(ROOT, "figures", "fig11_mmp_phago_lipid_CD11b.png"), dpi=150)

# SREBP2 / cholesterol-synthesis targets (a priori list; carry-over from astrocytes would push these UP in SPP1)
SREBP2 = "Hmgcr Hmgcs1 Sqle Fdft1 Cyp51 Dhcr7 Dhcr24 Srebf2 Insig1 Ldlr Idi1 Mvd Fdps Lss Msmo1 Nsdhl Mvk Sc5d Tmem97 Stard4".split()
lab = (ss.loc[cd, "treatment"] == "Spp1").values
g = [x for x in SREBP2 if x in cpm.index and cpm.loc[x, cd].mean() >= 5]
lc = np.log2(cpm.loc[g, cd] + 1); sc = lc.sub(lc.loc[:, ~lab].mean(1), axis=0).mean()
import itertools
obs = sc[lab].mean() - sc[~lab].mean()
null = [sc.values[list(c)].mean() - np.delete(sc.values, list(c)).mean() for c in itertools.combinations(range(6), 3)]
sub = de.reindex(g)
print(f"\nSREBP2 targets (n={len(g)}): set score SPP1 vs NTC {(2**obs-1)*100:+.1f}%, exact perm p={np.mean(np.abs(null) >= abs(obs)-1e-12):.2f}; "
      f"genes down {int((sub.log2FC < 0).sum())}/{len(g)}; carry-over predicts mean {sub.carryover_pred_log2FC.mean():+.3f} log2 (up)")
print("per-sample scores:", sc.round(3).to_dict())
