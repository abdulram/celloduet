"""10: 5xFAD pathology-associated programs, SPP1 vs NTC, per fraction.
Gene sets are fixed in advance from the literature (none selected from this dataset). Spp1 is removed from every set,
so a set score cannot move just because Spp1 was knocked down. Ex vivo dissociation-sensitive genes (Marsh 2022:
Ccl3, Ccl4, Tnf, Il1b, Fos, Jun, Egr1, Hspa1a ...) are kept out of the pathology sets.
Score per animal = mean over genes of log2CPM centred on the NTC mean (=> log2 ratio vs NTC).
Stats: Welch t 95% CI (df~4) and exact two-sided permutation over the 20 relabelings (min p = 0.1)."""
import os, sys, itertools
sys.path.insert(0, os.path.dirname(__file__))
from common import *
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

SETS = {  # (fraction it is read in, genes, source)
 "DAM_stage1":        ("CD11b", "Tyrobp Ctsb Ctsd Apoe B2m Fth1 Lyz2", "Keren-Shaul 2017"),
 "DAM_stage2":        ("CD11b", "Trem2 Axl Cst7 Ctsl Lpl Cd9 Csf1 Ccl6 Itgax Clec7a Lilrb4a Timp2", "Keren-Shaul 2017"),
 "MGnD_up":           ("CD11b", "Apoe Axl Clec7a Itgax Gpnmb Lgals3 Ccl2 Cst7 Lpl Itgb2 Cd9", "Krasemann 2017"),
 "ARM":               ("CD11b", "Cst7 Itgax Clec7a Dkk2 Gpnmb Cd74 H2-Aa H2-Ab1 H2-Eb1 Apoe Ank", "Sala Frigerio 2019"),
 "Homeostatic":       ("CD11b", "P2ry12 P2ry13 Tmem119 Cx3cr1 Sall1 Csf1r Selplg Gpr34 Siglech Olfml3 Hexb Tgfbr1", "Keren-Shaul/Krasemann"),
 "Plaque_induced_PIG":("CD11b", "C1qa C1qb C1qc C4b Ctss Cd68 Trem2 Tyrobp Gpnmb Ctsz Itgax Cd9 Clec7a B2m Cst7 H2-D1 H2-K1 Ctsd Prdx6", "Chen 2020 (microglia-expressed subset)"),
 "Phagocytosis":      ("CD11b", "Cd68 Axl Mertk Cd36 Msr1 Lrp1 Gas6 Trem2 Tyrobp", "curated"),
 "Lysosome":          ("CD11b", "Ctsb Ctsd Ctsl Ctss Ctsz Lamp1 Lamp2 Grn Hexa Hexb Npc2 Atp6v0d2 Gba", "curated"),
 "Complement":        ("CD11b", "C1qa C1qb C1qc C3ar1 C5ar1 Itgb2 C3", "curated"),
 "MHCII_antigen":     ("CD11b", "H2-Aa H2-Ab1 H2-Eb1 Cd74 H2-DMa H2-DMb1", "curated"),
 "Interferon_IRM":    ("CD11b", "Ifit1 Ifit2 Ifit3 Ifit3b Isg15 Irf7 Rtp4 Oasl2 Bst2 Usp18 Stat1 Ifi27l2a", "Sala Frigerio 2019"),
 "Proliferation":     ("CD11b", "Mki67 Top2a Ccnb2 Cdk1 Birc5 Cenpf Ube2c Pclaf", "curated"),
 "Lipid_handling":    ("CD11b", "Apoe Lpl Abca1 Abcg1 Plin2 Lipa Soat1 Cd36 Msr1", "curated (no Ch25h/Fabp5: flagged in pass 1)"),
 "DAA_astro":         ("ACSA2", "Gfap Vim Serpina3n Ctss C4b Cd9 Apoe Osmr Cst3 Clu", "Habib 2020"),
 "PanReactive_astro": ("ACSA2", "Lcn2 Steap4 S1pr3 Timp1 Hspb1 Cxcl10 Cd44 Osmr Cp Serpina3n Aspg Vim Gfap", "Zamanian 2012"),
 "A1_astro":          ("ACSA2", "C3 Serping1 H2-D1 Gbp2 Fkbp5 Ugt1a1 Psmb8 Srgn Amigo2", "Liddelow 2017"),
 "DOL_oligo":         ("ACSA2", "Serpina3n C4b H2-D1 H2-K1 B2m Il33 Klk6 Cdkn1a", "Kenigsbuch 2022 (approx.)"),
 "Homeostatic_astro": ("ACSA2", "Slc1a2 Slc1a3 Kcnj10 Glul Aqp4 Gja1 Sparcl1 Gpc4 Gpc6", "curated"),
 "Cholesterol_astro": ("ACSA2", "Hmgcr Hmgcs1 Sqle Fdft1 Srebf2 Idi1 Cyp51 Dhcr24 Ldlr Insig1", "curated"),
 "Interferon_glia":   ("ACSA2", "Ifit1 Ifit3 Ifit3b Isg15 Irf7 Rtp4 Oasl2 Bst2 Stat1 Gbp2 Iigp1 Igtp", "curated"),
}
ann, cnt = load_matrix(); ss = sample_sheet()
keep = ~ann.index.str.startswith("ERCC")
cpm = cnt[keep] / cnt[keep].sum() * 1e6; cpm.index = ann.gene_name[keep]; cpm = cpm.groupby(level=0).sum()
lc = np.log2(cpm + 1)

def perm_p(x, lab):
    obs = x[lab].mean() - x[~lab].mean(); null = []
    for c in itertools.combinations(range(len(x)), lab.sum()):
        m = np.zeros(len(x), bool); m[list(c)] = True; null.append(x[m].mean() - x[~m].mean())
    return (np.abs(null) >= abs(obs) - 1e-12).mean()

rows, persample = [], {}
for name, (frac, genes, src) in SETS.items():
    s = ss.index[ss.fraction == frac]; lab = (ss.loc[s, "treatment"] == "Spp1").values
    g = [x for x in genes.split() if x in lc.index and x != "Spp1" and cpm.loc[x, s].mean() >= 2]
    c = lc.loc[g, s].sub(lc.loc[g, s[~lab]].mean(1), axis=0)
    sc = c.mean(); persample[name] = sc
    a, b = sc[lab].values, sc[~lab].values
    d = a.mean() - b.mean()
    se = np.sqrt(a.var(ddof=1) / 3 + b.var(ddof=1) / 3)
    dfw = se**4 / ((a.var(ddof=1) / 3)**2 / 2 + (b.var(ddof=1) / 3)**2 / 2)
    t = stats.t.ppf(0.975, dfw)
    rows.append(dict(set=name, fraction=frac, source=src, n_genes=len(g), log2_SPP1_vs_NTC=d,
                     pct_change=(2**d - 1) * 100, ci95_lo_pct=(2**(d - t * se) - 1) * 100, ci95_hi_pct=(2**(d + t * se) - 1) * 100,
                     welch_p=stats.ttest_ind(a, b, equal_var=False).pvalue, perm_p=perm_p(sc.values, lab),
                     separated=(a.min() > b.max()) or (a.max() < b.min()),
                     per_gene_log2FC=";".join(f"{x}:{c.loc[x, s[lab]].mean():+.2f}" for x in g)))
res = pd.DataFrame(rows)
res.to_csv(os.path.join(ROOT, "09_pathology", "pathology_programs_SPP1_vs_NTC.tsv"), sep="\t", index=False, float_format="%.4g")
ps = pd.DataFrame(persample); ps.to_csv(os.path.join(ROOT, "09_pathology", "pathology_program_scores_per_sample.tsv"), sep="\t", float_format="%.3f")
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 200)
print(res.drop(columns="per_gene_log2FC").round(3).to_string(index=False))

# Knockdown index: Spp1 relative to the DAM2 program in the same animal (does Spp1 fall per unit of DAM response?)
cd = ss.index[ss.fraction == "CD11b"]
dam = lc.loc[SETS["DAM_stage2"][1].split(), cd].mean()
ki = pd.DataFrame({"treatment": ss.loc[cd, "treatment"], "log2_Spp1": lc.loc["Spp1", cd], "DAM2_mean_log2": dam,
                   "Spp1_minus_DAM2": lc.loc["Spp1", cd] - dam})
ki.to_csv(os.path.join(ROOT, "09_pathology", "spp1_knockdown_index.tsv"), sep="\t", float_format="%.3f")
print("\n", ki.round(3).to_string())

# Minimal detectable difference for a program score (80% power, alpha 0.05, n=3/3) from pooled SD
for name in ["DAM_stage2", "Plaque_induced_PIG", "Homeostatic", "DAA_astro", "DOL_oligo"]:
    f = SETS[name][0]; s = ss.index[ss.fraction == f]; lab = (ss.loc[s, "treatment"] == "Spp1").values
    sd = np.sqrt((ps.loc[s[lab], name].var(ddof=1) + ps.loc[s[~lab], name].var(ddof=1)) / 2)
    mdd = (stats.t.ppf(0.975, 4) + stats.t.ppf(0.8, 4)) * sd * np.sqrt(2 / 3)
    print(f"MDD {name}: +/-{(2**mdd-1)*100:.0f}% (pooled SD {sd:.3f} log2)")

# Figure: forest plot
fig, ax = plt.subplots(figsize=(7, 7))
r = res.iloc[::-1].reset_index(drop=True)
for i, x in r.iterrows():
    col = "#1baf7a" if x.fraction == "CD11b" else "#7a5cd6"
    ax.plot([x.ci95_lo_pct, x.ci95_hi_pct], [i, i], color=col, lw=2)
    ax.scatter(x["pct_change"], i, color=col, s=40, zorder=3, edgecolor="k")
ax.axvline(0, color="grey", lw=0.8)
ax.set_yticks(range(len(r))); ax.set_yticklabels([f"{a} ({b}, n={n})" for a, b, n in zip(r.set, r.fraction, r.n_genes)], fontsize=7)
ax.set_xlabel("SPP1 vs NTC, % change in program score (95% CI, n=3 vs 3)")
ax.set_title("5xFAD pathology programs (literature sets, Spp1 excluded)", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, "figures", "fig10_pathology_programs_forest.png"), dpi=150)

# ACSA2 sets re-expressed per unit of the host cell type (astrocyte or oligodendrocyte lineage score), because
# astrocyte RNA share is ~20% higher in SPP1 ACSA2 samples (deconvolution). Lineage genes are not in the state sets.
LIN = {"astro": "Aldh1l1 Atp1b2 Aldoc Sox9 Fgfr3 Gjb6 Slc4a4 Ntsr2", "oligo": "Plp1 Mbp Mog Mobp Mag Cldn11 Ermn Opalin"}
ac = ss.index[ss.fraction == "ACSA2"]; lab = (ss.loc[ac, "treatment"] == "Spp1").values
lin = {k: lc.loc[v.split(), ac].sub(lc.loc[v.split(), ac[~lab]].mean(1), axis=0).mean() for k, v in LIN.items()}
rows2 = []
for name, host in [("DAA_astro", "astro"), ("PanReactive_astro", "astro"), ("A1_astro", "astro"), ("Cholesterol_astro", "astro"),
                   ("Homeostatic_astro", "astro"), ("Interferon_glia", "astro"), ("DOL_oligo", "oligo")]:
    adj = ps.loc[ac, name] - lin[host]
    a, b = adj[lab], adj[~lab]
    rows2.append(dict(set=name, normalised_to=host, host_lineage_log2FC=lin[host][lab].mean() - lin[host][~lab].mean(),
                      raw_pct=(2**(ps.loc[ac, name][lab].mean() - ps.loc[ac, name][~lab].mean()) - 1) * 100,
                      per_host_pct=(2**(a.mean() - b.mean()) - 1) * 100, perm_p=perm_p(adj.values, lab),
                      welch_p=stats.ttest_ind(a, b, equal_var=False).pvalue))
r2 = pd.DataFrame(rows2)
r2.to_csv(os.path.join(ROOT, "09_pathology", "acsa2_programs_per_host_cell.tsv"), sep="\t", index=False, float_format="%.4g")
print("\nACSA2 programs normalised to host lineage:\n", r2.round(3).to_string(index=False))
