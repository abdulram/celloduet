"""14c - CD11b+ microglia, SPP1 vs NTC: pathway-link map, mechanistic genes next to dura, cross-tissue pathway
comparison, and upload kits for PANTHER and QIAGEN IPA. Mirrors dura step 13c.

Inputs: tables/DE_CD11b_SPP1_vs_NTC_{unadjusted,carry_adj}.tsv, tables/CD11b_pathway_camera_both_universes.tsv (14b),
        dura tables from ../composition_confounding_analysis/tables (step 13).
Link map: nodes = sets with conservative camera p < 0.15 in either universe; edges = Jaccard of leading genes >= 0.1
and >= 2 shared genes. Black ring = also p < 0.15, same direction, in the microglia-intrinsic universe.
"""
import itertools
import os
import re

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TAB = os.path.join(ROOT, "tables")
DURA = os.path.join(os.path.dirname(ROOT), "composition_confounding_analysis")
OUT = os.path.join(ROOT, "14_pathway_links")
KIT = os.path.join(ROOT, "14_pathway_kits")
for d in [OUT, os.path.join(KIT, "panther"), os.path.join(KIT, "ipa")]:
    os.makedirs(d, exist_ok=True)
exec(open(os.path.join(DURA, "scripts", "common.py")).read())  # palette + save()

u = pd.read_csv(os.path.join(TAB, "DE_CD11b_SPP1_vs_NTC_unadjusted.tsv"), sep="\t")
a = pd.read_csv(os.path.join(TAB, "DE_CD11b_SPP1_vs_NTC_carry_adj.tsv"), sep="\t")
cam = pd.read_csv(os.path.join(TAB, "CD11b_pathway_camera_both_universes.tsv"), sep="\t")
du = pd.read_csv(os.path.join(DURA, "tables", "DE_SPP1_vs_NTC_5xFADonly.tsv"), sep="\t")
df = pd.read_csv(os.path.join(DURA, "tables", "DE_SPP1_vs_NTC_fibroblast_adjusted.tsv"), sep="\t")
dcam = pd.read_csv(os.path.join(DURA, "tables", "spp1_pathway_camera_both_models.tsv"), sep="\t")

# ---------------------------------------------------------------- 1. upload kits
clean = lambda s: s.str.replace(r"\..*$", "", regex=True)
u["id"] = clean(u.gene_id)
counts = []
for univ, d in [("all", u), ("mg_intrinsic", u[u.microglia_intrinsic])]:
    d = d.dropna(subset=["pvalue"])
    d.id.to_csv(os.path.join(KIT, "panther", f"REFERENCE_{univ}_{len(d)}_EnsemblIDs.txt"), index=False, header=False)
    for cut, tag in [(0.01, "p01"), (0.05, "p05")]:
        up = d[(d.pvalue < cut) & (d.log2FoldChange > 0)].id
        dn = d[(d.pvalue < cut) & (d.log2FoldChange < 0)].id
        up.to_csv(os.path.join(KIT, "panther", f"{univ}_UP_{tag}_EnsemblIDs.txt"), index=False, header=False)
        dn.to_csv(os.path.join(KIT, "panther", f"{univ}_DOWN_{tag}_EnsemblIDs.txt"), index=False, header=False)
        counts.append({"universe": univ, "cutoff": f"p<{cut}", "n_up": len(up), "n_down": len(dn), "n_reference": len(d)})
    d[["id", "stat"]].dropna().to_csv(os.path.join(KIT, "panther", f"{univ}_RANKED_wald_stat.txt"), sep="\t", index=False, header=False)
pd.DataFrame(counts).to_csv(os.path.join(KIT, "panther", "list_sizes.csv"), index=False)
print(pd.DataFrame(counts).to_string())
ai = a.set_index("gene_id")
ipa = pd.DataFrame({"Ensembl_ID": u.id, "Symbol": u.gene_name,
                    "log2FC_unadjusted": u.log2FoldChange, "pvalue_unadjusted": u.pvalue, "padj_unadjusted": u.padj,
                    "log2FC_carryover_adjusted": ai.reindex(u.gene_id).log2FoldChange.values,
                    "pvalue_carryover_adjusted": ai.reindex(u.gene_id).pvalue.values,
                    "microglia_intrinsic": u.microglia_intrinsic.astype(int), "baseMean": u.baseMean})
ipa.to_csv(os.path.join(KIT, "ipa", "IPA_upload_SPP1_vs_NTC_5xFAD_CD11b_microglia.txt"), sep="\t", index=False, float_format="%.6g")

# ---------------------------------------------------------------- 2. mechanistic genes: dura vs microglia
MECH_GENES = {"SPP1 receptors / integrin-FAK": ["Spp1", "Cd44", "Itgav", "Itgb1", "Itgb3", "Itgb5", "Ptk2"],
              "DAM / lipid-associated": ["Gpnmb", "Apoe", "Ch25h", "Fabp5", "Igf1", "Lpl", "Cst7", "Itgax", "Trem2", "Tyrobp", "Clec7a"],
              "Homeostatic": ["P2ry12", "Tmem119", "Cx3cr1", "Sall1"],
              "YAP/TAZ (Hippo)": ["Yap1", "Wwtr1", "Tead1", "Amot", "Ccn1", "Ccn2"],
              "TGF-beta / SMAD": ["Tgfb1", "Tgfbr1", "Tgfbr2", "Smad3", "Smad7", "Skil", "Bhlhe40"],
              "Matrix output": ["Fn1", "Col3a1", "Lox", "Timp1", "Prg4", "Cldn11"],
              "Interferon / JAK-STAT": ["Rsad2", "Isg20", "Ifit3b", "Ifit1", "Ifit3", "Isg15", "Irf7", "Stat1", "Irf9", "Cxcl10"],
              "Cholesterol / SREBP": ["Insig1", "Hmgcs1", "Hmgcr", "Sqle", "Idi1", "Msmo1", "Cyp51", "Dhcr24", "Ldlr", "Srebf2"],
              "Eosinophil": ["Prg2", "Epx", "Ear2", "Alox15"]}
idx = lambda d: d.drop_duplicates("gene_name").set_index("gene_name")
U, A, DU, DF = idx(u), idx(a), idx(du), idx(df)
rows = []
for mod, gl in MECH_GENES.items():
    for g in gl:
        r = {"module": mod, "gene": g}
        for tag, d in [("dura_unadj", DU), ("dura_fibro_adj", DF), ("CD11b_unadj", U), ("CD11b_carry_adj", A)]:
            r[f"log2FC_{tag}"] = d.log2FoldChange.get(g, np.nan)
            r[f"p_{tag}"] = d.pvalue.get(g, np.nan)
        r["CD11b_baseMean"] = U.baseMean.get(g, np.nan)
        r["CD11b_microglia_intrinsic"] = U.microglia_intrinsic.get(g, np.nan)
        r["CD11b_excluded_because"] = U.reason_excluded.get(g, "not expressed") if g in U.index else "not expressed"
        rows.append(r)
mt = pd.DataFrame(rows)
mt.to_csv(os.path.join(OUT, "mechanistic_genes_dura_vs_CD11b.csv"), index=False)

cols = ["dura_unadj", "dura_fibro_adj", "CD11b_unadj", "CD11b_carry_adj"]
M = mt[[f"log2FC_{c}" for c in cols]].values
P = mt[[f"p_{c}" for c in cols]].values
fig, ax = plt.subplots(figsize=(9, 18))
im = ax.imshow(np.nan_to_num(M), cmap=DIV, vmin=-1, vmax=1, aspect="auto")
lab = []
for _, r in mt.iterrows():
    flag = "" if r.CD11b_microglia_intrinsic == True else (" (n.e.)" if r.CD11b_excluded_because == "not expressed" else " †")
    lab.append(f"{r.gene}{flag}  [{r.module}]")
ax.set_yticks(range(len(mt)), lab, fontsize=7.5)
ax.set_xticks(range(4), ["dura\nunadjusted", "dura\nfibroblast-adj", "CD11b microglia\nunadjusted", "CD11b microglia\ncarry-over-adj"], fontsize=8.5)
for i in range(len(mt)):
    for j in range(4):
        v, p = M[i, j], P[i, j]
        if np.isnan(v):
            ax.text(j, i, "n.e.", ha="center", va="center", fontsize=6.5, color=INK2)
            continue
        ax.text(j, i, f"{v:+.2f}{'**' if p < 0.01 else ('*' if p < 0.05 else '')}", ha="center", va="center", fontsize=7,
                color="white" if abs(v) > 0.7 else INK)
ax.axvline(1.5, color=INK, lw=1.5)
ax.grid(False)
fig.colorbar(im, ax=ax, shrink=0.3, label="log2FC SPP1 vs NTC (* p<0.05, ** p<0.01)")
ax.set_title("SPP1 vs NTC: mechanistic genes in dura (n=6v6) and sorted microglia (n=3v3)\n"
             "† = CD11b value may be carry-over (not microglia-intrinsic); n.e. = not expressed/tested", fontsize=9.5)
save(fig, os.path.join(OUT, "mechanistic_genes_dura_vs_CD11b.png"))

# ---------------------------------------------------------------- 3. cross-tissue pathway comparison
def signed(c, model):
    x = c[c.model == model].set_index("set")
    return (-np.log10(x.p_conservative) * np.where(x.direction == "Up", 1, -1)).rename(model)
X = pd.concat([signed(dcam, "unadjusted").rename("dura_unadj"), signed(dcam, "fibro_adj").rename("dura_fibro_adj"),
               signed(cam, "all_genes").rename("CD11b_all"), signed(cam, "mg_intrinsic").rename("CD11b_intrinsic")], axis=1)
X = X.dropna(subset=["dura_fibro_adj", "CD11b_intrinsic"])
X = X[~X.index.str.startswith("TECHNICAL")]
thr = -np.log10(0.1)
X["concordant_p<0.1"] = (np.sign(X.dura_fibro_adj) == np.sign(X.CD11b_intrinsic)) & (X.dura_fibro_adj.abs() > thr) & (X.CD11b_intrinsic.abs() > thr)
X.sort_values("CD11b_intrinsic").to_csv(os.path.join(OUT, "pathways_dura_vs_CD11b_signed_log10p.csv"), float_format="%.3f")
print("\nCross-tissue concordant (p<0.1 both, same direction):\n", X[X["concordant_p<0.1"]].round(2).to_string())
print("Spearman dura_fibro_adj vs CD11b_intrinsic:", X[["dura_fibro_adj", "CD11b_intrinsic"]].corr("spearman").iloc[0, 1].round(2),
      "| dura_unadj vs CD11b_all:", X[["dura_unadj", "CD11b_all"]].corr("spearman").iloc[0, 1].round(2))

fig, axs = plt.subplots(1, 2, figsize=(15, 7))
for ax, (xc, yc, ttl) in zip(axs, [("dura_unadj", "CD11b_all", "unadjusted (both tissues)"),
                                    ("dura_fibro_adj", "CD11b_intrinsic", "dura fibroblast-adjusted vs microglia-intrinsic genes")]):
    ax.scatter(X[xc], X[yc], s=18, c=GRID, edgecolor=INK2, lw=0.4)
    ax.axhline(0, color=INK2, lw=0.6); ax.axvline(0, color=INK2, lw=0.6)
    for v in (thr, -thr):
        ax.axhline(v, color=GRID, ls="--", lw=0.8); ax.axvline(v, color=GRID, ls="--", lw=0.8)
    hit = X[(X[xc].abs() > thr) | (X[yc].abs() > thr)]
    for s, r in hit.iterrows():
        both = abs(r[xc]) > thr and abs(r[yc]) > thr and np.sign(r[xc]) == np.sign(r[yc])
        ax.annotate(re.sub(r"^(H_|C_|M_|G_)", "", s), (r[xc], r[yc]), fontsize=6.5 if not both else 8,
                    color=INK if both else INK2, fontweight="bold" if both else "normal", xytext=(3, 2), textcoords="offset points")
        if both:
            ax.scatter(r[xc], r[yc], s=40, c="#e34948" if r[xc] > 0 else "#2a78d6", zorder=3)
    ax.set_xlabel(f"dura: signed -log10 conservative camera p ({xc})")
    ax.set_ylabel(f"CD11b microglia: signed -log10 p ({yc})")
    ax.set_title(ttl, fontsize=10)
fig.suptitle("SPP1 vs NTC pathway effects: dura vs sorted brain microglia (+ = up with SPP1; dashed = p 0.1; bold = same direction, p<0.1 in both)", fontsize=10)
save(fig, os.path.join(OUT, "pathways_dura_vs_CD11b.png"))

# ---------------------------------------------------------------- 4. link map (CD11b)
camw = cam.pivot(index="set", columns="model", values="p_conservative")
dirn = cam.pivot(index="set", columns="model", values="direction")
lead = cam.pivot(index="set", columns="model", values="leading_genes").fillna("")
keep = camw.index[camw.min(axis=1) < 0.15]
lead_sets = {s: set(";".join(lead.loc[s].values).split(";")) - {""} for s in keep}
G = nx.Graph()
for s in keep:
    best = camw.loc[s].idxmin()
    G.add_node(s, p=camw.loc[s].min(), model=best, direction=dirn.loc[s, best],
               intrinsic_robust=bool(camw.loc[s, "mg_intrinsic"] < 0.15))
for x1, x2 in itertools.combinations(keep, 2):
    A1, A2 = lead_sets[x1], lead_sets[x2]
    if A1 and A2:
        inter = A1 & A2
        jac = len(inter) / len(A1 | A2)
        if jac >= 0.1 and len(inter) >= 2:
            G.add_edge(x1, x2, weight=jac, shared=";".join(sorted(inter)))
edges = pd.DataFrame([{"set_a": x1, "set_b": x2, "jaccard": d["weight"], "shared_leading_genes": d["shared"]} for x1, x2, d in G.edges(data=True)])
edges.sort_values("jaccard", ascending=False).to_csv(os.path.join(OUT, "CD11b_pathway_link_edges.csv"), index=False)
nodes = pd.DataFrame([{"set": n, **d} for n, d in G.nodes(data=True)])
nodes["cluster"] = -1
for i, comp in enumerate(sorted(nx.connected_components(G), key=len, reverse=True)):
    nodes.loc[nodes.set.isin(comp), "cluster"] = i
nodes.to_csv(os.path.join(OUT, "CD11b_pathway_link_nodes.csv"), index=False)
pos = nx.spring_layout(G, weight="weight", seed=7, k=1.1)
fig, ax = plt.subplots(figsize=(15, 11))
for x1, x2, d in G.edges(data=True):
    ax.plot([pos[x1][0], pos[x2][0]], [pos[x1][1], pos[x2][1]], color="#c9c8c3", lw=0.6 + 5 * d["weight"], zorder=1)
for n, d in G.nodes(data=True):
    col = "#e34948" if d["direction"] == "Up" else "#2a78d6"
    size = 120 + 900 * min(-np.log10(d["p"]), 3) / 3
    ax.scatter(*pos[n], s=size, c=col, edgecolor=INK if d["intrinsic_robust"] else SURFACE, lw=2.2 if d["intrinsic_robust"] else 1, zorder=3)
    lab = re.sub(r"^H_", "H: ", re.sub(r"^C_", "", re.sub(r"^M_", "★ ", re.sub(r"^G_", "mg: ", n))))
    ax.annotate(lab, pos[n], xytext=(6, 5), textcoords="offset points", fontsize=8)
ax.axis("off")
ax.set_title("SPP1 vs NTC (5xFAD sorted CD11b+ microglia, n=3v3): pathways linked by shared driver genes\n"
             "red = up, blue = down with SPP1; size = -log10 conservative camera p; thick black ring = holds in microglia-intrinsic genes\n"
             "(i.e. not explainable by glial/BAM/vascular carry-over); ★ = mechanistic module; mg: = microglia module", fontsize=10)
save(fig, os.path.join(OUT, "CD11b_pathway_link_map.png"))
print(nodes.sort_values(["cluster", "p"]).to_string())
print(mt[["gene", "log2FC_dura_fibro_adj", "log2FC_CD11b_unadj", "p_CD11b_unadj", "log2FC_CD11b_carry_adj", "CD11b_excluded_because"]].round(3).to_string())
