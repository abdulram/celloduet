"""13c - SPP1 vs NTC (5xFAD): pathway-link map, mechanistic gene table, eosinophil sanity check, and upload kits
for PANTHER (pantherdb.org) and QIAGEN IPA.

Inputs: tables/DE_SPP1_vs_NTC_5xFADonly.tsv, tables/DE_SPP1_vs_NTC_fibroblast_adjusted.tsv,
        tables/spp1_pathway_camera_both_models.tsv, tables/spp1_pathway_sets_long.tsv
Outputs: 13_pathway_links/ (figures, tables), 13_pathway_kits/panther/, 13_pathway_kits/ipa/
Link map: nodes = gene sets with conservative camera p < 0.15 in either model; edge weight = Jaccard overlap of
leading genes (members with DESeq2 p < 0.05 in the set's direction), drawn when Jaccard >= 0.1 AND >= 2 shared genes.
"""
import itertools
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
CA = os.path.dirname(HERE)
TAB = os.path.join(CA, "tables")
OUT = os.path.join(CA, "13_pathway_links")
KIT = os.path.join(CA, "13_pathway_kits")
for d in [OUT, os.path.join(KIT, "panther"), os.path.join(KIT, "ipa")]:
    os.makedirs(d, exist_ok=True)
exec(open(os.path.join(HERE, "common.py")).read())

de = {"unadjusted": pd.read_csv(os.path.join(TAB, "DE_SPP1_vs_NTC_5xFADonly.tsv"), sep="\t"),
      "fibro_adj": pd.read_csv(os.path.join(TAB, "DE_SPP1_vs_NTC_fibroblast_adjusted.tsv"), sep="\t")}
cam = pd.read_csv(os.path.join(TAB, "spp1_pathway_camera_both_models.tsv"), sep="\t")
meta = pd.read_csv(os.path.join(TAB, "sample_metadata_with_composition.tsv"), sep="\t", index_col=0)
S = pd.read_csv(os.path.join(TAB, "composition_signature_scores_all.tsv"), sep="\t", index_col=0)
vst = pd.read_csv(os.path.join(CA, "01_QC", "vst_blind.tsv.gz"), sep="\t", index_col=0)
F5 = meta.index[meta.genotype == "5xFAD"]

# ---------------------------------------------------------------- 1. upload kits
def clean_id(x):
    return x.split(".")[0]


bg = de["unadjusted"].gene_id.map(clean_id).drop_duplicates()
bg.to_csv(os.path.join(KIT, "panther", "REFERENCE_expressed_genes_19622_EnsemblIDs.txt"), index=False, header=False)
counts = []
for m, d in de.items():
    d = d.dropna(subset=["pvalue"]).copy()
    d["id"] = d.gene_id.map(clean_id)
    for thr in [0.01, 0.05]:
        up = d[(d.pvalue < thr) & (d.log2FoldChange > 0)].id
        dn = d[(d.pvalue < thr) & (d.log2FoldChange < 0)].id
        tag = f"p{str(thr)[2:]}"
        up.to_csv(os.path.join(KIT, "panther", f"{m}_UP_{tag}_EnsemblIDs.txt"), index=False, header=False)
        dn.to_csv(os.path.join(KIT, "panther", f"{m}_DOWN_{tag}_EnsemblIDs.txt"), index=False, header=False)
        counts.append({"model": m, "cutoff": f"p<{thr}", "n_up": len(up), "n_down": len(dn)})
    # ranked list for PANTHER statistical enrichment test: ID <tab> DESeq2 Wald statistic (sign = direction)
    d[["id", "stat"]].dropna().to_csv(os.path.join(KIT, "panther", f"{m}_RANKED_wald_stat.txt"), sep="\t", index=False, header=False)
pd.DataFrame(counts).to_csv(os.path.join(KIT, "panther", "list_sizes.csv"), index=False)

u, f = de["unadjusted"].set_index("gene_id"), de["fibro_adj"].set_index("gene_id")
ipa = pd.DataFrame({"Ensembl_ID": [clean_id(i) for i in u.index], "Symbol": u.gene_name.values,
                    "log2FC_unadjusted": u.log2FoldChange.values, "pvalue_unadjusted": u.pvalue.values, "padj_unadjusted": u.padj.values,
                    "log2FC_fibroblast_adjusted": f.reindex(u.index).log2FoldChange.values, "pvalue_fibroblast_adjusted": f.reindex(u.index).pvalue.values,
                    "padj_fibroblast_adjusted": f.reindex(u.index).padj.values, "baseMean": u.baseMean.values})
ipa = ipa.dropna(subset=["pvalue_unadjusted", "pvalue_fibroblast_adjusted"], how="all")
ipa.to_csv(os.path.join(KIT, "ipa", "IPA_upload_SPP1_vs_NTC_5xFAD_dura.txt"), sep="\t", index=False, float_format="%.6g")
print(pd.DataFrame(counts).to_string())

# ---------------------------------------------------------------- 2. eosinophil sanity check (marrow carry-over?)
long = pd.read_csv(os.path.join(TAB, "spp1_pathway_sets_long.tsv"), sep="\t")
eos_ids = long[long.set == "M_EOSINOPHIL_GRANULE"].gene_id
ifn_ids = long[long.set == "M_TYPE_I_IFN_EXTENDED"].gene_id
z = vst.sub(vst.mean(1), axis=0).div(vst.std(1), axis=0)
eos = z.loc[eos_ids, F5].mean()
ifn = z.loc[ifn_ids, F5].mean()
chk = []
for s in ["neutrophil", "erythroid_blood", "bone_osteoblast", "osteoclast", "macrophage_BAM", "monocyte", "B_cell", "fibroblast_stromal", "olfactory_nasal"]:
    chk.append({"signature": s, "rho_with_eosinophil_score": spearmanr(eos, S.loc[F5, s])[0], "rho_with_IFN_score": spearmanr(ifn, S.loc[F5, s])[0]})
chk = pd.DataFrame(chk)
chk.loc[len(chk)] = {"signature": "IFN score", "rho_with_eosinophil_score": spearmanr(eos, ifn)[0], "rho_with_IFN_score": 1.0}
chk.to_csv(os.path.join(OUT, "eosinophil_and_IFN_vs_composition.csv"), index=False)
print(chk.round(2).to_string())
pd.DataFrame({"sample": meta.loc[F5, "sample"], "group": meta.loc[F5, "group"], "sex": meta.loc[F5, "sex"],
              "eosinophil_score": eos, "IFN_score": ifn}).to_csv(os.path.join(OUT, "eosinophil_IFN_scores_per_animal.csv"))

# ---------------------------------------------------------------- 3. mechanistic gene table
MECH_GENES = {"SPP1 receptors / integrin-FAK": ["Spp1", "Cd44", "Itgav", "Itgb1", "Itgb3", "Itgb5", "Ptk2", "Ilk"],
              "YAP/TAZ (Hippo)": ["Yap1", "Wwtr1", "Tead1", "Amot", "Sav1", "Lats2", "Ccn1", "Ccn2", "Ankrd1", "Amotl2"],
              "TGF-beta / SMAD": ["Tgfb1", "Tgfb2", "Smad3", "Smad7", "Skil", "Bmp2", "Serpine1", "Bhlhe40"],
              "Matrix output": ["Fn1", "Col3a1", "Col12a1", "Lox", "Loxl3", "Fbln1", "Timp1", "Cthrc1", "Prg4"],
              "Barrier": ["Cldn11", "Tjp1", "Cdh1"],
              "Interferon / JAK-STAT": ["Rsad2", "Isg20", "Ifit3b", "Ifit1", "Isg15", "Irf7", "Stat1", "Irf9", "Cxcl10"],
              "Cholesterol / SREBP": ["Insig1", "Sqle", "Idi1", "Msmo1", "Hmgcr", "Ldlr"],
              "Eosinophil": ["Prg2", "Prg3", "Epx", "Ear2", "Alox15"]}
rows = []
for mod, gl in MECH_GENES.items():
    for g in gl:
        r_u = u[u.gene_name == g]
        r_f = f[f.gene_name == g]
        if not len(r_u):
            continue
        rows.append({"module": mod, "gene": g, "baseMean": r_u.baseMean.iloc[0], "log2FC_unadjusted": r_u.log2FoldChange.iloc[0],
                     "p_unadjusted": r_u.pvalue.iloc[0], "log2FC_fibro_adj": r_f.log2FoldChange.iloc[0] if len(r_f) else np.nan,
                     "p_fibro_adj": r_f.pvalue.iloc[0] if len(r_f) else np.nan})
mt = pd.DataFrame(rows)
mt.to_csv(os.path.join(OUT, "mechanistic_axis_genes.csv"), index=False)

# ---------------------------------------------------------------- 4. link map
camw = cam.pivot(index="set", columns="model", values="p_conservative")
dirn = cam.pivot(index="set", columns="model", values="direction")
lead = cam.pivot(index="set", columns="model", values="leading_genes").fillna("")
keep = camw.index[(camw.min(axis=1) < 0.15) & ~camw.index.str.startswith("TECHNICAL")]
lead_sets = {s: set(";".join(lead.loc[s].values).split(";")) - {""} for s in keep}
G = nx.Graph()
for s in keep:
    best = camw.loc[s].idxmin()
    G.add_node(s, p=camw.loc[s].min(), model=best, direction=dirn.loc[s, best],
               fibro_robust=bool(camw.loc[s, "fibro_adj"] < 0.15 and dirn.loc[s, "fibro_adj"] == dirn.loc[s, "unadjusted"]))
for a, b in itertools.combinations(keep, 2):
    A, B = lead_sets[a], lead_sets[b]
    if not A or not B:
        continue
    inter = A & B
    jac = len(inter) / len(A | B)
    if jac >= 0.1 and len(inter) >= 2:
        G.add_edge(a, b, weight=jac, shared=";".join(sorted(inter)))
edges = pd.DataFrame([{"set_a": a, "set_b": b, "jaccard": d["weight"], "shared_leading_genes": d["shared"]} for a, b, d in G.edges(data=True)])
edges.sort_values("jaccard", ascending=False).to_csv(os.path.join(OUT, "pathway_link_edges.csv"), index=False)
nodes = pd.DataFrame([{"set": n, **d} for n, d in G.nodes(data=True)])
nodes["cluster"] = -1
for i, comp in enumerate(sorted(nx.connected_components(G), key=len, reverse=True)):
    nodes.loc[nodes.set.isin(comp), "cluster"] = i
nodes.to_csv(os.path.join(OUT, "pathway_link_nodes.csv"), index=False)

pos = nx.spring_layout(G, weight="weight", seed=7, k=1.1)
fig, ax = plt.subplots(figsize=(15, 11))
for a, b, d in G.edges(data=True):
    ax.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]], color="#c9c8c3", lw=0.6 + 5 * d["weight"], zorder=1)
for n, d in G.nodes(data=True):
    col = "#e34948" if d["direction"] == "Up" else "#2a78d6"
    size = 120 + 900 * min(-np.log10(d["p"]), 3) / 3
    ax.scatter(*pos[n], s=size, c=col, edgecolor=INK if d["fibro_robust"] else SURFACE, lw=2.2 if d["fibro_robust"] else 1, zorder=3)
    import re
    lab = re.sub(r"^H_", "H: ", re.sub(r"^C_", "", re.sub(r"^M_", "★ ", n)))
    ax.annotate(lab, pos[n], xytext=(6, 5), textcoords="offset points", fontsize=8)
ax.axis("off")
ax.set_title("SPP1 vs NTC (5xFAD dura): how pathways connect through shared driver genes\n"
             "red = up, blue = down with SPP1; size = -log10 conservative camera p; thick black ring = holds after fibroblast adjustment;\n"
             "edges = shared leading genes (Jaccard); ★ = mechanistic module", fontsize=10)
save(fig, os.path.join(OUT, "SPP1_pathway_link_map.png"))

# mechanistic heatmap (LFC both models)
fig, ax = plt.subplots(figsize=(7, 14))
M = mt.set_index("gene")[["log2FC_unadjusted", "log2FC_fibro_adj"]]
P = mt.set_index("gene")[["p_unadjusted", "p_fibro_adj"]]
im = ax.imshow(M.values, cmap=DIV, vmin=-1, vmax=1, aspect="auto")
ax.set_yticks(range(len(M)), [f"{g}  [{m}]" for g, m in zip(mt.gene, mt.module)], fontsize=7.5)
ax.set_xticks([0, 1], ["unadjusted", "fibroblast-adjusted"])
for i in range(len(M)):
    for j in range(2):
        p = P.values[i, j]
        if not np.isnan(M.values[i, j]):
            ax.text(j, i, f"{M.values[i, j]:+.2f}{'**' if p < 0.01 else ('*' if p < 0.05 else '')}", ha="center", va="center", fontsize=7,
                    color="white" if abs(M.values[i, j]) > 0.7 else INK)
ax.grid(False)
fig.colorbar(im, ax=ax, shrink=0.4, label="log2FC SPP1 vs NTC (* p<0.05, ** p<0.01)")
ax.set_title("Mechanistic axis genes")
save(fig, os.path.join(OUT, "SPP1_mechanistic_axis_genes.png"))
print(nodes.sort_values(["cluster", "p"]).to_string())
print(edges.sort_values("jaccard", ascending=False).head(25).to_string())
