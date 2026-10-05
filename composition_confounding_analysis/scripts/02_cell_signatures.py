"""02 - Cellular / anatomical composition: marker signatures (primary) + reference deconvolution (secondary).

A. Marker signature scores (primary composition measure; transparent, no reference assumptions)
   - curated broad-lineage marker sets (documented below; chosen to avoid markers shared across lineages)
   - per gene z-score of blind VST across samples, signature = mean z of detected markers
   - VALIDATION: markers of one population must move together. For each signature we report the mean pairwise
     Pearson r and each marker's item-rest correlation. A signature is "validated" if >= 3 markers are detected
     (median CPM >= 1) and the median item-rest correlation is >= 0.4. Only validated signatures enter the
     composition matrix. Non-validated signatures are reported but not used.
B. Reference-based deconvolution (secondary; RNA-contribution estimates, NOT cell counts)
   - dura-centric CELLxGENE pseudobulk reference (00_build_dura_reference.py)
   - trimmed NNLS ensemble (4 settings; markers >= 3-5x specific; mito/ribosomal/Ig-V genes excluded; markers the
     fit misses > 4-8x in a sample are dropped and refit). Stability = agreement across settings.
Outputs: 02_cell_composition/, tables/composition_signature_scores.tsv, tables/deconvolution_*.tsv
"""
import itertools
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.optimize import nnls
from scipy.stats import spearmanr

np.random.seed(7)
HERE = os.path.dirname(os.path.abspath(__file__))
CA = os.path.dirname(HERE)
ROOT = os.path.dirname(CA)
OUT = os.path.join(CA, "02_cell_composition")
TAB = os.path.join(CA, "tables")
exec(open(os.path.join(HERE, "common.py")).read())

meta = pd.read_csv(os.path.join(TAB, "sample_metadata_qc_pcs.tsv"), sep="\t", index_col=0)
vst = pd.read_csv(os.path.join(CA, "01_QC", "vst_blind.tsv.gz"), sep="\t", index_col=0)
raw = pd.read_csv(os.path.join(ROOT, "data", "VT3YXB-expression-matrix.tsv"), sep="\t")
raw = raw[~raw.gene_id.str.startswith("ERCC")].set_index("gene_id")
raw["gene_name"] = raw.gene_name.fillna(pd.Series(raw.index, index=raw.index))
libs = list(meta.index)
cpm = raw[[f"{l}_cpm" for l in libs]]
cpm.columns = libs
cnt = raw[[f"{l}_count" for l in libs]]
cnt.columns = libs
name2id = pd.Series(raw.index, index=raw.gene_name.values)
name2id = name2id[~name2id.index.duplicated()]

# ---------------------------------------------------------------- A. marker signatures
# Native dura lineages (spec) followed by adjacent / contaminating tissues identified earlier in this project.
SIGNATURES = {
    "fibroblast_stromal": ["Pdgfra", "Dcn", "Fbln1", "Mfap4", "Igfbp6", "Slc47a1", "Col14a1", "Smoc2"],
    "blood_endothelial": ["Pecam1", "Cdh5", "Cldn5", "Flt1", "Kdr", "Emcn", "Tek", "Esam"],
    "lymphatic_endothelial": ["Prox1", "Flt4", "Ccl21a", "Mmrn1", "Reln", "Tbx1"],
    "mural_pericyte_SMC": ["Myh11", "Cnn1", "Rgs5", "Kcnj8", "Abcc9", "Acta2", "Tagln"],
    "macrophage_BAM": ["Mrc1", "Cd163", "F13a1", "Pf4", "C1qa", "C1qb", "Adgre1", "Ms4a7"],
    "monocyte": ["Ly6c2", "Ccr2", "Plac8", "Ms4a4c", "F10"],
    "dendritic_cell": ["Flt3", "Xcr1", "Clec9a", "Cd209a", "Itgax", "Siglech"],
    "T_cell": ["Cd3e", "Cd3g", "Cd3d", "Cd2", "Cd5", "Cd28", "Lck"],
    "B_cell": ["Cd79a", "Cd79b", "Ms4a1", "Cd19", "Pax5", "Fcmr"],
    "plasma_cell": ["Jchain", "Mzb1", "Prdm1", "Tnfrsf17", "Sdc1"],
    "NK_cell": ["Nkg7", "Klrb1c", "Ncr1", "Gzma", "Klrd1", "Prf1"],
    "mast_cell": ["Cpa3", "Mcpt4", "Cma1", "Tpsb2", "Kit", "Fcer1a"],
    "neutrophil": ["S100a8", "S100a9", "Ngp", "Camp", "Ltf", "Retnlg", "Chil3"],
    "CNS_neural": ["Snap25", "Syt1", "Stmn2", "Plp1", "Mobp", "Aqp4"],
    "bone_osteoblast": ["Bglap", "Bglap2", "Ibsp", "Sp7", "Dmp1", "Sost"],
    "osteoclast": ["Ctsk", "Acp5", "Mmp9", "Atp6v0d2", "Dcstamp", "Ocstamp"],
    "olfactory_nasal": ["Omp", "Cyp2g1", "Ugt2a1", "Gnal", "Scgb1c1", "Bpifa1"],
    "choroid_plexus": ["Ttr", "Kl", "Folr1", "Kcnj13", "Clic6"],
    "erythroid_blood": ["Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Alas2"],
    "skeletal_muscle": ["Acta1", "Ckm", "Myh1", "Tnnt3", "Mylpf"],
}
NATIVE = ["fibroblast_stromal", "blood_endothelial", "lymphatic_endothelial", "mural_pericyte_SMC", "macrophage_BAM",
          "monocyte", "dendritic_cell", "T_cell", "B_cell", "plasma_cell", "NK_cell", "mast_cell", "neutrophil"]

zv = vst.sub(vst.mean(axis=1), axis=0).div(vst.std(axis=1).replace(0, np.nan), axis=0)
scores, val_rows, item_rows = {}, [], []
for sig, gl in SIGNATURES.items():
    ids = [name2id[g] for g in gl if g in name2id.index and name2id[g] in vst.index]
    det = [i for i in ids if cpm.loc[i].median() >= 1]
    use = det if len(det) >= 2 else ids
    Z = zv.loc[use]
    scores[sig] = Z.mean()
    rr = []
    for i in use:
        rest = Z.drop(i).mean()
        r = np.corrcoef(Z.loc[i], rest)[0, 1] if len(use) > 1 else np.nan
        rr.append(r)
        item_rows.append({"signature": sig, "gene": raw.gene_name[i], "median_CPM": round(cpm.loc[i].median(), 2), "item_rest_r": r})
    C = np.corrcoef(Z.values) if len(use) > 1 else np.array([[1.0]])
    pair = C[np.triu_indices_from(C, 1)] if len(use) > 1 else np.array([np.nan])
    valid = (len(det) >= 3) and (np.nanmedian(rr) >= 0.4)
    val_rows.append({"signature": sig, "category": "native dura" if sig in NATIVE else "adjacent/contaminant",
                     "markers_defined": len(gl), "markers_detected": len(det), "mean_pairwise_r": np.nanmean(pair),
                     "median_item_rest_r": np.nanmedian(rr), "validated": valid, "markers_used": ",".join(raw.gene_name[use])})
scores = pd.DataFrame(scores)
val = pd.DataFrame(val_rows)
val.to_csv(os.path.join(TAB, "signature_validation.tsv"), sep="\t", index=False)
pd.DataFrame(item_rows).to_csv(os.path.join(TAB, "signature_marker_item_rest.tsv"), sep="\t", index=False)
scores.to_csv(os.path.join(TAB, "composition_signature_scores_all.tsv"), sep="\t")
validated = val.signature[val.validated].tolist()
scores[validated].to_csv(os.path.join(TAB, "composition_signature_scores.tsv"), sep="\t")
print(val[["signature", "markers_detected", "mean_pairwise_r", "median_item_rest_r", "validated"]].round(2).to_string())

# ---------------------------------------------------------------- B. reference deconvolution
REF = os.path.join(CA, "reference", "dura_reference_pseudobulk_cp10k.csv.gz")
ref = pd.read_csv(REF, index_col=0).reindex(raw.index).dropna(how="all").fillna(0)
expressed = ((cnt >= 10).sum(axis=1) >= 6) & raw.gene_biotype.isin(["protein_coding", "lncRNA"])
gname = raw.gene_name.reindex(ref.index).astype(str)
bad = gname.str.match(r"^(mt-|Rp[ls]\d|Mrp[ls]\d|Igkv|Iglv|Ighv|Trav|Trbv|Trgv|Trdv)") | raw.gene_biotype.reindex(ref.index).astype(str).str.match(r"^(IG|TR|Mt)_")
ref_e = ref.loc[ref.index.intersection(expressed.index[expressed]).difference(ref.index[bad])]
EPS = 0.02


def build_markers(ntop, minratio):
    out = {}
    for c in ref.columns:
        others = ref_e.drop(columns=c).max(axis=1)
        r = (ref_e[c] + EPS) / (others + EPS)
        m = r[(ref_e[c] >= 1) & (r >= minratio)].sort_values(ascending=False).head(ntop).index
        if len(m) >= 5:
            out[c] = m
    return out


def trimmed_nnls(mk, thr, iters):
    genes = pd.Index(sorted(set().union(*mk.values())))
    S = ref.loc[genes, list(mk)] * 100
    B = cpm.loc[genes]
    w = 1 / np.sqrt(S.mean(axis=1).values + 1)
    coefs, fq = {}, {}
    for l in libs:
        keep = np.ones(len(genes), bool)
        for it in range(iters):
            coef, _ = nnls(S.values[keep] * w[keep, None], B[l].values[keep] * w[keep])
            pred = S.values @ coef
            miss = np.abs(np.log2((B[l].values + 1) / (pred + 1))) > thr
            if it < iters - 1:
                keep = ~miss
        coefs[l] = pd.Series(coef, index=list(mk))
        fq[l] = spearmanr(np.log2(B[l].values[keep] + 1), np.log2(pred[keep] + 1))[0]
    P = pd.DataFrame(coefs).T
    return P.div(P.sum(axis=1), axis=0), pd.Series(fq)


SETTINGS = [(25, 3, 2, 4), (25, 3, 1.5, 5), (50, 5, 2, 4), (50, 5, 1.5, 5)]
fits, fqs = [], []
for ntop, mr, thr, it in SETTINGS:
    P, fq = trimmed_nnls(build_markers(ntop, mr), thr, it)
    fits.append(P.reindex(columns=ref.columns).fillna(0))
    fqs.append(fq)
deconv = sum(fits) / len(fits)
deconv.to_csv(os.path.join(TAB, "deconvolution_proportions.tsv"), sep="\t")
# stability: per class, mean pairwise Spearman correlation of estimates across the 4 settings
stab = {}
for c in ref.columns:
    rs = [spearmanr(a[c], b[c])[0] for a, b in itertools.combinations(fits, 2) if a[c].std() > 0 and b[c].std() > 0]
    stab[c] = np.nanmean(rs) if rs else np.nan
MAP = {"fibroblast_stromal": "fibroblast_stromal", "blood_endothelial": "blood_endothelial", "lymphatic_endothelial": "lymphatic_endothelial",
       "mural_pericyte_SMC": "mural_pericyte_SMC", "macrophage_BAM": "macrophage_BAM", "monocyte": "monocyte", "dendritic_cell": "dendritic_cell",
       "T_cell": "T_cell", "B_cell": "B_cell", "plasma_cell": "plasma_cell", "NK_ILC": "NK_cell", "mast_cell": "mast_cell",
       "neutrophil": "neutrophil", "CNS_neural_glial": "CNS_neural", "osteoblast": "bone_osteoblast", "osteoclast": "osteoclast",
       "olfactory_epithelium": "olfactory_nasal", "choroid_plexus": "choroid_plexus", "erythroid": "erythroid_blood", "skeletal_muscle": "skeletal_muscle"}
dv = pd.DataFrame({"reference_class": ref.columns,
                   "median_pct": (deconv.median() * 100).round(1).values,
                   "max_pct": (deconv.max() * 100).round(1).values,
                   "stability_across_settings_rho": [stab[c] for c in ref.columns],
                   "agreement_with_marker_signature_rho": [spearmanr(deconv[c], scores[MAP[c]])[0] if deconv[c].std() > 0 else np.nan for c in ref.columns]})
dv["reliable"] = (dv.stability_across_settings_rho >= 0.7) & (dv.agreement_with_marker_signature_rho >= 0.5)
dv.to_csv(os.path.join(TAB, "deconvolution_validation.tsv"), sep="\t", index=False)
pd.DataFrame({"fit_quality_rho": sum(fqs) / len(fqs)}).to_csv(os.path.join(TAB, "deconvolution_fit_quality.tsv"), sep="\t")
print(dv.round(2).to_string())

# ---------------------------------------------------------------- figures
# provisional ordering by first PC of validated signatures (formal composition PCA is in step 03)
Zs = scores[validated]
Zs = (Zs - Zs.mean()) / Zs.std()
u, s, vt = np.linalg.svd(Zs.values, full_matrices=False)
pc1 = pd.Series(u[:, 0] * s[0], index=libs)
if np.corrcoef(pc1, scores.get("fibroblast_stromal", pc1))[0, 1] < 0:
    pc1 = -pc1
order = pc1.sort_values().index

# Figure 2: signature heatmap (all signatures, validated marked)
fig, ax = plt.subplots(figsize=(13, 7))
M = scores.loc[order].T
im = ax.imshow(M.values, cmap=DIV, vmin=-2, vmax=2, aspect="auto")
ylab = [f"{s}{'' if s in validated else '  (not validated)'}" for s in M.index]
ax.set_yticks(range(len(M)), ylab, fontsize=8)
for t, s in zip(ax.get_yticklabels(), M.index):
    if s not in validated:
        t.set_color("#9a9994")
ax.set_xticks(range(len(order)), [f"{meta.animal[l]} {meta.group[l]} {meta.sex[l]}" for l in order], rotation=90, fontsize=7.5)
for t, l in zip(ax.get_xticklabels(), order):
    t.set_color(GROUP_COLORS[meta.group[l]])
ax.axhline(len(NATIVE) - 0.5, color=INK, lw=1)
ax.grid(False)
fig.colorbar(im, ax=ax, shrink=0.6, label="signature score (mean marker z, VST)")
ax.set_title("Broad cellular / anatomical signatures per animal (ordered by provisional composition axis)")
save(fig, os.path.join(OUT, "signature_heatmap.png"))

# Marker heatmap: every marker, grouped by signature, samples ordered by composition axis
rows, labels, seps = [], [], []
for sig in SIGNATURES:
    ids = [name2id[g] for g in SIGNATURES[sig] if g in name2id.index and name2id[g] in vst.index]
    rows += ids
    labels += [f"{raw.gene_name[i]}" for i in ids]
    seps.append(len(rows))
fig, ax = plt.subplots(figsize=(13, 26))
im = ax.imshow(zv.loc[rows, order].values, cmap=DIV, vmin=-2.5, vmax=2.5, aspect="auto")
ax.set_yticks(range(len(rows)), labels, fontsize=6.5)
start = 0
for sig, end in zip(SIGNATURES, seps):
    ax.text(len(order) + 0.3, (start + end - 1) / 2, sig + ("" if sig in validated else " (n.v.)"), va="center", fontsize=8,
            color=INK if sig in validated else "#9a9994")
    ax.axhline(end - 0.5, color=SURFACE, lw=2)
    start = end
ax.set_xticks(range(len(order)), [f"{meta.animal[l]} {meta.group[l]}" for l in order], rotation=90, fontsize=7)
for t, l in zip(ax.get_xticklabels(), order):
    t.set_color(GROUP_COLORS[meta.group[l]])
ax.grid(False)
ax.set_title("Lineage marker genes (row z-score of VST); samples ordered by provisional composition axis")
fig.colorbar(im, ax=ax, shrink=0.2, label="z")
save(fig, os.path.join(OUT, "lineage_marker_heatmap.png"))

# Figure 3: signature scores by genotype/group, individual animals
sig_plot = validated
ncol = 4
nrow = int(np.ceil(len(sig_plot) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(15, 3.2 * nrow), sharex=True)
for ax, s in zip(axes.flat, sig_plot):
    strip(ax, scores[s], meta)
    ax.set_title(s, fontsize=9.5)
for ax in list(axes.flat)[len(sig_plot):]:
    ax.axis("off")
fig.suptitle("Validated composition signatures by group (each point = one animal; circle F, triangle M)", fontweight="bold")
fig.tight_layout()
save(fig, os.path.join(OUT, "signature_scores_by_group.png"))

# Deconvolution stacked bars
grp = {"Fibroblast/stromal": ["fibroblast_stromal"], "Blood endothelial": ["blood_endothelial"], "Lymphatic endothelial": ["lymphatic_endothelial"],
       "Mural": ["mural_pericyte_SMC"], "Myeloid (BAM, mono, DC, neutrophil, mast)": ["macrophage_BAM", "monocyte", "dendritic_cell", "neutrophil", "mast_cell"],
       "Lymphoid (T, B, plasma, NK)": ["T_cell", "B_cell", "plasma_cell", "NK_ILC"], "Bone (osteoblast, osteoclast)": ["osteoblast", "osteoclast"],
       "Other adjacent (CNS, CP, olfactory, blood, muscle)": ["CNS_neural_glial", "choroid_plexus", "olfactory_epithelium", "erythroid", "skeletal_muscle"]}
gd = pd.DataFrame({k: deconv[v].sum(axis=1) for k, v in grp.items()})
gd.to_csv(os.path.join(TAB, "deconvolution_grouped.tsv"), sep="\t")
cols = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
fig, ax = plt.subplots(figsize=(12, 7.5))
for i, l in enumerate(order):
    left = 0
    for j, c in enumerate(gd.columns):
        v = gd.loc[l, c] * 100
        ax.barh(i, v, left=left, color=cols[j], height=0.8, edgecolor=SURFACE, lw=2)
        if v >= 8:
            ax.text(left + v / 2, i, f"{v:.0f}", ha="center", va="center", fontsize=7.5, color="white" if j in (0, 5, 6) else INK)
        left += v
ax.set_yticks(range(len(order)), [f"{meta.animal[l]} {meta.group[l]} {meta.sex[l]}" for l in order], fontsize=8)
for t, l in zip(ax.get_yticklabels(), order):
    t.set_color(GROUP_COLORS[meta.group[l]])
ax.set_xlim(0, 100)
ax.invert_yaxis()
ax.grid(axis="y", visible=False)
ax.set_xlabel("estimated share of marker-gene RNA signal (%) - NOT cell counts")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=c, label=n) for n, c in zip(gd.columns, cols)], ncol=2, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.08))
ax.set_title("Reference deconvolution (dura-centric CELLxGENE reference, trimmed NNLS ensemble)")
save(fig, os.path.join(OUT, "deconvolution_stacked.png"))
print("validated signatures:", validated)
