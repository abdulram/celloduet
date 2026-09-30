"""Build a mouse cell-type pseudobulk reference from the CELLxGENE Census (release 2025-01-30).

For each reference class, up to N cells are sampled (capped per dataset so no single study dominates),
raw counts are normalised to counts-per-10k per cell, averaged per dataset, and then averaged across datasets.
Output: data/reference/census_pseudobulk_cp10k.csv (genes x classes) + class/cell provenance tables.
"""
import os

import numpy as np
import pandas as pd
import scipy.sparse as sp
import tiledbsoma as soma

from census_ctx import open_census

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "data", "reference")
os.makedirs(OUT, exist_ok=True)
N_PER_CLASS, N_PER_DATASET, SEED = 1200, 300, 7

DURA = ["dura mater", "meningeal dura mater"]
# class -> (cell_type list, optional tissue list, is_fibroblast)
CLASSES = {
    "FIB_dura": (["stromal cell", "fibroblast"], DURA, True),
    "FIB_leptomeningeal": (["vascular leptomeningeal cell", "leptomeningeal cell", "vascular leptomeningeal cell (Mmus)"], None, True),
    "FIB_brain_perivascular": (["fibroblast"], "BRAIN", True),
    "osteoblast": (["osteoblast"], None, False),
    "osteoclast": (["osteoclast"], None, False),
    "chondrocyte": (["chondrocyte"], None, False),
    "choroid_plexus_epithelial": (["choroid plexus epithelial cell"], None, False),
    "ependymal": (["ependymal cell"], None, False),
    "olfactory_receptor_neuron": (["olfactory receptor cell"], None, False),
    "olfactory_epithelial": (["olfactory epithelial cell"], None, False),
    "olfactory_ensheathing": (["olfactory ensheathing cell"], None, False),
    "erythrocyte": (["erythrocyte", "erythroblast", "proerythroblast"], None, False),
    "neutrophil": (["neutrophil", "granulocyte"], None, False),
    "monocyte": (["monocyte"], None, False),
    "macrophage": (["macrophage", "meningeal macrophage", "perivascular macrophage"], None, False),
    "microglia": (["microglial cell"], None, False),
    "B_cell": (["B cell"], None, False),
    "plasma_cell": (["plasma cell"], None, False),
    "T_cell": (["T cell", "CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell"], None, False),
    "NK_ILC": (["natural killer cell", "innate lymphoid cell"], None, False),
    "dendritic_cell": (["conventional dendritic cell", "dendritic cell"], None, False),
    "mast_cell": (["mast cell"], None, False),
    "endothelial": (["endothelial cell"], None, False),
    "pericyte": (["pericyte"], None, False),
    "smooth_muscle": (["smooth muscle cell", "vascular associated smooth muscle cell"], None, False),
    "skeletal_muscle": (["skeletal muscle satellite cell", "skeletal muscle myoblast", "skeletal muscle fiber"], None, False),
    "neuron": (["glutamatergic neuron", "GABAergic neuron", "neuron"], "BRAIN", False),
    "astrocyte": (["astrocyte"], None, False),
    "oligodendrocyte": (["oligodendrocyte"], None, False),
    "OPC": (["oligodendrocyte precursor cell"], None, False),
    "schwann": (["Schwann cell"], None, False),
    "keratinocyte": (["keratinocyte"], None, False),
    "adipocyte": (["adipocyte"], None, False),
}

# Present in the Census only (or almost only) in the embryonic organogenesis atlas
EMBRYO_OK = {"osteoblast", "osteoclast", "choroid_plexus_epithelial", "olfactory_receptor_neuron",
             "olfactory_epithelial", "olfactory_ensheathing"}
with open_census() as census:
    mm = census["census_data"]["mus_musculus"]
    obs = mm.obs.read(value_filter="is_primary_data == True",
                      column_names=["soma_joinid", "cell_type", "tissue_general", "tissue", "dataset_id", "assay"]).concat().to_pandas()
    is_embryo = obs.tissue_general.astype(str) == "embryo"
    var = mm.ms["RNA"].var.read(column_names=["soma_joinid", "feature_id", "feature_name"]).concat().to_pandas()

    picked = []
    for cls, (cts, tissues, is_fib) in CLASSES.items():
        # adult cells only, except classes that exist in the Census only in the embryonic atlas
        sub = obs[obs.cell_type.isin(cts) & (~is_embryo | (cls in EMBRYO_OK))]
        if tissues == "BRAIN":
            sub = sub[sub.tissue_general.astype(str).isin(["brain", "central nervous system", "nervous system", "spinal cord"])]
            sub = sub[~sub.tissue.astype(str).isin(DURA)]
        elif tissues is not None:
            sub = sub[sub.tissue.astype(str).isin(tissues)]
        # cap per dataset (lifted when a class comes from a single study), then per class
        cap = N_PER_CLASS if sub.dataset_id.nunique() == 1 else N_PER_DATASET
        parts = []
        for ds, g in sub.groupby("dataset_id", observed=True):
            parts.append(g.sample(min(len(g), cap), random_state=SEED))
        sub = pd.concat(parts) if parts else sub
        if len(sub) > N_PER_CLASS:
            sub = sub.sample(N_PER_CLASS, random_state=SEED)
        sub = sub.assign(ref_class=cls, is_fibroblast=is_fib)
        picked.append(sub)
        print(f"{cls:28s} cells={len(sub):5d} datasets={sub.dataset_id.nunique()}", flush=True)
    picked = pd.concat(picked)
    picked.to_csv(os.path.join(OUT, "census_cells_used.csv.gz"), index=False)

    ids = np.sort(picked.soma_joinid.unique())
    with mm.axis_query("RNA", obs_query=soma.AxisQuery(coords=(ids,))) as q:
        X = q.X("raw").tables()
        rows, cols, vals = [], [], []
        for t in X:
            rows.append(t["soma_dim_0"].to_numpy())
            cols.append(t["soma_dim_1"].to_numpy())
            vals.append(t["soma_data"].to_numpy())
    rows, cols, vals = np.concatenate(rows), np.concatenate(cols), np.concatenate(vals).astype(np.float32)

row_pos = pd.Series(np.arange(len(ids)), index=ids)
M = sp.csr_matrix((vals, (row_pos.loc[rows].values, cols)), shape=(len(ids), int(var.soma_joinid.max()) + 1))
lib = np.asarray(M.sum(axis=1)).ravel()
M = sp.diags(1e4 / np.maximum(lib, 1)) @ M  # CP10K per cell

picked = picked.drop_duplicates("soma_joinid").set_index("soma_joinid").loc[ids]
pb = {}
for cls, g in picked.groupby("ref_class", sort=False):
    per_ds = []
    for ds, gd in g.groupby("dataset_id", observed=True):
        per_ds.append(np.asarray(M[row_pos.loc[gd.index].values].mean(axis=0)).ravel())
    pb[cls] = np.mean(per_ds, axis=0)
pb = pd.DataFrame(pb)
pb.index = var.set_index("soma_joinid").reindex(pb.index)["feature_id"].values
pb = pb[pb.index.notna()]
pb.to_csv(os.path.join(OUT, "census_pseudobulk_cp10k.csv.gz"))
var.to_csv(os.path.join(OUT, "census_genes.csv.gz"), index=False)
summary = picked.groupby("ref_class").agg(cells=("dataset_id", "size"), datasets=("dataset_id", "nunique"),
                                          is_fibroblast=("is_fibroblast", "first"))
summary.to_csv(os.path.join(OUT, "census_class_summary.csv"))
print(summary)
