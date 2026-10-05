"""00 - Build a dura-centric single-cell pseudobulk reference from the CELLxGENE Census (release 2025-01-30).

Immune / stromal / endothelial classes come from the adult mouse DURA dataset in the Census
(dataset 58b01044-c5e5-4b0f-8a2d-6ebf951e01ff, tissues 'dura mater' / 'meningeal dura mater'), so that their
profiles reflect dura cells. Classes not present in that dataset (lymphatic endothelium, mural cells) and
the contaminating tissues identified in this project (bone, osteoclast, olfactory, choroid plexus, erythrocyte,
skeletal muscle, CNS) are taken from other adult (or, where only embryonic data exist, embryonic) Census data.

Per class: <= 1200 cells, capped per dataset; per-cell CP10K; mean per dataset; mean across datasets.
Output: composition_confounding_analysis/reference/
"""
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp
import tiledbsoma as soma

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "analysis", "reference"))
from census_ctx import open_census  # noqa: E402

OUT = os.path.join(ROOT, "composition_confounding_analysis", "reference")
os.makedirs(OUT, exist_ok=True)
SEED, N_CLASS, N_DS = 7, 1200, 300
DURA_DS = "58b01044-c5e5-4b0f-8a2d-6ebf951e01ff"
DURA_T = ["dura mater", "meningeal dura mater"]
BRAIN = ["brain", "central nervous system", "nervous system", "spinal cord"]

# class: (cell types, source) ; source = "dura" | "adult" | "adult_brain" | "any" (embryo allowed)
CLASSES = {
    # --- native dura (from the dura dataset)
    "fibroblast_stromal": (["stromal cell", "fibroblast"], "dura"),
    "blood_endothelial": (["endothelial cell"], "dura"),
    "macrophage_BAM": (["macrophage"], "dura"),
    "monocyte": (["monocyte", "classical monocyte", "non-classical monocyte"], "dura"),
    "dendritic_cell": (["conventional dendritic cell", "plasmacytoid dendritic cell"], "dura"),
    "neutrophil": (["neutrophil"], "dura"),
    "B_cell": (["B cell"], "dura"),
    "plasma_cell": (["plasma cell"], "dura"),
    "T_cell": (["CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell", "T cell",
                "gamma-delta T cell"], "dura"),
    "NK_ILC": (["natural killer cell", "innate lymphoid cell"], "dura"),
    "mast_cell": (["mast cell"], "dura"),
    # --- native but absent from the dura dataset
    "lymphatic_endothelial": (["lymphatic endothelial cell of medulla ceiling", "lymphatic endothelial cell of subcapsular sinus floor",
                               "lymphatic endothelial cell of subcapsular sinus ceiling", "lymph node lymphatic vessel endothelial cell",
                               "endothelial cell of lymphatic vessel", "lymphatic endothelial cell"], "adult"),
    "mural_pericyte_SMC": (["pericyte", "smooth muscle cell", "vascular associated smooth muscle cell"], "adult"),
    # --- contaminating / adjacent tissues
    "osteoblast": (["osteoblast"], "any"),
    "osteoclast": (["osteoclast"], "any"),
    "erythroid": (["erythrocyte", "erythroblast", "proerythroblast"], "adult"),
    "olfactory_epithelium": (["olfactory epithelial cell", "olfactory receptor cell"], "any"),
    "choroid_plexus": (["choroid plexus epithelial cell"], "any"),
    "skeletal_muscle": (["skeletal muscle satellite cell", "skeletal muscle myoblast", "skeletal muscle fiber"], "adult"),
    "CNS_neural_glial": (["glutamatergic neuron", "GABAergic neuron", "neuron", "astrocyte", "oligodendrocyte"], "adult_brain"),
}

with open_census() as census:
    mm = census["census_data"]["mus_musculus"]
    obs = mm.obs.read(value_filter="is_primary_data == True",
                      column_names=["soma_joinid", "cell_type", "tissue_general", "tissue", "dataset_id", "assay"]).concat().to_pandas()
    embryo = obs.tissue_general.astype(str) == "embryo"
    var = mm.ms["RNA"].var.read(column_names=["soma_joinid", "feature_id", "feature_name"]).concat().to_pandas()
    picked = []
    for cls, (cts, src) in CLASSES.items():
        sub = obs[obs.cell_type.isin(cts)]
        if src == "dura":
            sub = sub[(sub.dataset_id == DURA_DS) & sub.tissue.astype(str).isin(DURA_T)]
        elif src == "adult":
            sub = sub[~embryo.loc[sub.index]]
        elif src == "adult_brain":
            sub = sub[~embryo.loc[sub.index] & sub.tissue_general.astype(str).isin(BRAIN)]
        cap = N_CLASS if sub.dataset_id.nunique() == 1 else N_DS
        parts = [g.sample(min(len(g), cap), random_state=SEED) for _, g in sub.groupby("dataset_id", observed=True)]
        sub = pd.concat(parts) if parts else sub
        if len(sub) > N_CLASS:
            sub = sub.sample(N_CLASS, random_state=SEED)
        picked.append(sub.assign(ref_class=cls, source=src))
        print(f"{cls:24s} cells={len(sub):5d} datasets={sub.dataset_id.nunique()}", flush=True)
    picked = pd.concat(picked)
    picked.to_csv(os.path.join(OUT, "reference_cells_used.csv.gz"), index=False)
    ids = np.sort(picked.soma_joinid.unique())
    with mm.axis_query("RNA", obs_query=soma.AxisQuery(coords=(ids,))) as q:
        rows, cols, vals = [], [], []
        for t in q.X("raw").tables():
            rows.append(t["soma_dim_0"].to_numpy())
            cols.append(t["soma_dim_1"].to_numpy())
            vals.append(t["soma_data"].to_numpy())

rows, cols, vals = np.concatenate(rows), np.concatenate(cols), np.concatenate(vals).astype(np.float32)
pos = pd.Series(np.arange(len(ids)), index=ids)
M = sp.csr_matrix((vals, (pos.loc[rows].values, cols)), shape=(len(ids), int(var.soma_joinid.max()) + 1))
M = sp.diags(1e4 / np.maximum(np.asarray(M.sum(axis=1)).ravel(), 1)) @ M
picked = picked.drop_duplicates("soma_joinid").set_index("soma_joinid").loc[ids]
pb = {}
for cls, g in picked.groupby("ref_class", sort=False):
    per_ds = [np.asarray(M[pos.loc[gd.index].values].mean(axis=0)).ravel() for _, gd in g.groupby("dataset_id", observed=True)]
    pb[cls] = np.mean(per_ds, axis=0)
pb = pd.DataFrame(pb)
pb.index = var.set_index("soma_joinid").reindex(pb.index)["feature_id"].values
pb = pb[pb.index.notna()][list(CLASSES)]
pb.to_csv(os.path.join(OUT, "dura_reference_pseudobulk_cp10k.csv.gz"))
summ = picked.groupby("ref_class").agg(cells=("dataset_id", "size"), datasets=("dataset_id", "nunique"), source=("source", "first"))
summ.to_csv(os.path.join(OUT, "dura_reference_class_summary.csv"))
print(summ)
