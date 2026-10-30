"""08a: adult mouse brain pseudobulk reference from CELLxGENE Census 2025-01-30 (public S3; adapted from the dura project).
Up to 1200 cells/class, <=300 per dataset; CP10K per cell, averaged per dataset then across datasets."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd, scipy.sparse as sp, tiledbsoma as soma
from census_ctx import open_census
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "reference"); os.makedirs(OUT, exist_ok=True)
N_PER_CLASS, N_PER_DATASET, SEED = 1200, 300, 7
BRAIN = ["brain", "central nervous system", "nervous system"]
CLASSES = {
 "microglia": ["microglial cell"],
 "BAM": ["perivascular macrophage", "border associated macrophage", "macrophage"],
 "monocyte": ["monocyte", "classical monocyte"],
 "neutrophil": ["neutrophil"],
 "T_NK": ["T cell", "natural killer cell", "CD8-positive, alpha-beta T cell", "CD4-positive, alpha-beta T cell"],
 "B_cell": ["B cell"],
 "astrocyte": ["astrocyte", "astrocyte of the cerebral cortex", "Bergmann glial cell"],
 "oligodendrocyte": ["oligodendrocyte"],
 "OPC": ["oligodendrocyte precursor cell"],
 "neuron": ["glutamatergic neuron", "GABAergic neuron", "neuron"],
 "endothelial": ["endothelial cell", "brain endothelial cell", "capillary endothelial cell"],
 "pericyte": ["pericyte", "brain pericyte"],
 "smooth_muscle": ["smooth muscle cell", "vascular associated smooth muscle cell"],
 "VLMC_fibroblast": ["vascular leptomeningeal cell", "fibroblast", "leptomeningeal cell"],
 "ependymal": ["ependymal cell"],
 "choroid_plexus": ["choroid plexus epithelial cell"],
}
BLOODOK = {"neutrophil", "monocyte", "B_cell", "T_NK"}
with open_census() as census:
    mm = census["census_data"]["mus_musculus"]
    obs = mm.obs.read(value_filter="is_primary_data == True",
                      column_names=["soma_joinid", "cell_type", "tissue_general", "tissue", "dataset_id", "development_stage"]).concat().to_pandas()
    obs = obs[obs.tissue_general.astype(str) != "embryo"]
    var = mm.ms["RNA"].var.read(column_names=["soma_joinid", "feature_id", "feature_name"]).concat().to_pandas()
    picked = []
    for cls, cts in CLASSES.items():
        sub = obs[obs.cell_type.isin(cts)]
        sb = sub[sub.tissue_general.astype(str).isin(BRAIN)]
        if len(sb) >= 100 or cls not in BLOODOK: sub = sb
        cap = N_PER_CLASS if sub.dataset_id.nunique() == 1 else N_PER_DATASET
        parts = [g.sample(min(len(g), cap), random_state=SEED) for _, g in sub.groupby("dataset_id", observed=True)]
        sub = pd.concat(parts) if parts else sub
        if len(sub) > N_PER_CLASS: sub = sub.sample(N_PER_CLASS, random_state=SEED)
        picked.append(sub.assign(ref_class=cls))
        print(f"{cls:18s} cells={len(sub):5d} datasets={sub.dataset_id.nunique()} tissues={sub.tissue_general.astype(str).unique()[:3]}", flush=True)
    picked = pd.concat(picked).drop_duplicates("soma_joinid")
    ids = np.sort(picked.soma_joinid.unique())
    with mm.axis_query("RNA", obs_query=soma.AxisQuery(coords=(ids,))) as q:
        rows, cols, vals = [], [], []
        for t in q.X("raw").tables():
            rows.append(t["soma_dim_0"].to_numpy()); cols.append(t["soma_dim_1"].to_numpy()); vals.append(t["soma_data"].to_numpy())
rows, cols, vals = np.concatenate(rows), np.concatenate(cols), np.concatenate(vals).astype(np.float32)
pos = pd.Series(np.arange(len(ids)), index=ids)
M = sp.csr_matrix((vals, (pos.loc[rows].values, cols)), shape=(len(ids), int(var.soma_joinid.max()) + 1))
M = sp.diags(1e4 / np.maximum(np.asarray(M.sum(axis=1)).ravel(), 1)) @ M
picked = picked.set_index("soma_joinid").loc[ids]
pb = {}
for cls, g in picked.groupby("ref_class", sort=False):
    pb[cls] = np.mean([np.asarray(M[pos.loc[gd.index].values].mean(axis=0)).ravel() for _, gd in g.groupby("dataset_id", observed=True)], axis=0)
pb = pd.DataFrame(pb); pb.index = var.set_index("soma_joinid").reindex(pb.index)["feature_id"].values
pb = pb[pb.index.notna()]
pb.to_csv(os.path.join(OUT, "brain_pseudobulk_cp10k.csv.gz"))
picked.groupby("ref_class").agg(cells=("dataset_id", "size"), datasets=("dataset_id", "nunique")).to_csv(os.path.join(OUT, "brain_class_summary.csv"))
print("done", pb.shape)
