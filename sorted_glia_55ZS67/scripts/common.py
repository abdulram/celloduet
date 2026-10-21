"""Shared loaders for the 55ZS67 sorted CD11b+/ACSA-2+ bulk RNA-seq analysis."""
import os, re
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
COLORS = {"NTC": "#eb6834", "Spp1": "#1baf7a"}  # 5xFAD_NTC / 5xFAD_SPP1 palette from dura project
np.random.seed(1)


def load_matrix():
    m = pd.read_csv(os.path.join(DATA, "55ZS67-expression-matrix.tsv.gz"), sep="\t")
    ann = m[["gene_id", "gene_name", "gene_biotype"]].set_index("gene_id")
    cnt = m.set_index("gene_id")[[c for c in m.columns if c.endswith("_count")]]
    cnt.columns = [c.replace("_count", "") for c in cnt.columns]
    order = sorted(cnt.columns, key=lambda s: int(s.split("_")[1]))
    return ann, cnt[order]


def sample_sheet():
    """Sample metadata parsed from the facility's sample names (rerun_<trt>_<fraction>_<rep>_<n>)."""
    mp = pd.read_csv(os.path.join(DATA, "55ZS67-mapping-stats-reads.csv"), index_col=0)
    rows = []
    for name in mp.index:
        _, trt, frac, rep, n = name.split("_")
        rows.append(dict(sample=f"55ZS67_{n}", facility_name=name, treatment=trt,
                         fraction={"microglia": "CD11b", "astrocyte": "ACSA2"}[frac],
                         replicate=int(rep), group=f"{trt}_{frac}"))
    return pd.DataFrame(rows).set_index("sample")
