"""09 - Pathology associations.

Searches every metadata source in the project for mouse-level pathology measurements (plaque burden, Abeta40/42,
Iba1, GFAP, cognition, vascular / lymphatic measurements, age). None exist: the only sample annotations are the
facility sample names (animal, genotype, sex, treatment) and library QC. Pathology analyses therefore cannot be
performed; this script documents the search and writes a placeholder Figure 11.
"""
import glob
import os
import re

import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CA = os.path.dirname(HERE)
ROOT = os.path.dirname(CA)
OUT = os.path.join(CA, "09_pathology_associations")
TERMS = r"plaque|abeta|aβ|ab40|ab42|amyloid|iba1|gfap|cogniti|morris|y-maze|pathology|lymphatic_measure|age|weeks|months"
hits = []
for f in glob.glob(os.path.join(ROOT, "data", "**", "*"), recursive=True):
    if os.path.isfile(f) and f.endswith((".csv", ".tsv", ".txt", ".xlsx")) and "reference" not in f and "genesets" not in f:
        try:
            head = open(f, errors="ignore").readline().lower()
        except Exception:
            continue
        cols = [c for c in re.split(r"[\t,]", head) if re.search(TERMS, c)]
        hits.append({"file": os.path.relpath(f, ROOT), "pathology_like_columns": ";".join(cols) or "none"})
pd.DataFrame(hits).to_csv(os.path.join(OUT, "pathology_metadata_search.tsv"), sep="\t", index=False)
print(pd.DataFrame(hits).to_string())
fig, ax = plt.subplots(figsize=(8, 2.2))
ax.axis("off")
ax.text(0.5, 0.5, "Figure 11 - Pathology associations: NOT POSSIBLE\nNo plaque, Aβ, gliosis, cognitive, vascular/lymphatic or age data\nwere found in the project (see pathology_metadata_search.tsv).",
        ha="center", va="center", fontsize=11)
fig.savefig(os.path.join(OUT, "Fig11_pathology_NOT_AVAILABLE.png"), dpi=150, bbox_inches="tight")
