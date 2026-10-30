"""11: RNA-quality ("rerun") checks without BAMs.
(a) Does the SPP1-vs-NTC log2FC track gene length (degradation / fragmentation signature)?
(b) Per-sample long-vs-short gene index: log2 CPM ratio of the longest vs shortest quintile of expressed protein-coding genes.
(c) Does this index separate treatment within fraction, and does Spp1's own effect survive adjustment for it?"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from common import *
from scipy import stats
ann, cnt = load_matrix(); ss = sample_sheet()
L = pd.read_csv(os.path.join(DATA, "reference", "census_mouse_gene_lengths.csv.gz")).set_index("feature_id")["feature_length"]
keep = ~ann.index.str.startswith("ERCC")
cpm = cnt[keep] / cnt[keep].sum() * 1e6
pc = ann.index[(ann.gene_biotype == "protein_coding") & ~ann.gene_name.fillna("").str.startswith("mt-")]
out = []
for frac, d in [("CD11b", "05_DE_CD11b"), ("ACSA2", "06_DE_ACSA2")]:
    s = ss.index[ss.fraction == frac]
    de = pd.read_csv(os.path.join(ROOT, d, f"DE_Spp1_vs_NTC_{frac}_annotated.tsv"), sep="\t").set_index("gene_id")
    de = de[de.index.isin(pc) & de.index.isin(L.index)]
    de["len"] = np.log10(L.reindex(de.index).values)
    rho, p = stats.spearmanr(de.len, de.log2FC)
    # long vs short index
    e = cpm.loc[de.index, s]; ex = e[e.min(1) >= 10]
    q = pd.qcut(np.log10(L.reindex(ex.index)), 5, labels=False)
    idx = np.log2(ex[q == 4].sum() / ex[q == 0].sum())
    print(f"{frac}: Spearman(log10 length, log2FC) rho={rho:+.3f} (p={p:.1e}, n={len(de)})")
    for x in s:
        out.append(dict(sample=x, fraction=frac, treatment=ss.loc[x, "treatment"], long_short_log2=idx[x]))
    r = pd.DataFrame(out); r = r[r.fraction == frac]
    print(r.round(3).to_string(index=False))
    # length-binned median log2FC
    de["bin"] = pd.qcut(de.len, 5, labels=["Q1 short", "Q2", "Q3", "Q4", "Q5 long"])
    print(de.groupby("bin", observed=True).log2FC.median().round(3).to_string())
pd.DataFrame(out).to_csv(os.path.join(ROOT, "03_artifacts", "gene_length_quality_index.tsv"), sep="\t", index=False, float_format="%.4f")
