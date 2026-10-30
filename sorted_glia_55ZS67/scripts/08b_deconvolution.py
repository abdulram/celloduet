"""08b: NNLS deconvolution (RNA shares) against the Census brain pseudobulk reference.
Markers: per class, top-N genes by specificity (class / max other class >= 4), excluding mt-, ribosomal, Ig/TCR-V and
Hb genes. Each gene row is scaled by its reference maximum; two rounds of trimming the 10% worst-fit genes.
Stability is reported over N = 20, 30, 50 markers. Profiles sum to 1e4 per class, so coefficients ~ RNA shares."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from common import *
from scipy.optimize import nnls
ann, cnt = load_matrix(); ss = sample_sheet()
ref = pd.read_csv(os.path.join(DATA, "reference", "brain_pseudobulk_cp10k.csv.gz"), index_col=0)
ref = ref.loc[:, ref.sum() > 0]
ref.index = ref.index.str.split(".").str[0]
ref = ref.groupby(level=0).sum()
keep = ~ann.index.str.startswith("ERCC")
cpm = (cnt[keep] / cnt[keep].sum() * 1e4)  # per 1e4 to match reference scale
common = ref.index.intersection(cpm.index)
nm = ann.gene_name.reindex(common).fillna("")
bad = nm.str.match(r"^(mt-|Rp[sl]\d|Rplp|Igh|Igk|Igl|Trav|Trbv|Trgv|Trdv|Hb[ab]-)")
R = ref.loc[common[~bad]]; Y = cpm.loc[R.index]
R = R / R.sum() * 1e4  # renormalise on the shared gene universe
def markers(n):
    out = []
    for c in R.columns:
        other = R.drop(columns=c).max(1)
        spec = (R[c] + 0.01) / (other + 0.01)
        cand = spec[(spec >= 4) & (R[c] >= 0.5)].sort_values(ascending=False)
        out += list(cand.index[:n])
    return sorted(set(out))
def fit(y, M, genes):
    g = list(genes)
    for _ in range(3):
        w = 1 / M.loc[g].max(1)
        A = M.loc[g].mul(w, axis=0).values; b = y[g].values * w.values
        coef, _ = nnls(A, b)
        resid = np.abs(A @ coef - b)
        cut = np.quantile(resid, 0.9)
        g = [x for x, r in zip(g, resid) if r <= cut]
    return coef
res = []
for n in [20, 30, 50]:
    mk = markers(n)
    for s in Y.columns:
        c = fit(Y[s], R, mk)
        res.append(pd.Series(c / c.sum() * 100, index=R.columns, name=(n, s)))
D = pd.DataFrame(res)
D.index = pd.MultiIndex.from_tuples(D.index, names=["n_markers", "sample"])
D.to_csv(os.path.join(ROOT, "02_purity_celltyping", "deconvolution_rna_share_pct_all_settings.tsv"), sep="\t", float_format="%.2f")
med = D.groupby(level="sample").median().loc[ss.index]
rng = (D.groupby(level="sample").max() - D.groupby(level="sample").min()).loc[ss.index]
med.to_csv(os.path.join(ROOT, "02_purity_celltyping", "deconvolution_rna_share_pct_median.tsv"), sep="\t", float_format="%.2f")
pd.set_option("display.width", 250)
print("median RNA share % across N=20/30/50 markers:"); print(med.round(1).assign(trt=ss.treatment, frac=ss.fraction).to_string())
print("\nmax-min range across settings:"); print(rng.round(1).to_string())
# group comparison per fraction
for f in ["CD11b", "ACSA2"]:
    s = ss.index[ss.fraction == f]
    m = med.loc[s].groupby(ss.loc[s, "treatment"]).mean().T
    m["diff"] = m.Spp1 - m.NTC
    print(f"\n{f} mean RNA share by treatment:"); print(m[m.max(1) > 0.3].round(2).to_string())
