"""07: are CD11b and ACSA2 samples with the same treatment/replicate number from the same animal?
Residual (within treatment x fraction) log-expression profiles are correlated across fractions over genes
expressed in both; matched pairs are compared with all 3x3 within-treatment pairings."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from common import *
ann, cnt = load_matrix(); ss = sample_sheet()
keep = ~ann.index.str.startswith("ERCC")
cpm = cnt[keep] / cnt[keep].sum() * 1e6
lc = np.log2(cpm + 1)
cd = ss.index[ss.fraction == "CD11b"]; ac = ss.index[ss.fraction == "ACSA2"]
expr = (cpm[cd].min(1) >= 20) & (cpm[ac].min(1) >= 20)
res = lc[expr].copy()
for g, s in ss.groupby("group").groups.items():
    res[s] = lc.loc[expr, s].sub(lc.loc[expr, s].mean(1), axis=0)
# exclude genes that differ between fractions strongly (carry-over-dominated) to reduce cross-contamination coupling
ratio = (lc.loc[expr, cd].mean(1) - lc.loc[expr, ac].mean(1)).abs()
genes = ratio.index[ratio < 1]
R = pd.DataFrame(index=cd, columns=ac, dtype=float)
for a in cd:
    for b in ac:
        R.loc[a, b] = np.corrcoef(res.loc[genes, a], res.loc[genes, b])[0, 1]
lab = lambda s: f"{ss.loc[s,'treatment']}_r{ss.loc[s,'replicate']}"
R.index = [f"{i} ({lab(i)})" for i in R.index]; R.columns = [f"{i} ({lab(i)})" for i in R.columns]
print(f"genes used: {len(genes)}")
pd.set_option("display.width", 250); print(R.round(3).to_string())
for t in ["NTC", "Spp1"]:
    rows = [i for i in R.index if f"({t}_" in i]; cols = [c for c in R.columns if f"({t}_" in c]
    sub = R.loc[rows, cols].values
    diag = np.mean(np.diag(sub))
    import itertools
    null = [np.mean([sub[i, p[i]] for i in range(3)]) for p in itertools.permutations(range(3))]
    best = max(null)
    print(f"{t}: matched-pair mean r={diag:.3f}; rank among 6 pairings={sorted(null, reverse=True).index(diag)+1}; "
          f"argmax each row matches replicate: {[int(np.argmax(sub[i])==i) for i in range(3)]}")
R.to_csv(os.path.join(ROOT, "08_crosstalk", "pairing_residual_correlation.tsv"), sep="\t", float_format="%.3f")
