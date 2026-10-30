"""13d - Is the SPP1 effect on the IFN and eosinophil scores explained by cell composition?

Per-animal scores from step 13c (mean z-scored VST of M_TYPE_I_IFN_EXTENDED / M_EOSINOPHIL_GRANULE members,
12 5xFAD animals). Each score is standardised and regressed on sex + treatment (+ composition signatures, each
standardised). SPP1_effect_SD = treatment coefficient in SD units. perm_p = exact within-sex permutation of the
treatment labels (C(6,3)^2 = 400 relabelings), two-sided on |t|.
Output: 13_pathway_links/IFN_eosinophil_robustness.csv
"""
import itertools
import os

import numpy as np
import pandas as pd
import statsmodels.api as sm

HERE = os.path.dirname(os.path.abspath(__file__))
CA = os.path.dirname(HERE)
TAB = os.path.join(CA, "tables")
OUT = os.path.join(CA, "13_pathway_links")
sc = pd.read_csv(os.path.join(OUT, "eosinophil_IFN_scores_per_animal.csv"), index_col=0)
meta = pd.read_csv(os.path.join(TAB, "sample_metadata_with_composition.tsv"), sep="\t", index_col=0)
S = pd.read_csv(os.path.join(TAB, "composition_signature_scores_all.tsv"), sep="\t", index_col=0)
F5 = sc.index
std = lambda x: (x - x.mean()) / x.std()
spp1 = (meta.loc[F5, "treatment"] == "SPP1").astype(float).values
male = (sc["sex"] == "M").astype(float).values

# all within-sex relabelings (3 SPP1 of 6 per sex)
idx_f, idx_m = np.where(male == 0)[0], np.where(male == 1)[0]
perms = []
for cf in itertools.combinations(idx_f, 3):
    for cm in itertools.combinations(idx_m, 3):
        v = np.zeros(len(F5)); v[list(cf) + list(cm)] = 1; perms.append(v)

def fit(y, trt, covs):
    X = sm.add_constant(np.column_stack([male, trt] + [std(S.loc[F5, c]).values for c in covs]))
    r = sm.OLS(y, X).fit()
    return r.params[2], r.tvalues[2], r.pvalues[2]

COVS = [[], ["fibroblast_stromal"], ["monocyte"], ["neutrophil"], ["monocyte", "neutrophil"], ["fibroblast_stromal", "monocyte"]]
rows = []
for outcome in ["IFN_score", "eosinophil_score"]:
    y = std(sc[outcome]).values
    for covs in COVS:
        b, t, p = fit(y, spp1, covs)
        tnull = np.array([abs(fit(y, v, covs)[1]) for v in perms])
        rows.append({"outcome": outcome, "covariates": "sex only" if not covs else "+ " + " + ".join(covs),
                     "SPP1_effect_SD": b, "p_ols": p, "perm_p": np.mean(tnull >= abs(t) - 1e-12)})
res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "IFN_eosinophil_robustness.csv"), index=False)
print(res.round(3).to_string(index=False))
