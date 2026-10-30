"""02/03: lineage panel scores, cross-fraction carry-over, dissociation score, and
technical/purity metrics vs treatment (per fraction). n=3 vs 3 -> exact permutation over 20 labelings."""
import os, sys, itertools
sys.path.insert(0, os.path.dirname(__file__))
from common import *
from scipy import stats
from gene_sets_loader import load_sets

ann, cnt = load_matrix()
ss = sample_sheet()
q = pd.read_csv(os.path.join(ROOT, "01_QC", "qc_metrics.tsv"), sep="\t", index_col=0)
keep = ~ann.index.str.startswith("ERCC")
cpm = cnt[keep] / cnt[keep].sum() * 1e6
cpm.index = ann.gene_name[keep]
cpm = cpm.groupby(level=0).sum()
lc = np.log2(cpm + 1)
sets = load_sets()

def exact_perm(x, lab):
    """two-sided exact permutation p for difference in means, all C(6,3)=20 relabelings."""
    x = np.asarray(x, float); lab = np.asarray(lab)
    obs = x[lab == 1].mean() - x[lab == 0].mean()
    n1 = lab.sum(); idx = range(len(x)); null = []
    for c in itertools.combinations(idx, n1):
        m = np.zeros(len(x), bool); m[list(c)] = True
        null.append(x[m].mean() - x[~m].mean())
    null = np.array(null)
    return obs, (np.abs(null) >= abs(obs) - 1e-12).mean()

def auc(x, lab):
    x = np.asarray(x); lab = np.asarray(lab)
    a = x[lab == 1]; b = x[lab == 0]
    return np.mean([(i > j) + 0.5 * (i == j) for i in a for j in b])

rows, scores, valid = [], {}, []
for frac in ["CD11b", "ACSA2"]:
    s = ss.index[ss.fraction == frac]
    lab = (ss.loc[s, "treatment"] == "Spp1").astype(int).values
    z = lc[s].sub(lc[s].mean(1), axis=0).div(lc[s].std(1).replace(0, np.nan), axis=0)
    for k, genes in sets.items():
        g = [x for x in genes if x in lc.index and cpm.loc[x, s].mean() >= 1]
        if len(g) < 2: continue
        sc = z.loc[g].mean()
        scores[(frac, k)] = sc
        ir = []
        for x in g:
            rest = z.loc[[y for y in g if y != x]].mean()
            ir.append(np.corrcoef(z.loc[x], rest)[0, 1])
        valid.append(dict(fraction=frac, set=k, n_detected=len(g), median_item_rest_r=np.nanmedian(ir),
                          mean_log2cpm=lc.loc[g, s].mean().mean(),
                          passes=(len(g) >= 3) and (np.nanmedian(ir) >= 0.4)))
sc = pd.DataFrame(scores).T
val = pd.DataFrame(valid)
val.to_csv(os.path.join(ROOT, "02_purity_celltyping", "panel_validation.tsv"), sep="\t", index=False, float_format="%.3f")

# Cross-fraction carry-over: sample CPM / reference-fraction mean CPM, median over specific markers
mg_spec = ["P2ry12", "Tmem119", "Siglech", "Gpr34", "Cx3cr1", "Hexb", "Olfml3", "Itgam", "Ptprc", "Trem2", "Tyrobp", "C1qa", "Csf1r"]
as_spec = ["Aldh1l1", "Slc1a2", "Slc1a3", "Gja1", "Aqp4", "Atp1b2", "Aldoc", "Sox9", "Gjb6", "Fgfr3"]
ol_spec = ["Plp1", "Mbp", "Mog", "Mobp", "Mag", "Cldn11"]
cd = ss.index[ss.fraction == "CD11b"]; ac = ss.index[ss.fraction == "ACSA2"]
pur = pd.DataFrame(index=ss.index)
pur["microglial_RNA_pct_vs_CD11b"] = (cpm.loc[mg_spec].div(cpm.loc[mg_spec, cd].mean(1), axis=0)).median() * 100
pur["astro_marker_pct_vs_ACSA2"] = (cpm.loc[as_spec].div(cpm.loc[as_spec, ac].mean(1), axis=0)).median() * 100
pur["oligo_marker_pct_vs_ACSA2"] = (cpm.loc[ol_spec].div(cpm.loc[ol_spec, ac].mean(1), axis=0)).median() * 100
pur["Plp1_cpm"] = cpm.loc["Plp1"]; pur["Atp1b2_cpm"] = cpm.loc["Atp1b2"]
pur["Plp1_over_Atp1b2"] = cpm.loc["Plp1"] / cpm.loc["Atp1b2"]
pur["Ttr_cpm"] = cpm.loc["Ttr"]
pur["neutrophil_cpm_sum"] = cpm.loc[["S100a8", "S100a9", "Ngp", "Camp", "Retnlg", "Ltf", "Lcn2"]].sum()
pur["hemoglobin_cpm"] = cpm.loc[["Hbb-bs", "Hba-a1", "Hbb-bt", "Hba-a2"]].sum()
pur["dissociation_IEG_log2cpm_mean"] = lc.loc[[g for g in sets["X_Dissociation_IEG"] if g in lc.index]].mean()
pur = ss[["treatment", "fraction"]].join(pur)
pur.to_csv(os.path.join(ROOT, "02_purity_celltyping", "purity_carryover_per_sample.tsv"), sep="\t", float_format="%.3f")

# Tests: every technical / purity metric and every panel score vs treatment within fraction
metrics = ["pct_unique", "pct_multi", "pct_unmapped", "pct_assigned_of_unique", "assigned_counts", "genes_ge5",
           "pct_mito", "pct_ribo", "pct_Mt_rRNA", "pct_rRNA", "pct_lncRNA", "pct_top20_genes"]
for frac, s in [("CD11b", cd), ("ACSA2", ac)]:
    lab = (ss.loc[s, "treatment"] == "Spp1").astype(int).values
    for m in metrics:
        d, p = exact_perm(q.loc[s, m], lab)
        rows.append(dict(fraction=frac, kind="technical", metric=m, NTC_mean=q.loc[s, m][lab == 0].mean(),
                         Spp1_mean=q.loc[s, m][lab == 1].mean(), diff=d, AUC_Spp1_gt_NTC=auc(q.loc[s, m], lab), perm_p=p))
    for m in pur.columns[2:]:
        d, p = exact_perm(pur.loc[s, m], lab)
        rows.append(dict(fraction=frac, kind="purity/artifact", metric=m, NTC_mean=pur.loc[s, m][lab == 0].mean(),
                         Spp1_mean=pur.loc[s, m][lab == 1].mean(), diff=d, AUC_Spp1_gt_NTC=auc(pur.loc[s, m], lab), perm_p=p))
    for (f, k), v in sc.iterrows():
        if f != frac: continue
        v = v[s]
        d, p = exact_perm(v, lab)
        rows.append(dict(fraction=frac, kind="panel_score_z", metric=k, NTC_mean=v[lab == 0].mean(),
                         Spp1_mean=v[lab == 1].mean(), diff=d, AUC_Spp1_gt_NTC=auc(v, lab), perm_p=p))
res = pd.DataFrame(rows)
res.to_csv(os.path.join(ROOT, "03_artifacts", "metric_vs_treatment_tests.tsv"), sep="\t", index=False, float_format="%.4g")
sc.T.to_csv(os.path.join(ROOT, "02_purity_celltyping", "panel_scores_z_within_fraction.tsv"), sep="\t", float_format="%.3f")

pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
print(pur.round(2).to_string())
print("\nNote: with 3 vs 3 the smallest attainable two-sided exact p is 0.1 (AUC 0 or 1).")
print(res[(res.AUC_Spp1_gt_NTC.isin([0, 1]))].round(3).to_string())
print("\nPanel validation (fails):")
print(val[~val.passes & val.set.str.startswith("L_")].round(2).to_string())
