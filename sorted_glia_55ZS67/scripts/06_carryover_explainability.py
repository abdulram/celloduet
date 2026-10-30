"""06: can cross-fraction carry-over explain each DE gene?
Model: sample CPM = intrinsic + a_sample * (other-fraction mean CPM), with a_sample the measured carry-over share.
Intrinsic is estimated from the NTC arm; predicted log2FC is what carry-over alone would produce.
Also a contamination-free PCA (genes >= 4x higher in the other fraction removed)."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from common import *

ann, cnt = load_matrix()
ss = sample_sheet()
pur = pd.read_csv(os.path.join(ROOT, "02_purity_celltyping", "purity_carryover_per_sample.tsv"), sep="\t", index_col=0)
keep = ~ann.index.str.startswith("ERCC")
cpm = cnt[keep] / cnt[keep].sum() * 1e6
cd = ss.index[ss.fraction == "CD11b"]; ac = ss.index[ss.fraction == "ACSA2"]
out = {}
for frac, s, other, acol, d in [("CD11b", cd, ac, "astro_marker_pct_vs_ACSA2", "05_DE_CD11b"),
                                ("ACSA2", ac, cd, "microglial_RNA_pct_vs_CD11b", "06_DE_ACSA2")]:
    de = pd.read_csv(os.path.join(ROOT, d, f"DE_Spp1_vs_NTC_{frac}.tsv"), sep="\t").set_index("gene_id")
    A = cpm.loc[de.index, other].mean(1)
    a = pur.loc[s, acol] / 100
    ntc = [x for x in s if ss.loc[x, "treatment"] == "NTC"]; sp = [x for x in s if ss.loc[x, "treatment"] == "Spp1"]
    obs_ntc = cpm.loc[de.index, ntc].mean(1)
    I = (obs_ntc - a[ntc].mean() * A).clip(lower=0.05 * obs_ntc)
    pred = np.log2((I + a[sp].mean() * A + 0.5) / (I + a[ntc].mean() * A + 0.5))
    de["other_fraction_cpm"] = A.values
    de["ratio_other_over_this"] = (A / cpm.loc[de.index, s].mean(1).replace(0, np.nan)).values
    de["carryover_pred_log2FC"] = pred.values
    de["carryover_share_of_obs"] = (pred / de.log2FC).values
    de["carryover_explains"] = (np.sign(pred) == np.sign(de.log2FC)) & (de.carryover_share_of_obs.abs() >= 0.5)
    de.reset_index().to_csv(os.path.join(ROOT, d, f"DE_Spp1_vs_NTC_{frac}_annotated.tsv"), sep="\t", index=False, float_format="%.4g")
    sig = de[de.padj < 0.1]
    print(f"== {frac}: carry-over a NTC={a[ntc].mean()*100:.2f}%  Spp1={a[sp].mean()*100:.2f}%  | padj<0.1: {len(sig)}; "
          f"explained by carry-over: {int(sig.carryover_explains.sum())}; "
          f"genes >=4x higher in other fraction: {int((sig.ratio_other_over_this >= 4).sum())}")
    nonc = sig[~sig.carryover_explains & (sig.ratio_other_over_this < 4)]
    print("   not explained by carry-over and not other-fraction-enriched:")
    print(nonc[["gene", "baseMean", "log2FC", "padj", "carryover_pred_log2FC", "ratio_other_over_this", "loo_max_pvalue"]].round(3).to_string())
    out[frac] = de

# Rank of Spp1 effect among highly expressed CD11b genes
de = out["CD11b"]
hi = de[(de.baseMean > 1000)]
print("\nCD11b genes baseMean>1000:", len(hi), "| Spp1 rank by |log2FC|:", int((hi.log2FC.abs() >= abs(hi.loc[hi.gene == 'Spp1', 'log2FC'].iloc[0])).sum()))
print(hi.reindex(hi.log2FC.abs().sort_values(ascending=False).index)[["gene", "baseMean", "log2FC", "padj", "carryover_pred_log2FC"]].head(12).round(3).to_string())
