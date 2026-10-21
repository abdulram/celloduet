"""01: sample sheet, technical QC metrics, sex check, ERCC spike-ins."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from common import *

ann, cnt = load_matrix()
ss = sample_sheet()
mp = pd.read_csv(os.path.join(DATA, "55ZS67-mapping-stats-reads.csv"), index_col=0)
mp.index = [f"55ZS67_{s.split('_')[-1]}" for s in mp.index]

name = ann.gene_name
is_mt = name.str.startswith("mt-")
is_ribo = name.str.match(r"^Rp[sl]\d+[a-z]?\d*$|^Rplp[012]$|^Rpsa$")
is_ercc = ann.index.str.startswith("ERCC")
lib = cnt.sum()
q = pd.DataFrame(index=cnt.columns)
q["total_reads"] = mp.sum(1)
q["uniquely_mapped"] = mp["Uniquely Mapped"]
q["pct_unique"] = 100 * mp["Uniquely Mapped"] / q.total_reads
q["pct_multi"] = 100 * mp["Multi-mapped"] / q.total_reads
q["pct_unmapped"] = 100 * mp["Unmapped"] / q.total_reads
q["assigned_counts"] = lib
q["pct_assigned_of_unique"] = 100 * lib / mp["Uniquely Mapped"]
q["genes_ge5"] = (cnt[~is_ercc] >= 5).sum()
q["pct_mito"] = 100 * cnt[is_mt].sum() / lib
q["pct_ribo"] = 100 * cnt[is_ribo].sum() / lib
q["pct_ERCC"] = 100 * cnt[is_ercc].sum() / lib
for bt in ["protein_coding", "lncRNA", "Mt_rRNA", "rRNA", "misc_RNA", "snoRNA"]:
    q[f"pct_{bt}"] = 100 * cnt[ann.gene_biotype == bt].sum() / lib
# top-gene dominance (library complexity)
top = cnt[~is_ercc].apply(lambda c: c.sort_values(ascending=False).iloc[:20].sum() / c.sum() * 100)
q["pct_top20_genes"] = top
cpm = cnt / lib * 1e6
sexg = ["Xist", "Tsix", "Ddx3y", "Uty", "Kdm5d", "Eif2s3y"]
for g in sexg:
    ids = name[name == g].index
    q[f"cpm_{g}"] = cpm.loc[ids].sum()
q["sex_call"] = np.where((q.cpm_Xist > 10) & (q[["cpm_Ddx3y", "cpm_Uty", "cpm_Kdm5d", "cpm_Eif2s3y"]].sum(1) < 2), "F",
                np.where((q.cpm_Xist < 2) & (q[["cpm_Ddx3y", "cpm_Uty", "cpm_Kdm5d", "cpm_Eif2s3y"]].sum(1) > 5), "M", "ambiguous"))
q = ss.join(q)
q.to_csv(os.path.join(ROOT, "01_QC", "qc_metrics.tsv"), sep="\t", float_format="%.4g")
q[["facility_name", "treatment", "fraction", "replicate", "sex_call"]].to_csv(os.path.join(ROOT, "tables", "sample_sheet.tsv"), sep="\t")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print(q.drop(columns=["facility_name", "group"]).round(2).to_string())

# ERCC: dose response sanity (detected spikes only)
er = cnt[is_ercc]
print("\nERCC detected (>=10 counts in any sample):", int((er.max(1) >= 10).sum()), "of", len(er))
print("ERCC total counts per sample:\n", er.sum().round(0).to_string())

# Fractional counts = multimapper splitting?
frac = (cnt.values % 1 != 0).mean()
print(f"\nfraction of non-integer count entries: {frac:.4f}")
