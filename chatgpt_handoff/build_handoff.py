"""Assemble the self-contained handoff package (data/ + figures/) from the project outputs.

Run from the repository root:  python3 chatgpt_handoff/build_handoff.py
All outputs are CSV (gzipped when large). See README.md for the data dictionary.
"""
import glob
import os
import shutil

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H = os.path.join(ROOT, "chatgpt_handoff")
D = os.path.join(H, "data")
FIG = os.path.join(H, "figures")
CA = os.path.join(ROOT, "composition_confounding_analysis")
T = os.path.join(CA, "tables")
R = os.path.join(ROOT, "results")
for d in [D, FIG]:
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)


def rd(path, **kw):
    sep = "\t" if path.endswith((".tsv", ".tsv.gz")) else ","
    return pd.read_csv(path, sep=sep, **kw)


def out(df, section, name, index=False):
    os.makedirs(os.path.join(D, section), exist_ok=True)
    p = os.path.join(D, section, name)
    df.to_csv(p, index=index)
    return p


def cp(src, section, name):
    df = rd(src)
    if df.columns[0].startswith("Unnamed"):
        df = df.rename(columns={df.columns[0]: "id"})
    return out(df, section, name)


# ------------------------------------------------------------ 00 core
raw = pd.read_csv(os.path.join(ROOT, "data", "VT3YXB-expression-matrix.tsv"), sep="\t")
raw = raw[~raw.gene_id.str.startswith("ERCC")]
raw["gene_name"] = raw.gene_name.fillna(raw.gene_id)
libs = [f"VT3YXB_{i}" for i in range(1, 19)]
cm = rd(os.path.join(T, "sample_metadata_with_composition.tsv"), index_col=0)
lab = dict(zip(cm.index, cm["sample"]))
counts = raw[["gene_id", "gene_name", "gene_biotype"] + [f"{l}_count" for l in libs]].rename(columns={f"{l}_count": lab[l] for l in libs})
cpm = raw[["gene_id", "gene_name", "gene_biotype"] + [f"{l}_cpm" for l in libs]].rename(columns={f"{l}_cpm": lab[l] for l in libs})
out(counts, "00_core", "counts_raw_all_genes.csv.gz")
out(cpm, "00_core", "cpm_all_genes.csv.gz")
vst = rd(os.path.join(CA, "01_QC", "vst_blind.tsv.gz"), index_col=0)
g2n = raw.set_index("gene_id").gene_name
v = vst.rename(columns=lab)
v.insert(0, "gene_name", g2n.reindex(v.index).values)
out(v.reset_index().rename(columns={"index": "gene_id"}), "00_core", "vst_blind_filtered_genes.csv.gz")

# master per-sample table
m = cm.copy()
old = rd(os.path.join(R, "qc_and_composition.csv"), index_col=0)
m["Spp1_cpm"] = old["Spp1_cpm"]
m["genotype_signature_score_initial"] = old["genotype_signature_score"]
for c in [c for c in old if c.startswith("mod_")]:
    m[f"initial_marker_module__{c[4:]}"] = old[c]
sig = rd(os.path.join(T, "composition_signature_scores_all.tsv"), index_col=0)
m = m.join(sig.add_prefix("signature__"))
dg = rd(os.path.join(T, "deconvolution_grouped.tsv"), index_col=0)
m = m.join(dg.add_prefix("deconv_dura_ref__"))
fgc = rd(os.path.join(R, "fibroblast", "deconvolution_grouped_categories.csv"), index_col=0)
m = m.join(fgc.add_prefix("deconv_fibro_ensemble__"))
fsc = rd(os.path.join(R, "fibroblast", "fibroblast_scores_per_sample.csv"), index_col=0)
m["fibroblast_score_strict_genes"] = fsc["fibroblast_score"]
m["fibroblast_fraction_deconv"] = fsc["fibroblast_fraction_deconv"]
ecm = rd(os.path.join(R, "fibroblast", "ecm_score_per_sample.csv"), index_col=0)
m["fibroblast_ECM_score"] = ecm["fibroblast_ECM_leading_edge_score"]
tg = rd(os.path.join(T, "tgfb_panel_scores.tsv"), index_col=0)
m = m.join(tg)
m = m.drop(columns=["age", "batch", "dissection"])
m.insert(0, "library_id", m.index)
out(m, "00_core", "sample_master_table.csv")

# ------------------------------------------------------------ 01 QC / PCA
cp(os.path.join(T, "transcriptome_PC_variance.tsv"), "01_qc_pca", "transcriptome_PC_variance.csv")
cp(os.path.join(CA, "01_QC", "transcriptome_PC_top_genes.tsv"), "01_qc_pca", "transcriptome_PC_top_loading_genes.csv")
cp(os.path.join(CA, "01_QC", "transcriptome_PCA_loadings.tsv"), "01_qc_pca", "transcriptome_PCA_loadings_top2000.csv")
top = vst.var(axis=1).sort_values(ascending=False).index[:2000]
corr = vst.loc[top].rename(columns=lab).corr()
out(corr, "01_qc_pca", "sample_correlation_pearson_top2000_vst.csv", index=True)
ms = rd(os.path.join(ROOT, "data", "VT3YXB-mapping-stats.csv"))
out(ms.rename(columns={ms.columns[0]: "sample"}), "01_qc_pca", "mapping_stats_reads.csv")

# ------------------------------------------------------------ 02 initial whole-tissue analysis (results/)
for f in glob.glob(os.path.join(R, "DE_*.csv")) + glob.glob(os.path.join(R, "GSEA_*.csv")):
    cp(f, "02_initial_whole_tissue", os.path.basename(f))
for f in ["pc_covariate_spearman.csv", "genotype_DE_sensitivity.csv", "genotype_signature_gene_classes.csv", "spp1_stats.csv", "DE_summary.csv"]:
    cp(os.path.join(R, f), "02_initial_whole_tissue", f)

# ------------------------------------------------------------ 03 composition
for f, n in [("composition_signature_scores_all.tsv", "signature_scores_all20.csv"), ("signature_validation.tsv", "signature_validation.csv"),
             ("signature_marker_item_rest.tsv", "signature_marker_genes_item_rest.csv"), ("genotype_vs_composition.tsv", "genotype_vs_composition.csv"),
             ("composition_PC_loadings.tsv", "composition_PC_loadings.csv"), ("composition_PC_variance.tsv", "composition_PC_variance.csv"),
             ("txPC_vs_composition_spearman_rho.tsv", "txPC_vs_composition_spearman_rho.csv"), ("txPC_vs_composition_spearman_p.tsv", "txPC_vs_composition_spearman_p.csv"),
             ("txPC_variance_explained_by_composition.tsv", "txPC_variance_explained_by_composition.csv"),
             ("deconvolution_proportions.tsv", "deconv_dura_reference_proportions_20class.csv"), ("deconvolution_validation.tsv", "deconv_dura_reference_validation.csv"),
             ("deconvolution_grouped.tsv", "deconv_dura_reference_grouped.csv")]:
    cp(os.path.join(T, f), "03_composition", n)

# ------------------------------------------------------------ 04 genotype: baseline / adjusted / matched / position test
for f in glob.glob(os.path.join(T, "DE_*_5xFAD_vs_WT.tsv")):
    df = rd(f)
    out(df, "04_genotype_confounding", os.path.basename(f).replace(".tsv", ".csv.gz"))
for f, n in [("genotype_effects_baseline_vs_adjusted.tsv.gz", "genotype_effects_baseline_vs_adjusted.csv.gz"),
             ("adjusted_models_summary.tsv", "adjusted_models_summary.csv"), ("adjusted_models_collinearity.tsv", "adjusted_models_collinearity.csv"),
             ("common_support_cohorts_animals.tsv", "common_support_cohorts_animals.csv"), ("common_support_DE_summary.tsv", "common_support_DE_summary.csv"),
             ("technical_axis_animal_position_test.tsv", "technical_axis_animal_position_summary.csv")]:
    cp(os.path.join(T, f), "04_genotype_confounding", n)
# per-gene positions of each WT animal on the 6,787 baseline genes (Fig. animal-position histogram)
base = rd(os.path.join(T, "DE_baseline_5xFAD_vs_WT.tsv"))
sg = base.gene_id[base.padj < 0.05].values
WT = cm.index[cm.group == "WT_NTC"]
FX = cm.index[cm.group == "5xFAD_NTC"]
pos = {}
for a in WT:
    rest = WT.drop(a)
    pos[f"WT_{cm.animal[a]}"] = (vst.loc[sg, a] - vst.loc[sg, rest].mean(axis=1)) / (vst.loc[sg, FX].mean(axis=1) - vst.loc[sg, rest].mean(axis=1))
pos = pd.DataFrame(pos)
pos.insert(0, "gene_name", g2n.reindex(pos.index).values)
out(pos.reset_index().rename(columns={"index": "gene_id"}), "04_genotype_confounding", "technical_axis_animal_position_per_gene.csv.gz")

# ------------------------------------------------------------ 05 pathways + cell state
for f in glob.glob(os.path.join(T, "pathways_camera_*.tsv")):
    cp(f, "05_pathways_cell_state", os.path.basename(f).replace(".tsv", ".csv"))
for f, n in [("pathway_robustness_matrix.tsv", "pathway_robustness_matrix.csv"), ("pathway_module_scores.tsv", "pathway_module_scores_per_sample.csv"),
             ("cell_state_vs_abundance.tsv", "cell_state_vs_abundance.csv"), ("cell_state_scores.tsv", "cell_state_scores_per_sample.csv"),
             ("gene_sets.tsv", "gene_sets_long.csv")]:
    cp(os.path.join(T, f), "05_pathways_cell_state", n)

# ------------------------------------------------------------ 06 fibroblast-specific analysis
F = os.path.join(R, "fibroblast")
for f in ["fibroblast_specific_genes.csv", "fibroblast_specificity_all_genes.csv.gz", "deconvolution_proportions.csv", "deconvolution_grouped_categories.csv",
          "deconvolution_fit_quality.csv", "deconvolution_markers.csv", "fibroblast_scores_per_sample.csv", "ecm_score_per_sample.csv",
          "ecm_score_spp1_robustness.csv", "ecm_score_vs_spp1.csv", "fibroblast_genes_hallmark_overlap.csv", "DE_summary.csv"] + \
         [os.path.basename(x) for x in glob.glob(os.path.join(F, "DE_*fibnorm*.csv")) + glob.glob(os.path.join(F, "GSEA_*.csv")) + glob.glob(os.path.join(F, "compare_*.csv"))]:
    cp(os.path.join(F, f), "06_fibroblast_specific", f.replace("deconvolution_", "deconv_fibro_ensemble_"))
fv = rd(os.path.join(F, "fibroblast_normalised_vst.csv"), index_col=0).rename(columns=lab)
fv.insert(0, "gene_name", g2n.reindex(fv.index).values)
out(fv.reset_index().rename(columns={"index": "gene_id"}), "06_fibroblast_specific", "fibroblast_normalised_vst.csv")
for src, n in [(os.path.join(ROOT, "data", "reference", "census_pseudobulk_cp10k.csv.gz"), "reference_atlas_33class_pseudobulk_cp10k.csv.gz"),
               (os.path.join(CA, "reference", "dura_reference_pseudobulk_cp10k.csv.gz"), "reference_dura_20class_pseudobulk_cp10k.csv.gz")]:
    ref = pd.read_csv(src, index_col=0)
    ref.insert(0, "gene_name", g2n.reindex(ref.index).values)
    out(ref.reset_index().rename(columns={"index": "gene_id"}), "06_fibroblast_specific", n)
cp(os.path.join(ROOT, "data", "reference", "census_class_summary.csv"), "06_fibroblast_specific", "reference_atlas_33class_summary.csv")
cp(os.path.join(CA, "reference", "dura_reference_class_summary.csv"), "06_fibroblast_specific", "reference_dura_20class_summary.csv")

# ------------------------------------------------------------ 07 SPP1 on the 5xFAD signature (fibroblast genes)
for f in glob.glob(os.path.join(R, "fibroblast_5xfad_spp1", "*.csv")):
    cp(f, "07_spp1_vs_5xfad_signature", os.path.basename(f))

# ------------------------------------------------------------ 08 TGF-beta panel
for f in glob.glob(os.path.join(T, "tgfb_*.tsv")):
    cp(f, "08_tgfb_panel", os.path.basename(f).replace(".tsv", ".csv"))

# ------------------------------------------------------------ 09 SPP1 matrix-gene summary (consolidated)
G = ["Cldn11", "Fbln1", "Thbd", "Fn1", "Matn4", "Ifit3b", "Rgs4", "Prg4", "Col3a1", "Lum", "Mmp2", "Col16a1", "Loxl3", "Spp1"]


def pick(df, g, c):
    r = df[df.gene_name == g]
    return r[c].iloc[0] if len(r) else np.nan


w = rd(os.path.join(R, "DE_5xFAD_SPP1_vs_5xFAD_NTC_base.csv"))
wa = rd(os.path.join(R, "DE_5xFAD_SPP1_vs_5xFAD_NTC_adjusted.csv"))
f1 = rd(os.path.join(F, "DE_5xFAD_SPP1_vs_5xFAD_NTC_fibnorm.csv"))
f2 = rd(os.path.join(F, "DE_5xFAD_SPP1_vs_5xFAD_NTC_fibnorm_mitoadj.csv"))
tp = rd(os.path.join(T, "tgfb_panel_effects.tsv")).set_index("readout")
cpmg = cpm.set_index("gene_name")
lg = np.log2(cpmg[[lab[l] for l in libs]] + 1)
ntc = [lab[l] for l in cm.index[cm.group == "5xFAD_NTC"]]
spp = [lab[l] for l in cm.index[cm.group == "5xFAD_SPP1"]]
rows = []
for g in G:
    x = lg.loc[g] if g in lg.index else None
    if isinstance(x, pd.DataFrame):
        x = x.iloc[0]
    sx = {}
    for s in "FM":
        a = [c for c in ntc if f"_{s}_" in c]
        b = [c for c in spp if f"_{s}_" in c]
        sx[s] = np.median(x[b]) - np.median(x[a]) if x is not None else np.nan
    rows.append({"gene": g, "median_CPM_5xFAD_NTC": np.median(cpmg.loc[g, ntc]) if g in cpmg.index else np.nan,
                 "median_CPM_5xFAD_SPP1": np.median(cpmg.loc[g, spp]) if g in cpmg.index else np.nan,
                 "whole_tissue_log2FC": pick(w, g, "log2FoldChange"), "whole_tissue_padj": pick(w, g, "padj"),
                 "whole_tissue_composition_adjusted_padj": pick(wa, g, "padj"),
                 "fibroblast_normalised_log2FC": pick(f1, g, "log2FoldChange"), "fibroblast_normalised_padj": pick(f1, g, "padj"),
                 "fibroblast_normalised_mito_adjusted_padj": pick(f2, g, "padj"),
                 "median_log2CPM_diff_females": sx["F"], "median_log2CPM_diff_males": sx["M"],
                 "permutation_p_within_sex": tp.loc[g, "spp1_perm_p"] if g in tp.index else np.nan})
out(pd.DataFrame(rows), "09_spp1_matrix_genes", "spp1_vs_ntc_matrix_gene_summary.csv")

# ------------------------------------------------------------ 10 PRIMARY: SPP1 vs NTC within 5xFAD
for f, n in [("SPP1_vs_NTC_all_models.tsv.gz", "SPP1_vs_NTC_all_genes_all_models.csv.gz"), ("SPP1_vs_NTC_DE_summary.tsv", "SPP1_vs_NTC_DE_summary.csv"),
             ("SPP1_vs_NTC_key_genes.tsv", "SPP1_vs_NTC_key_genes.csv"), ("SPP1_vs_NTC_module_effects.tsv", "SPP1_vs_NTC_module_effects.csv"),
             ("SPP1_vs_NTC_composition_balance.tsv", "SPP1_vs_NTC_composition_balance.csv"), ("pathways_camera_SPP1_vs_NTC.tsv", "SPP1_vs_NTC_camera_pathways.csv")]:
    cp(os.path.join(T, f), "10_SPP1_vs_NTC_PRIMARY", n)
for f in glob.glob(os.path.join(T, "DE_SPP1_vs_NTC_*.tsv")):
    out(rd(f), "10_SPP1_vs_NTC_PRIMARY", os.path.basename(f).replace(".tsv", ".csv.gz"))

# ------------------------------------------------------------ figures
FIGS = {
    "A01_initial_PCA.png": "results/figures/02_pca.png",
    "A02_initial_composition_modules.png": "results/figures/04_composition_modules.png",
    "A03_spp1_vs_osteoclast.png": "results/figures/05_spp1_vs_osteoclast.png",
    "A04_initial_GSEA.png": "results/figures/08_gsea.png",
    "A05_genotype_signature_vs_technical.png": "results/figures/10_genotype_signature_vs_technical.png",
    "B01_fibroblast_reference_specificity_heatmap.png": "results/fibroblast/F1_reference_specificity_heatmap.png",
    "B02_fibroblast_deconvolution_heatmap.png": "results/fibroblast/F2_deconvolution.png",
    "B03_cell_representation_per_sample.png": "results/fibroblast/F2b_cell_representation_per_sample.png",
    "B04_spp1_fibroblast_ECM_score.png": "results/fibroblast/F7_spp1_fibroblast_ecm.png",
    "C01_genotype_vs_technical_616.png": "results/fibroblast_5xfad_spp1/G1_genotype_vs_technical.png",
    "C02_top_genotype_genes_616.png": "results/fibroblast_5xfad_spp1/G2_top_genotype_genes_616.png",
    "C03_spp1_on_5xfad_signature.png": "results/fibroblast_5xfad_spp1/G3_spp1_on_5xfad_signature.png",
}
for f in sorted(glob.glob(os.path.join(CA, "10_final_figures", "*.png"))):
    FIGS["D_" + os.path.basename(f)] = os.path.relpath(f, ROOT)
for f in sorted(glob.glob(os.path.join(CA, "11_tgfb_panel", "*.png"))):
    FIGS["E_" + os.path.basename(f)] = os.path.relpath(f, ROOT)
for f in sorted(glob.glob(os.path.join(CA, "12_spp1_vs_ntc", "*.png"))):
    FIGS["F_" + os.path.basename(f)] = os.path.relpath(f, ROOT)
for dst, src in FIGS.items():
    shutil.copy(os.path.join(ROOT, src), os.path.join(FIG, dst))

# ------------------------------------------------------------ manifest (every file, rows, columns)
man = []
for p in sorted(glob.glob(os.path.join(D, "**", "*.csv*"), recursive=True)):
    df = pd.read_csv(p, nrows=5)
    nrows = len(pd.read_csv(p, usecols=[0]))
    cols = list(df.columns)
    man.append({"file": os.path.relpath(p, H), "rows": nrows, "n_columns": len(cols),
                "columns": "; ".join(map(str, cols[:40])) + (f"; ... (+{len(cols) - 40} more)" if len(cols) > 40 else "")})
pd.DataFrame(man).to_csv(os.path.join(H, "FILE_MANIFEST.csv"), index=False)
print(f"{len(man)} data files, {len(FIGS)} figures")
