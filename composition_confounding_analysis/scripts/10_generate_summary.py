"""10 - Assemble final figures, machine-readable summary (analysis_summary.tsv) and package versions.

Reads only upstream outputs (tables/ and figure folders). analysis_summary.tsv rows are generated from the
tables so they stay in sync with the analysis; the narrative README / EXECUTIVE_SUMMARY are written by hand
from these same tables.
"""
import os
import platform
import shutil
import subprocess

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CA = os.path.dirname(HERE)
TAB = os.path.join(CA, "tables")
FF = os.path.join(CA, "10_final_figures")

# ---------------------------------------------------------------- final figures
FIGS = {
    "Fig01_transcriptome_PCA_by_genotype.png": "01_QC/PCA_by_genotype.png",
    "Fig02_cell_signature_heatmap.png": "02_cell_composition/signature_heatmap.png",
    "Fig02b_lineage_marker_heatmap.png": "02_cell_composition/lineage_marker_heatmap.png",
    "Fig03_signature_scores_by_genotype.png": "02_cell_composition/signature_scores_by_group.png",
    "Fig04_composition_PCA.png": "03_composition_PCA/Fig4_composition_PCA.png",
    "Fig05_txPC_vs_composition.png": "03_composition_PCA/Fig5_txPC_vs_composition.png",
    "Fig06_common_support.png": "03_composition_PCA/Fig6_common_support.png",
    "Fig07_baseline_vs_adjusted_LFC.png": "05_composition_adjusted_DE/Fig7_baseline_vs_adjusted_LFC.png",
    "Fig08_pathway_effects_by_model.png": "08_pathways/Fig8_pathway_effects_by_model.png",
    "Fig09_state_vs_abundance.png": "07_cell_state_analysis/Fig9_state_vs_abundance.png",
    "Fig10a_common_support_selection.png": "06_common_support_matching/Fig10a_common_support_selection.png",
    "Fig10b_common_support_LFC.png": "06_common_support_matching/Fig10b_common_support_LFC.png",
    "Fig10c_technical_axis_animal_position.png": "06_common_support_matching/Fig10c_animal_position_test.png",
    "Fig11_pathology_NOT_AVAILABLE.png": "09_pathology_associations/Fig11_pathology_NOT_AVAILABLE.png",
    "Fig12_pathway_robustness_summary.png": "08_pathways/Fig12_pathway_robustness_summary.png",
}
for dst, src in FIGS.items():
    shutil.copy(os.path.join(CA, src), os.path.join(FF, dst))

# ---------------------------------------------------------------- machine-readable summary
gc = pd.read_csv(os.path.join(TAB, "genotype_vs_composition.tsv"), sep="\t").set_index("measure")
adj = pd.read_csv(os.path.join(TAB, "adjusted_models_summary.tsv"), sep="\t").set_index("model")
cs = pd.read_csv(os.path.join(TAB, "common_support_DE_summary.tsv"), sep="\t").set_index("cohort")
pos = pd.read_csv(os.path.join(TAB, "technical_axis_animal_position_test.tsv"), sep="\t").set_index("animal")
pw = pd.read_csv(os.path.join(TAB, "pathway_robustness_matrix.tsv"), sep="\t")
st = pd.read_csv(os.path.join(TAB, "cell_state_vs_abundance.tsv"), sep="\t")
base = pd.read_csv(os.path.join(TAB, "DE_baseline_5xFAD_vs_WT.tsv"), sep="\t")
val = pd.read_csv(os.path.join(TAB, "signature_validation.tsv"), sep="\t")
pos.index = pos.index.astype(str)

rows = []
nb = int((base.padj < 0.05).sum())
rows.append(dict(finding="Genome-wide baseline 5xFAD vs WT signature (6.8k genes) tracks the technical/processing axis",
                 category="global DE", cell_type_or_pathway="all genes",
                 baseline_effect=f"{nb} genes padj<0.05",
                 adjusted_effect=f"+CompPC1: {adj.loc['comp1','adjusted_sig']} sig, median retention {adj.loc['comp1','median_retention']}; +%mito: {adj.loc['tech','adjusted_sig']} sig, retention {adj.loc['tech','median_retention']}",
                 matched_effect=f"matched 2v2: direction concordance {cs.loc['matched','direction_concordance_baseline_sig']}, 0 sig; WT616 position {pos.loc['616','median_position']:.2f} (other WT {pos.drop('616').median_position.min():.2f}..{pos.drop('616').median_position.max():.2f})",
                 pathology_association="not available", confidence_tier="Tier 4",
                 interpretation="WT 616 (only WT with 5xFAD-like library profile) looks 5xFAD on 84% of genes; genotype not separable from processing batch (or 616 mis-genotyped)",
                 supporting_figure="Fig07; Fig10c", supporting_table="tables/adjusted_models_summary.tsv; tables/technical_axis_animal_position_test.tsv"))
for m, desc in [("Composition_PC1", "dura-purity axis (fibroblast/endothelial vs olfactory/CNS/muscle)"), ("Composition_PC3", "bone axis (osteoblast/osteoclast)")]:
    r = gc.loc[m]
    rows.append(dict(finding=f"{m} differs between genotypes ({desc})", category="composition", cell_type_or_pathway=m,
                     baseline_effect=f"Hedges g {r.hedges_g:.2f}, AUC {r.AUC_5xFAD_gt_WT:.2f}, p {r.MWU_p:.3f}", adjusted_effect="n/a",
                     matched_effect="common support poor", pathology_association="not available", confidence_tier="Tier 3/4",
                     interpretation="cannot distinguish disease remodeling from dissection differences; poor overlap between genotypes",
                     supporting_figure="Fig04; Fig06", supporting_table="tables/genotype_vs_composition.tsv"))
for s_ in ["blood_endothelial", "bone_osteoblast", "dendritic_cell", "B_cell"]:
    r = gc.loc[s_]
    rows.append(dict(finding=f"{s_} signature differs by genotype (nominal)", category="composition", cell_type_or_pathway=s_,
                     baseline_effect=f"g {r.hedges_g:.2f}, AUC {r.AUC_5xFAD_gt_WT:.2f}, p {r.MWU_p:.3f}", adjusted_effect="n/a", matched_effect="n/a",
                     pathology_association="not available", confidence_tier="Tier 3",
                     interpretation="abundance difference; biology vs dissection indistinguishable", supporting_figure="Fig03; Fig06",
                     supporting_table="tables/genotype_vs_composition.tsv"))
for _, r in pw[pw.baseline_p < 0.1].iterrows():
    rows.append(dict(finding=f"{r.pathway} module", category="pathway", cell_type_or_pathway=r.pathway,
                     baseline_effect=f"{r.baseline_effect_SD:.2f} SD (p {r.baseline_p:.3f})",
                     adjusted_effect=f"+CompPC1 {r.comp_adjusted_effect_SD:.2f} (p {r.comp_adjusted_p:.2f}); +PC1+PC3 {r.comp13_adjusted_effect_SD:.2f}; +%mito {r.tech_adjusted_effect_SD:.2f} (p {r.tech_adjusted_p:.2f})",
                     matched_effect=f"{r.matched_effect_SD:.2f} SD (2v2); WT616 position {r.position_616:.2f}",
                     pathology_association="not available", confidence_tier=r.confidence_tier,
                     interpretation={"technical": "negative control: tracks RNA quality / library prep"}.get(r.associated_lineage, ""),
                     supporting_figure="Fig08; Fig12", supporting_table="tables/pathway_robustness_matrix.tsv"))
for _, r in st[st.abundance_signature_validated & (st.perm_p_given_abundance < 0.1)].iterrows():
    rows.append(dict(finding=f"{r.state_module} in {r.lineage} beyond abundance", category="cell state", cell_type_or_pathway=f"{r.lineage}:{r.state_module}",
                     baseline_effect=f"unadj beta {r.genotype_beta_unadjusted:.2f} (p {r.p_unadjusted:.2f})",
                     adjusted_effect=f"given abundance {r.genotype_beta_given_abundance:.2f} (perm p {r.perm_p_given_abundance:.2f}); +%mito {r.genotype_beta_given_abundance_and_mito:.2f} (p {r.p_given_abundance_and_mito:.2f})",
                     matched_effect="n/a", pathology_association="not available", confidence_tier="Tier 4",
                     interpretation="nominal only (no multiple-testing correction across 15 tests) and lost after technical adjustment",
                     supporting_figure="Fig09", supporting_table="tables/cell_state_vs_abundance.tsv"))
for _, r in val[~val.validated].iterrows():
    rows.append(dict(finding=f"{r.signature} cannot be quantified (markers incoherent)", category="limitation", cell_type_or_pathway=r.signature,
                     baseline_effect="n/a", adjusted_effect="n/a", matched_effect="n/a", pathology_association="not available",
                     confidence_tier="Tier 4", interpretation=f"median item-rest r {r.median_item_rest_r:.2f}; abundance and state of this lineage not estimable",
                     supporting_figure="Fig02b", supporting_table="tables/signature_validation.tsv"))
summary = pd.DataFrame(rows)
summary.to_csv(os.path.join(CA, "analysis_summary.tsv"), sep="\t", index=False)
print(summary[["finding", "confidence_tier"]].to_string())

# ---------------------------------------------------------------- versions
import importlib
vers = {"python": platform.python_version()}
for p in ["numpy", "pandas", "scipy", "statsmodels", "matplotlib", "tiledbsoma"]:
    try:
        vers[p] = importlib.import_module(p).__version__
    except Exception:
        vers[p] = "n/a"
r = subprocess.run(["Rscript", "-e", 'cat(R.version.string, "\\n"); for (p in c("DESeq2","limma","MatchIt","car","ggplot2","pheatmap")) cat(p, as.character(packageVersion(p)), "\\n")'],
                   capture_output=True, text=True)
with open(os.path.join(CA, "logs", "package_versions.txt"), "w") as fh:
    fh.write("\n".join(f"{k} {v}" for k, v in vers.items()) + "\n" + r.stdout)
print(open(os.path.join(CA, "logs", "package_versions.txt")).read())
