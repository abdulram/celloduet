# Figure → data map
| Figure | Shows | Data | Script |
|---|---|---|---|
| figures/fig04_pca_per_fraction.png | PCA (VST, top 2,000 genes) per fraction | 04_exploration/pca_scores_*.tsv, pc_covariate_correlations.tsv | scripts/05_pca_variance.py |
| figures/fig05_spp1_and_candidates_CD11b.png | Spp1 and co-regulated microglial genes, per sample | 05_DE_CD11b/DE_Spp1_vs_NTC_CD11b_annotated.tsv | scripts/09_figures_key_genes.py |
| figures/fig06_DAM_controls_CD11b.png | Core DAM/homeostatic genes (unchanged) | same | same |
| figures/fig07_spp1_axis_ACSA2.png | Spp1, microglial carry-over markers, Spp1 receptors, reactive genes in ACSA-2 | 06_DE_ACSA2/... | same |
| figures/fig08_glucocorticoid_both.png | Glucocorticoid-response genes in both fractions | raw CPM | same |
| figures/fig09_carryover_CD11b.png | Carry-over genes in CD11b (oligo/astro/neutrophil/BAM) | 02_purity_celltyping/purity_carryover_per_sample.tsv | same |
| figures/fig10_pathology_programs_forest.png | 20 pathology programs, SPP1 vs NTC % change with 95% CI | 09_pathology/pathology_programs_SPP1_vs_NTC.tsv | scripts/10_pathology_programs.py |
