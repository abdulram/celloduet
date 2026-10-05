"""03b - Define pathway / module gene sets used in steps 04-08 (written to tables/gene_sets.tsv).

Sources
  MSIGDB_HALLMARK_*  MSigDB Hallmark v7.0 (human symbols; mapped to mouse by case-insensitive symbol match).
                     File: data/genesets/h.all.v7.0.symbols.gmt (mirrored from the GSEApy test data on GitHub).
  CURATED_*          Mouse gene lists curated for this analysis from canonical markers / well-established
                     literature programs (e.g. DAM: Keren-Shaul et al. 2017 Cell; BAM: Van Hove et al. 2019 Nat
                     Neurosci; meningeal lymphatics: Louveau 2015 / Da Mesquita 2018 Nature; dural immune trafficking:
                     Rustenhoven et al. 2021 Cell). They are small, hand-picked and should be read as modules, not
                     exhaustive pathways.
  TECHNICAL_*        Negative-control "technical sentinel" sets (mitochondrial-encoded, cytosolic ribosomal):
                     expected to track RNA quality / library prep, not 5xFAD biology.
Lineage-identity genes used for composition signatures (step 02) are deliberately NOT reused in the state modules
used by the cell-state analysis (step 07).
"""
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CA = os.path.dirname(HERE)
ROOT = os.path.dirname(CA)

CURATED = {
    "CURATED_COMPLEMENT_CLASSICAL": ["C1qa", "C1qb", "C1qc", "C1ra", "C1s1", "C2", "C3", "C4b", "C5ar1", "C3ar1", "Cfh", "Serping1", "Cd55", "Itgam"],
    "CURATED_TYPE_I_IFN_ISG": ["Ifit1", "Ifit2", "Ifit3", "Ifit3b", "Isg15", "Irf7", "Oasl2", "Rsad2", "Usp18", "Stat1", "Bst2", "Ifi27l2a", "Mx1", "Oas1a", "Rtp4", "Xaf1", "Ifi44"],
    "CURATED_MHCII_ANTIGEN_PRESENTATION": ["H2-Aa", "H2-Ab1", "H2-Eb1", "Cd74", "H2-DMa", "H2-DMb1", "Ciita", "Cd86", "Cd80", "Ctss"],
    "CURATED_CHEMOKINES": ["Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl7", "Ccl8", "Ccl12", "Ccl19", "Ccl21a", "Cxcl1", "Cxcl2", "Cxcl9", "Cxcl10", "Cxcl12", "Cxcl13", "Cx3cl1"],
    "CURATED_BAM_PROGRAM": ["Mrc1", "Cd163", "Lyve1", "Pf4", "F13a1", "Ms4a7", "Cbr2", "Folr2", "Stab1", "Dab2", "Gas6"],
    "CURATED_MONOCYTE_RECRUITMENT": ["Ccl2", "Ccl7", "Ccl12", "Ccr2", "Ly6c2", "Sell", "Itga4", "Vcam1", "Icam1", "Cx3cr1"],
    "CURATED_T_CELL_ACTIVATION": ["Cd69", "Il2ra", "Icos", "Ctla4", "Pdcd1", "Cd40lg", "Tnfrsf4", "Tnfrsf9", "Il7r", "Cd44"],
    "CURATED_T_CELL_CYTOTOXICITY": ["Gzma", "Gzmb", "Gzmk", "Prf1", "Ifng", "Nkg7", "Ccl5", "Fasl", "Cd8a", "Cd8b1"],
    "CURATED_B_PLASMA": ["Cd79a", "Cd79b", "Ms4a1", "Cd19", "Jchain", "Mzb1", "Igha", "Ighm", "Xbp1", "Prdm1", "Tnfrsf17"],
    "CURATED_PHAGOCYTOSIS_LYSOSOME": ["Cd68", "Lamp1", "Lamp2", "Ctsb", "Ctsd", "Ctsl", "Mertk", "Axl", "Gas6", "Msr1", "Marco", "Fcgr1", "Fcgr3", "Lgals3"],
    "CURATED_LIPID_HANDLING": ["Apoe", "Lpl", "Abca1", "Abcg1", "Lipa", "Cd36", "Trem2", "Plin2", "Fabp5", "Soat1", "Npc2"],
    "CURATED_DAM_KerenShaul2017": ["Cst7", "Itgax", "Clec7a", "Lpl", "Spp1", "Gpnmb", "Cd9", "Apoe", "Trem2", "Tyrobp", "Axl", "Ank", "Csf1", "Ctsb", "Ctsd", "Ctsl", "Cd63", "Lgals3", "Igf1", "Ch25h", "Fabp5", "Lilrb4a"],
    "CURATED_ECM_REMODELING": ["Mmp2", "Mmp3", "Mmp14", "Timp1", "Timp2", "Timp3", "Lox", "Loxl1", "Loxl2", "Col1a1", "Col3a1", "Col4a1", "Col5a1", "Fn1", "Sparc", "Postn", "Tnc"],
    "CURATED_FIBROBLAST_ACTIVATION": ["Postn", "Cthrc1", "Acta2", "Fn1", "Tnc", "Timp1", "Serpine1", "Thbs1", "Ctgf", "Lox", "Col12a1", "Fap"],
    "CURATED_INFLAMMATORY_FIBROBLAST": ["Il6", "Ccl2", "Cxcl1", "Cxcl12", "Ccl11", "Il33", "Lif", "C3", "Cfb", "Saa3", "Lcn2", "Ptgs2"],
    "CURATED_VASCULAR_INFLAMMATION_ADHESION": ["Sele", "Selp", "Vcam1", "Icam1", "Madcam1", "Ackr1", "Lrg1", "Plvap", "Cd34", "Il1r1"],
    "CURATED_ANGIOGENESIS_ENDOTHELIAL": ["Vegfa", "Kdr", "Flt1", "Angpt2", "Tek", "Dll4", "Esm1", "Apln", "Aplnr", "Nos3", "Notch1"],
    "CURATED_LYMPHATIC_ENDOTHELIAL": ["Prox1", "Flt4", "Lyve1", "Ccl21a", "Pdpn", "Vegfc", "Mmrn1", "Reln", "Foxc2", "Itga9"],
    "CURATED_MENINGEAL_IMMUNE_TRAFFICKING": ["Ccl19", "Ccl21a", "Ccr7", "S1pr1", "Cxcl12", "Cxcr4", "Madcam1", "Itga4", "Itgb7", "Sell", "Glycam1"],
    "CURATED_INNATE_IMMUNE_ACTIVATION": ["Tlr2", "Tlr4", "Myd88", "Nfkbia", "Tnf", "Il1b", "Il1a", "Nlrp3", "Casp1", "Cd14", "Tnfaip3"],
    "TECHNICAL_MITO_ENCODED": ["mt-Nd1", "mt-Nd2", "mt-Co1", "mt-Co2", "mt-Atp8", "mt-Atp6", "mt-Co3", "mt-Nd3", "mt-Nd4l", "mt-Nd4", "mt-Nd5", "mt-Nd6", "mt-Cytb"],
}
RIBO = ["Rpl3", "Rpl4", "Rpl5", "Rpl7", "Rpl8", "Rpl10", "Rpl11", "Rpl13", "Rpl18", "Rpl19", "Rpl23", "Rpl27a", "Rpl32", "Rpl37",
        "Rps2", "Rps3", "Rps5", "Rps6", "Rps8", "Rps9", "Rps11", "Rps12", "Rps14", "Rps16", "Rps18", "Rps19", "Rps23", "Rps27a"]
CURATED["TECHNICAL_CYTOSOLIC_RIBOSOME"] = RIBO

HALLMARK_KEEP = ["HALLMARK_COMPLEMENT", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                 "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                 "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_ANGIOGENESIS",
                 "HALLMARK_COAGULATION", "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
                 "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_HYPOXIA", "HALLMARK_APOPTOSIS", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"]

raw = pd.read_csv(os.path.join(ROOT, "data", "VT3YXB-expression-matrix.tsv"), sep="\t", usecols=["gene_id", "gene_name"])
raw = raw.dropna()
up = pd.Series(raw.gene_id.values, index=raw.gene_name.str.upper().values)
up = up[~up.index.duplicated()]
mouse = pd.Series(raw.gene_id.values, index=raw.gene_name.values)
mouse = mouse[~mouse.index.duplicated()]
rows = []
with open(os.path.join(ROOT, "data", "genesets", "h.all.v7.0.symbols.gmt")) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[0] in HALLMARK_KEEP:
            for g in f[2:]:
                if g.upper() in up.index:
                    rows.append({"set": "MSIGDB_" + f[0], "gene_id": up[g.upper()], "source": "MSigDB Hallmark v7.0 (human->mouse by symbol)"})
for s, gl in CURATED.items():
    for g in gl:
        if g in mouse.index:
            rows.append({"set": s, "gene_id": mouse[g], "source": "curated (see script docstring)" if s.startswith("CURATED") else "technical sentinel"})
gs = pd.DataFrame(rows).drop_duplicates(["set", "gene_id"])
gs["gene_name"] = raw.set_index("gene_id").gene_name.reindex(gs.gene_id).values
gs.to_csv(os.path.join(CA, "tables", "gene_sets.tsv"), sep="\t", index=False)
print(gs.groupby("set").size().to_string())
