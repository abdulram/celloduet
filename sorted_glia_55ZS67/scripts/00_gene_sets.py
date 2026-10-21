"""00: curated lineage panels and state modules (written to data/modules.gmt for R)."""
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LINEAGE = {
 "L_Microglia": "P2ry12 Tmem119 Sall1 Cx3cr1 Hexb Siglech Gpr34 Olfml3",
 "L_BAM": "Mrc1 Lyve1 Cd163 Pf4 F13a1 Ms4a7",
 "L_Monocyte": "Ly6c2 Ccr2 Plac8 Chil3",
 "L_Neutrophil": "S100a8 S100a9 Ly6g Retnlg Ngp Camp",
 "L_DC": "Flt3 Xcr1 Clec9a Cd209a",
 "L_Tcell": "Cd3e Cd3g Cd2",
 "L_Bcell": "Cd79a Ms4a1 Cd19",
 "L_NK": "Nkg7 Ncr1 Klrb1c",
 "L_Astrocyte": "Aldh1l1 Slc1a3 Slc1a2 Gja1 Aqp4 Atp1b2 Sox9 Aldoc",
 "L_Neuron": "Snap25 Rbfox3 Syt1 Stmn2",
 "L_Oligodendrocyte": "Mbp Plp1 Mog Mobp",
 "L_OPC": "Pdgfra Cspg4 Olig2",
 "L_Endothelial": "Cldn5 Pecam1 Flt1 Ly6c1",
 "L_Pericyte_SMC": "Pdgfrb Rgs5 Kcnj8 Acta2",
 "L_Fibroblast": "Dcn Col1a2 Slc47a1",
 "L_Ependymal_CP": "Foxj1 Ttr",
 "L_Erythrocyte": "Hbb-bs Hba-a1",
}
STATE = {
 "M_Homeostatic": "P2ry12 Tmem119 Cx3cr1 Sall1",
 "M_DAM1": "Tyrobp Apoe B2m Ctsb Ctsd Lyz2",
 "M_DAM2": "Cst7 Itgax Clec7a Lpl Spp1 Gpnmb Axl Cd9 Igf1 Csf1",
 "M_DAM2_noSpp1": "Cst7 Itgax Clec7a Lpl Gpnmb Axl Cd9 Igf1 Csf1",
 "M_MGnD": "Apoe Axl Clec7a Itgax Spp1 Gpnmb Lgals3 Ccl2 Cst7 Lpl Itgb2 Cd9",
 "M_ARM": "Cst7 Itgax Clec7a Lpl Spp1 Gpnmb Dkk2 Cd74 H2-Aa H2-Ab1 Apoe Ank",
 "M_IRM_IFN": "Ifit1 Ifit2 Ifit3 Ifit3b Isg15 Irf7 Rtp4 Oasl2 Bst2 Usp18 Stat1 Ifi27l2a",
 "M_MHCII": "H2-Aa H2-Ab1 H2-Eb1 Cd74",
 "M_Proliferation": "Mki67 Top2a Ccnb2 Cdk1 Birc5 Cenpf",
 "M_Lipid": "Apoe Lpl Abca1 Ch25h Abcg1 Plin2 Fabp5 Lipa Soat1",
 "M_Lysosome_phago": "Ctsb Ctsd Ctsl Ctss Ctsz Lamp1 Cd68 Grn Hexa Npc2 Atp6v0d2",
 "M_Complement": "C1qa C1qb C1qc C3ar1 C5ar1 C4b",
 "M_Inflammatory": "Tnf Il1b Nfkbia Il1a Ccl2 Cxcl2 Ptgs2",
 "A_DAA": "Gfap Vim Serpina3n Ctss C4b Cd9 Apoe",
 "A_Neurotoxic_A1": "C3 Serping1 H2-D1 Gbp2 Fkbp5 Ugt1a1 Psmb8 Srgn Amigo2",
 "A_PanReactive": "Lcn2 Steap4 S1pr3 Timp1 Hspb1 Cxcl10 Cd44 Osmr Cp Serpina3n Aspg Vim Gfap",
 "A_A2": "Clcf1 Tgm1 Ptx3 S100a10 Sphk1 Cd109 Ptgs2 Emp1 Slc10a6 Tm4sf1 B3gnt5 Cd14",
 "A_IFN": "Ifit1 Ifit3 Isg15 Irf7 Rtp4 Oasl2 Bst2 Stat1 Gbp2 Iigp1 Igtp",
 "A_Synapse": "Sparcl1 Gpc4 Gpc6 Thbs1 Thbs2 Sparc Hevin",
 "A_Glut_K": "Slc1a2 Slc1a3 Kcnj10 Glul",
 "A_Cholesterol": "Hmgcr Hmgcs1 Sqle Fdft1 Srebf2 Idi1 Cyp51 Dhcr24 Ldlr Insig1",
 "A_Endfoot": "Aqp4 Gja1",
 "X_Dissociation_IEG": "Fos Fosb Jun Junb Egr1 Egr2 Atf3 Ier2 Ier3 Dusp1 Hspa1a Hspa1b Klf2 Klf4 Rgs1 Ccl3 Ccl4 Nfkbia Zfp36",
 "X_Spp1_receptors": "Cd44 Itgav Itgb1 Itgb3 Itgb5",
}
ALL = {**LINEAGE, **STATE}
if __name__ == "__main__":
    with open(os.path.join(ROOT, "data", "modules.gmt"), "w") as f:
        for k, v in ALL.items():
            f.write("\t".join([k, "curated"] + v.split()) + "\n")
    print(len(ALL), "sets written")
