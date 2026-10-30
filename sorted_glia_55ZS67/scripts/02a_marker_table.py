"""02a: CPM of lineage, state and Spp1-axis genes per sample (raw look before modelling)."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from common import *
ann, cnt = load_matrix()
cpm = cnt / cnt.sum() * 1e6
cpm.index = ann.gene_name
cpm = cpm.groupby(level=0).sum()
genes = ["Itgam","Ptprc","P2ry12","Tmem119","Sall1","Cx3cr1","Hexb","Fcrls","Siglech","Gpr34","Olfml3",
 "Mrc1","Lyve1","Cd163","Pf4","F13a1","Ms4a7","Ly6c2","Ccr2","Plac8","Chil3","S100a8","S100a9","Ly6g","Retnlg","Ngp","Camp",
 "Flt3","Xcr1","Clec9a","Cd209a","Cd3e","Cd3g","Cd2","Cd79a","Ms4a1","Cd19","Nkg7","Ncr1","Klrb1c",
 "Aldh1l1","Slc1a3","Slc1a2","Gja1","Aqp4","Atp1b2","Sox9","Gfap","Aldoc","Snap25","Rbfox3","Syt1","Stmn2",
 "Mbp","Plp1","Mog","Mobp","Pdgfra","Cspg4","Olig2","Cldn5","Pecam1","Flt1","Ly6c1","Pdgfrb","Rgs5","Kcnj8","Acta2",
 "Dcn","Col1a2","Slc47a1","Foxj1","Ttr","Hbb-bs","Hba-a1",
 "Spp1","Cd44","Itgav","Itgb1","Itgb3","Itgb5","Cst7","Itgax","Clec7a","Lpl","Gpnmb","Axl","Cd9","Igf1","Csf1","Apoe","Trem2","Tyrobp","Lgals3","Ccl3","Ccl4",
 "Fos","Junb","Egr1","Hspa1a","Ifit3","Irf7","Isg15","Cd74","H2-Aa","Mki67","Serpina3n","C4b","Vim","Lcn2","C3","App","Psen1","Thy1"]
t = cpm.reindex(genes).round(1)
t.columns = [f"{c.split('_')[1]}" for c in t.columns]
t.to_csv(os.path.join(ROOT, "02_purity_celltyping", "marker_cpm_raw.tsv"), sep="\t")
pd.set_option("display.width", 250)
print("cols: 1-3 NTC CD11b | 4-6 Spp1 CD11b | 7-9 NTC ACSA2 | 10-12 Spp1 ACSA2")
print(t.to_string())
