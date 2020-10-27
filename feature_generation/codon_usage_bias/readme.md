This script generates codon usage bias features. Detailed explanations can be found within the code. 


**Files:**
1. The script: "CUB_calculation-FINAL.ipynb" 
2. Output results: "CUB_results_final.xlsx"


**Input data files:**

1. file name: 4932-GPM_2012_09_Saccharomyces_cerevisiae.txt

  source: https://pax-db.org/species/4932. S.cerevisiae - Whole organism, SC (GPM,Oct,2012). 

  direct link: https://pax-db.org/dataset/4932/450/

2. file name: GSE75897_RiboZero_RPKMs.txt.gz 

  source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75897

  another link (same file, different source): ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75897/suppl/

3. file name: 1-s2.0-S0092867410003193-mmc2.xlsx

  source: https://www.sciencedirect.com/science/article/pii/S0092867410003193, Supplemental Information, Table S1. tRNA Copy Numbers, tRNA Abundance, and tAI of Codons in Various Organisms, Related to Figure 1.

  direct link: https://ars.els-cdn.com/content/image/1-s2.0-S0092867410003193-mmc2.xls

  The code also uses "orf_coding_all.fasta.gz" in the main directory of the feature generation project. 
