In this directory, features related to codon usage bias are generated. 

**Code:**
The script: "CUB_calculation-FINAL.ipynb". This script generates codon usage bias related features. 

**Input data files:**

1. file name: "4932-GPM_2012_09_Saccharomyces_cerevisiae.txt"

  source: https://pax-db.org/species/4932. S.cerevisiae - Whole organism, SC (GPM,Oct,2012). 

  direct link: https://pax-db.org/dataset/4932/450/
  
  It contains the protein abundance of approximately 90% of the genes in yeast. This data is used as a seperate feature, and also as a way to select the top 2% highly expressed genes to serve as reference set for CAI and RCA calculation (see below)

2. file name: "GSE75897_RiboZero_RPKMs.txt.gz" 

  source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75897

  another link (same file, different source): ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75897/suppl/
  
  It contains the mRNA abundance of the genes in yeast. This data is used as a seperate feature, and also for calculation of NTE feature (see below). 

3. file name: "1-s2.0-S0092867410003193-mmc2.xlsx"

  source: https://www.sciencedirect.com/science/article/pii/S0092867410003193, Supplemental Information, Table S1. tRNA Copy Numbers, tRNA Abundance, and tAI of Codons in Various Organisms, Related to Figure 1.
  
  It containts the tAI weights of each codon, adapted for yeast. 

  direct link: https://ars.els-cdn.com/content/image/1-s2.0-S0092867410003193-mmc2.xls

4. file name: "orf_coding_all.fasta.gz" in the main directory of the feature generation project. 
It contains the sequences of the yeast's genes. Since codon usage bias refers to these sequences, this file is highly relevant. 

**output:** 
Output results: "CUB_results_final.xlsx"

The calculated features:
1. gene length
2. mRNA and protein levels
3. ENC measure - effective number of codons (Wright, 1990), measures the deviation from uniform use of synonymous codons. 
4. CAI measure - codon adaptation index. It compares the codon usage of a given gene to a reference set. This feature is calculated once with the whole genome as reference, and once only with highly-expressed genes as reference. 
5. tAI - tRNA adaptation index (the adaptation to the tRNA pool). 
6. NTE - normalized translation efficiency, defined as the supply/demand ratio of tRNA molecules. 
7. RCBS - relative codon usage bias. 
8. RCA -  relative codon adaptation. It compares the nucleotide usage in a given gene to a reference set. This feature is calculated once with the whole genome as reference, and once only with highly-expressed genes as reference. 

The equations for the calculation of these measures are presented within the script. 

