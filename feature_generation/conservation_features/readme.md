# Table of Contents
1. [Input](#input-files)
2. [Output](#output-files)
3. [Features](#Features)

In this directory, we generate features related to conservation of amino acid sequences in orthologous genes.
We thank Michael Peeri, a Ph.D. student in Laboratory of Computational Systems and Synthetic Biology headed by Prof. Tamir Tuller, for preparing the input files.
Since the data does not belong to us, we cannot upload the input files, and will only describe them.

## Input files
1. Clw files - These files contais the Multiple Sequence Alignments (MSAs) of genes in S. cerevisiae and their orthologous genes from up to 5 Orthologous Groups (OGs)  
(Saccharomycetaceae, Saccharomycetes, Ascomycota, Fungi, Eukaryota) in ClustalW format. OGs for S. cerevisiae genes were obtained from OrthoDB [1].
Amino-acid MSA for each OG was performed using Muscle with the option - maxiters 4.
The coding sequences for each genome were obtained from Ensembl FTP. For each amino acid position in each S. cerevisiae protein, we examined the MSAs for all OGs which
include it and analyzed the codons in that position.

2. 5310 csv files based on the MSA files (see data files of the first part of the conservation project). Each file contains the frequency of each of the amino acids
that make up the gene and the frequency of each of the nucleotides in the codon. In addition, the maximal value found in any of the OGs in a specific position was recorded
for each amino acid and for each nucleotide position in the codon.

3. ensembelsld.pkl - a pickle file used by python for saving objects. We used this file in order to store a mapping from OrthoDB ids to Ensemble ids.
Please note that this file is sampled due to the size of the original data.

## Output files
1. entropy features gene names converted.csv. This table contains 417 columns: the first column is the gene name, while the other 416 columns are the extracted features-see details below. In the output file there are 4691 genes of S. cerevisiae.

2. conservation features gene names converted.csv, this table contains 25 columns: the first column contains the gene name, while the other 24 columns are the extracted features-  see details below. In the output file there are 5310 genes of S. cerevisiae.

## Features
### Features based on the Clw files
The first feature family is a calculation of the mean and median percentage of indels (gaps) in each position (column) within a multiple sequence alignment.

The second feature family is a calculation of entropy measures. The Shannon Entropy (SE) score was calculated for the amino acid distribution of each column within a MSA.
This calculation of the SE score was performed using two different methods: the first time we omitted the indels from the SE calculation, i.e. we took into 
account the 20 possible amino acids only, while in the second time we recalculated the SE score including indels as one of the possible options, i.e. we considered 21 options overall. 
In each of the different methods for SE calculation, the mean and median values of SE were calculated based on the entropy values of all the gene's sequence.
In addition, we calculated the SE in the first and the last 100 sliding windows of the MSA, where each window contains 30 codons / positions. Afterwards, we added more features
that calculate the indices of the first window with a minimal and maximal entropy score, based on the values of the SE score in the first and the last 100 windows within the MSA. 

### Features based on the csv files
The calculation of the genetic preservation measures was performed on two sets of organisms's genomes: one set is the data from all creatures that have orthological
genes for  S. cerevisiae, and the other set is only the data of organisms closest evolutionarily to S. cerevisiae. 
In each of the sets, the geometric mean and median values of the frequencies were calculated for both amino acids and nucleotides. Each time, we focused on different areas of the sequence for each gene: the beginning, the end and the whole sequence of the gene. The beginning and end of the gene were defined as the first or last 50 amino acids respectively.


[1] E. M. Zdobnov et al., “OrthoDB v9.1: Cataloging evolutionary and functional annotations for animal, fungal, plant, archaeal, bacterial and viral orthologs
,” Nucleic Acids Research, vol. 45, no. D1, pp. D744–D749, Jan. 2017, doi: 10.1093/nar/gkw1119.
