Extraction of features related to conservation of amino acid sequences in orthologous genes.
We thank Michael Peeri, Ph.D. student in Laboratory of Computational Systems and Synthetic Biology headed by Prof. Tamir Tuller, for preparing the input files.
Since the data is not belong to us we will only describe it and will not upload the input files.
Input / data files:
1.	Clw files - This kind of file contains the MSAs of genes in S. cerevisiae and their orthologous genes from up to 5 Orthologous Groups (OGs)  
(Saccharomycetaceae, Saccharomycetes, Ascomycota, Fungi, Eukaryota) in ClustalW format. OGs for S. cerevisiae genes were obtained from OrthoDB [1].
Amino-acid MSA for each OG was performed using Muscle with the option - maxiters 4.
The coding sequences for each genome were obtained from Ensembl FTP. For each amino acid position in each S. cerevisiae protein, we examined the MSAs for all OGs which
include it and analyzed the codons in that position.

2.  5310 csv files based on the MSA files (see data files of the first part of the conservation project).Each file contains the frequency of each of the amino acids 
that make up the gene and the frequency of each of the nucleotides in the codon. In addition, the maximal value found in any of the OGs in a specific position was recorded
for each amino acid and for each nucleotide position in the codon.

Featues based on the Clw files:
The first kind of feature was a calculation of the mean and median of the percentage of indels (gaps) in each position (column) within a multiple sequence alignment.
The second kind of feature was a calculation of entropy measures. The Shannon Entropy (SE) score was calculated for the amino acid distribution of each column within a MSA,
The calculation of the SE score (the features) was performed using two different methods: the first time we omitted the indels from the SE calculation, i.e. we took into 
account the 20 amino acids only, while in second time we recalculated the SE score including indels as one of the possible options, i.e. we considered 21 options overall in
this SE calculations. In each of the different methods for SE calculation, the mean and median values of SE were calculated based on the entropy values of all the gene's sequence.
In addition, we calculated the SE in the first and the last 100 sliding windows of the MSA, where each window contains 30 codons / positions. Afterwards, we added more features
that calculate the indices of the first window with a minimal and maximal entropy score, based on the values of the SE score in the first and the last 100 windows within the MSA. 

Featues based on the csv files:
The calculations of the genetic preservation measures (the features) were performed on two sets of creatures: one set is based on data from all creatures that have orthological
genes for  S. cerevisiae, and the other set contains creatures that are evolutionarily closest to S. cerevisiae. 
In each of the sets, the geometric mean and median values of the frequencies were calculated for both amino acids and nucleotides. Each time we focused on different areas of the 
sequence for each gene: the beginning, the end and the whole sequence of the gene. The beginning and end of the gene were defined as the first or last 50 amino acids respectively.


[1]        E. M. Zdobnov et al., “OrthoDB v9.1: Cataloging evolutionary and functional annotations for animal, fungal, plant, archaeal, bacterial and viral orthologs
           ,” Nucleic Acids Research, vol. 45, no. D1, pp. D744–D749, Jan. 2017, doi: 10.1093/nar/gkw1119.
