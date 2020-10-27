Target Gene Features

We calculated distances features between the target gene and each of the yeast genes, and from the mean of the 5 strongest connections of each yeast gene. All features described below are calculated based on profiles we created for each subcategory. For each profile, we applied several distance functions between it and the target gene’s profile and the mean of the 5 strongest connections' profiles. 

The 5 strongest connection for each yeast gene were calculated using the following function:
max_5_ind_neighbours
For each gene, we used the string PPI to create a vector containing its 5 strongest connections (5 highest scores from the string PPI).
inputs - 4932.protein.links.v11.0.txt, seq_orf.csv
output - max_5_ind.csv

All distances features were calculated using the final_func functions and include:
a.	Codon relative frequency 
For each gene, we calculated the relative frequency of each codon – the number of appearances of a codon divided by the number of appearances of the amino acid that the codon encodes. The result is a vector of length 64 that contains each relative frequency.
This vector was also calculated for the target gene.
We calculated the following distances between the target gene’s vector and each yeast gene’s vector:
i.	Euclidian distance
ii.	L1 distance 
iii.	Spearman correlation
iv.	Pearson correlation 
v.	KS distance
vi.	Cartesian product

b.	Amino acid frequency
For each gene, we calculated the relative frequency of each amino acid – the number of appearances of an amino acid divided by the length of the amino-acid sequence. 
This vector was also calculated for the target gene.
We calculated the following distances between the target gene’s vector and each yeast gene’s vector:
i.	Euclidian distance
ii.	L1 distance 
iii.	Spearman correlation
iv.	Pearson correlation 
v.	KS distance
vi.	Cartesian product

c.	GC content profiles
We created a GC content profile based on 10 windows of 30 nucleotides (nt) in a 300 nt window (the first 300 nt). For each sub-window, the GC content was calculated. 
As for the target gene – the same calculation was done using the 300 nt from the end of the sequence.
We calculated the following distances between the target gene’s vector and each yeast gene’s vector:
- Euclidian distance
- L1 distance 
- Spearman correlation
- Pearson correlation 
- KS distance
- Cartesian product
* for genes with a sequence length less than 300, the rest of the profile vector was filled with the average of the existing GC content in windows of 30 nt. For example – for a gene with sequence length of 100 nt – a 1X3 vector was created which contains 3 GC content results of 3 windows of 30 nt. The rest of the vector was filled with the average of this 1X3 vector in order to create a 1X10 vector.

d.	String PPI
For each gene, we used the string PPI to generate a vector containing its 5 strongest connections (5 highest scores from the string PPI).
For each gene, we calculated and compared the mean codon relative frequency, the mean amino acid frequency and the mean GC content profile For it and its strongest connections, using the following distance measures: 
- Euclidian distance
- L1 distance 
- Spearman correlation
- Pearson correlation 
- KS distance
- Cartesian product
* for genes with no connections the distance is NaN.
* if a gene has less than 5 connections, then the mean is calculated based on the connections it has available.

e.	mRNA folding energy
We created mRNA folding energy profiles using the function RNAfold from the package ViennaRNA.  These profiles were calculated for each conjugated gene, using 12 windows of 40 nucleotides (nt) with a 10 nt shift, resulting in a total window of the first 150 nt window. For each sub-window, the folding energy was calculated. 
The target gene had the same calculation using the 150 nt from the end of the sequence.
We calculated the following distances between the target gene’s vector and each yeast gene’s vector:
- Euclidian distance
- L1 distance 
- Spearman correlation
- Pearson correlation 
- KS distance
- Cartesian product
* for genes with a sequence length less than 150, the rest of the profile vector was filled with the mean on the existing folding energy value. For example – for a gene with sequence length of 100 nt, a vector was created containing 7 folding energy values resulting from 7 windows of 40 nt with a 10 nt shift each time. The rest of the vector was filled with the average of this 1X7 vector in order to create a 1X12 vector.

One more feature is the longest common substring (between the target gene and each yeast gene) which can also be found in the final_func function.

Inputs - target gene sequence, "max_5_ind.csv, seq.csv
