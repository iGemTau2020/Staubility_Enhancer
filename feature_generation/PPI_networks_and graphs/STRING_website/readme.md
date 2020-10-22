 extraction of features related to different types of interactions between genes from STRING website.
Files: 
Input files: 4932.protein.links.v11.0.txt, this is a table contains 2 columns of genes exist in S. cerevisiae and also a column of combined score which is actually the strength of interaction between 2 genes.
The code file: extract_ranks.m
Output file: rank_data_updated.csv, this is a table contains the 3 column: First column is gene name, while the other 2 columns are the extracted features-  see details below. In the output file there are 6574 genes of S. cerevisiae.

Features:
	Rank of gene- number of interactions that each gene creates with other genes. 
The rank of the gene is an important feature because genes that create a larger number of interactions are probably evolutionarily older, and also more conserved because if they 
are related to more genes, then most likely they tend to develop fewer mutations.
	Anoter feature is Weighted rank_i which takes into account not only the number of interactions in which each gene is involved, but also the weight of each one of its edges
  and creates a kind of average measure of the strength of interaction in which the gene is involved.
 
