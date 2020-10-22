Input / data files :SGA_ExE.txt, SGA_NxN.txt, SGA_ExN_NxE.txt tables.
Each table contains 2 columns for the pair of mutant genes that create interaction between them, one more column for genetic interaction score which is the strength of 
or negative interaction between the two genes (see details below), additional column of the fitness value of each gene and also a column of the double mutant fitness of interaction between pairs of genes.

After combining and sorting the 3 networks there we have collected information about 5707 genes in yeast.
The code file : Extract_rank_and_fitness.ipynb 
Output / The results : not_normalized_smf_dmf_ranks_feature.csv, this is a table contains one column of genes name and other column of the following features:

Features:
	Rank of gene - number of interactions that each gene creates with other genes.
	negative rank, Weighted rank positive,  Weighted rank negative, Absolute weighted rank, Smf - the average of the fitness of single mutant per gene and Weighted dmf.
