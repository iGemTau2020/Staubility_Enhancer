# Table of Contents
1. [Input files](#input-files)
2. [The code file ](#the-code-file)
3. [Output](#output)
4. [Features](#features)

## Input files
Input files: SGA_ExE.txt, SGA_NxN.txt, SGA_ExN_NxE.txt tables.  
Please note that all the files here are sampled due to their size and can be used only as sample data.  
Each table contains 2 columns for each pair of interacting mutant genes, one more column for genetic interaction score which is the strength of their mutual interaction (see details below), an additional column of the fitness value of each gene and also a column of the double mutant fitness of interaction between pairs of genes.
After combining and sorting the 3 networks there we have collected information about 5707 genes in yeast.

## The code file 
Extract_rank_and_fitness.ipynb 

## Output 
not_normalized_smf_dmf_ranks_feature.csv, this table contains one column of genes name and other columns with the following features:

## Features
Rank of gene - number of interactions that each gene creates with other genes.
negative rank, Weighted rank positive,  Weighted rank negative, Absolute weighted rank, Smf - the mean fitness of single mutant per gene and Weighted dmf.
