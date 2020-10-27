extraction of features related to different types of interactions between genes from STRING website.

# Table of Contents
1. [Input](#input-files)
2. [Code](#code-file)
3. [Output](#output-file)
4. [Features](#Features)

## Input files
4932.protein.links.v11.0.txt, this table contains 2 columns of every gene in S. cerevisiae and a column of their interaction score.

## Code file
extract_ranks.m

## Output file
rank_data_updated.csv, this table contains 3 columns: First column is gene name, while the other 2 columns are the extracted features-  see details below. In the output file there are 6574 genes of S. cerevisiae.

## Features
Rank of gene- number of interactions that each gene creates with other genes. 
The rank of the gene is an important feature because genes that create a larger number of interactions are probably evolutionarily older, and also more conserved because if they 
are related to more genes, then most likely they tend to develop fewer mutations.
Another feature is Weighted rank_i which takes into account not only the number of interactions in which each gene is involved, but also the weight of each one of its edges and creates a kind of average measure of the strength of interaction in which the gene is involved.
 
