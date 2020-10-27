The scripts that merge the data to a single table and then filter the features are presented in this directory.

The script "Merge_Tables_and_Filter_Genes_And_Features-NO-PARALOGS.ipynb": 
* merges all the feature tables generated for this project, using the genes' ORF name. 
* filters genes with too many missing values, and features with too many missing values. 
* Another important filtering is of paralog genes. We remove them from analysis because they will not contribute to the target gene's stability when conjugated. The list "Gene_YeastParalogs_only_ORF.csv" contains the ORF names of paralog genes, extracted from the yeastmine website (http://yeastmine.yeastgenome.org:8080/yeastmine/template.do?name=Gene_YeastParalogs&scope=global).  

The output is the full features' table, normalized. 

This output table then goes to the "feature selection" module for further filtering. 

The script "Filter_flourocensce_tables.ipynb" sorts the fluorescence tables of GFP and RFP in the same order as the features' table, allowing easier filtration in the next scripts.
