The code merges all the features tables that were generated for this project according to the ORF name. 

It also filters genes with too many missing values, and features with too many missing values. 

Another important filtering is of paralogs genes, that we remove from analysis because they will not contribute to the target gene's stability when conjugated. The list "Gene_YeastParalogs_only_ORF.csv" contains the ORF names of paralogs genes, derived from yeastmine website (http://yeastmine.yeastgenome.org:8080/yeastmine/template.do?name=Gene_YeastParalogs&scope=global).  

The output is the final features' table, normalized. 

This output table then goes to the "feature selection" module for further filtering. 
