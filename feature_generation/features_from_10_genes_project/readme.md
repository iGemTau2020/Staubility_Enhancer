As part of the POC experiment, we calculated many initial features. These features were calculated using the scripts in this directory. 
They were generated for all 6,000 genes in the S. cerevisiae, except for several features like dn / ds and CAI, which had many missing values. 

We recalculate the CAI feature for all genes, and instead of the dn / ds feature, we generate several evolutionary conservation features, 
as we know this type of feature is important for our predictor.


These features include:
 1. features related to length and location of gene, such as distance from end of chromosome. These are generated using length_and_loc.ipynb.
 2. fluorescence features, related to the empirical data gathered by Prof. Schuldiner. These are generated using flourescence.ipynb.
 3. bio features, including many features derived by the biopython library. These include GC content, aromaticity, alpha sheets, and many others. These are generated using bio_features.ipynb.
 
 these features are combined to a single table using table_merger.ipynb, and enumerated using enumeration_full_table.ipynb.
