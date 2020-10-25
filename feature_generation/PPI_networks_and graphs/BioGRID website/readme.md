Source: BioGRID website.

Pre-processing:
Input: BIOGRID-ALL-3.5.184.tab3.zip (Can be easily downloaded from the website)
Code: In the preprocessing we found out that there are columns in the table named  Organism ID Interactor A/B and we found in the website that the ID of yeast is 559292.c (https://wiki.thebiogrid.org/doku.php/biogridrest#organisms_-_fetching_supported_organisms_list)
We used firstrun and secondrun python codes in order to create a csv file that contains only yesast-yeast and yeast-other creatures' interactions.
The number of interactions in the file is identical to the one detailed in BIOGRID website  (https://wiki.thebiogrid.org/doku.php/build_3.5.184- in row Saccharomyces cerevisiae (S288c) , Raw interactions , combined) .
Afterwards, we separated into two different csv files the interactions of yeast-yeast and interactions of yeast with other creatures by using the ID column we mentioned before.
The next step was to choose the relevant columns from the tables we created. The columns' names are detailed here- https://wiki.thebiogrid.org/doku.php/biogrid_tab_version_3.0. 
This page was updated after we downloaded the data , columns 17 and 19 were added to it. 
The following columns were found relevant to our future processing: Systematic name Interactor A, Systematic name Interactor B, Official symbol Interactor A, Official symbol Interactor B, Experimental System Name, Experimental System Type, Organism ID Interactor A, Organism ID Interactor B, Interaction Throughput.
Output: We created a new csv file named 'biogrid-yeastyeast.csv' that contains the data in these columns and appears here.

Features' creation:
Input: 'biogrid-yeastyeast.csv'
Code: We used MATLAB script – 'features.m' in order to calculate some features from the data. The features we decided to compute are:
i.	Rank- how many interactions does each protein have with other proteins. We describe in the code how do we count the interactions and which interactions we consider as duplications. We noticed that sometimes we would have the same interaction multiple times and we decided to count it once. Moreover, we found out that sometimes an interaction is counted twice- once A and B are the interactors and once B and A, we removed these duplicates as well.
ii.	Interaction type- How many of the protein’s interactions are 'genetic' and how many are 'physical', as seen in the column Experimental System Type.
iii.	Experimental System Name – How many of the protein’s interactions have each one of the possible values in the column Experimental System Name in the original data ('Affinity Capture-MS','Affinity Capture-Western', 'Dosage Rescue'	, 'Synthetic Growth Defect', 'Synthetic Lethality', 'Synthetic Rescue', 'Two-hybrid', 'Reconstituted Complex', 'Biochemical Activity', 'Co-crystal Structure', 'Far Western', 'Protein-peptide', 'FRET', 'Co-localization', 'Affinity Capture-RNA', 'Protein-RNA', 'PCA', 'Co-purification', 'Co-fractionation', 'Dosage Lethality', 'Phenotypic Enhancement', 'Phenotypic Suppression', 'Dosage Growth Defect', 'Negative Genetic', 'Positive Genetic', 'Synthetic Haploinsufficiency', 'Affinity Capture-Luminescence', 'Proximity Label-MS)'. 
iv.	Interaction Throughput- How many of the protein’s interactions have each one of the possible values in the column Interaction Throughput (High, Low, High and Low).
The features I calculated in the script were copied to a new csv file manually and it is the final product of this project ('biogrid final').
Output:'biogrid final.xlsx'
