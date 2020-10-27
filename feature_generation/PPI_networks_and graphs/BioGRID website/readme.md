This script aims to create PPI features from BIOGRID's data.

**Inputs:**

1.	BIOGRID-ALL-3.5.184.tab3.zip (Can easily be downloaded from the BIOGRID website) 

2.	biogrid-yeastyeast.csv – an edited version of BIOGRID-ALL-3.5.184.tab3.zip.

**Codes:**

1.	'firstrun.py' and 'secondrun.py'- these scripts create a CSV file that contains only yeast-yeast interactions from the original file that contains interactions between major model organism species (Input=input1 , output=input2).
*Please run firstrun before seconrun.

2.	'features.m'- this script calculates some features from this data (Input=input2,Output=output1).

**Output: 'biogrid final.csv' (columns=features, rows= yeast's genes). Please visit our wiki,model page for further explanations.**  

*The features:*

a)	Rank- how many interactions does each protein have with other proteins. We describe in the code how do we count the interactions and which interactions we consider as duplications. 

b)	Interaction type- How many of the protein’s interactions are 'genetic' and how many are 'physical', as seen in the column Experimental System Type.

c)	iExperimental System Name – How many of the protein’s interactions have each one of the possible values in the column Experimental System Name in the original data ('Affinity Capture-MS','Affinity Capture-Western', 'Dosage Rescue', 'Synthetic Growth Defect', 'Synthetic Lethality', 'Synthetic Rescue', 'Two-hybrid', 'Reconstituted Complex', 'Biochemical Activity', 'Co-crystal Structure', 'Far Western', 'Protein-peptide', 'FRET', 'Co-localization', 'Affinity Capture-RNA', 'Protein-RNA', 'PCA', 'Co-purification', 'Co-fractionation', 'Dosage Lethality', 'Phenotypic Enhancement', 'Phenotypic Suppression', 'Dosage Growth Defect', 'Negative Genetic', 'Positive Genetic', 'Synthetic Haploinsufficiency', 'Affinity Capture-Luminescence', 'Proximity Label-MS'). 

d)	Interaction Throughput- How many of the protein’s interactions have each one of the possible values in the column Interaction Throughput (High, Low, High and Low).

**NOTE: The features calculated in the script were copied to the output file manually.**
