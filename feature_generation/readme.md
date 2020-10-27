In order to predict the best candidates for conjugate genes, we built a predictor model. However, for it to have meaningful predictions, we had to supply it with meaningful data.

In order to do so, we generated many descriptive features for the genes within the host genome, each with a potential connection to the gene's stability or expression level.
For this, we used a variety of data sources (academic publications, empiric measurements, etc.).

The models, data and code for generation of these features are presented within this directory.

Many features used the file "orf_coding_all.fasta.gz". It contains the names and sequences of all the ORFs in yeast (from SGD website). This file appears in this directory. 
