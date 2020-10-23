import os
import pandas as pd

RELEVANT_DIR = r'../inputs'

mapper = pd.read_csv(os.path.join(RELEVANT_DIR, "mapper of genes id.csv")) # Reading the mapper
entropy_genes_df = pd.read_csv(os.path.join(RELEVANT_DIR, 'Ensembl_statistics.csv'))  # Reading my table
conservation_genes_df = pd.read_csv(os.path.join(RELEVANT_DIR, 'conservation_features.csv'))  # Reading my table
conservation_genes_df['Gene'] = [g.split('.')[0] for g in conservation_genes_df.Gene.values]  # Removing the dot and the number after the dot from the protein id in conservation_genes_df

# Merging tables
entropy_genes_df_converted = entropy_genes_df.merge(mapper, how='left', on='Gene')
conservation_genes_df_converted = conservation_genes_df.merge(mapper, how='left', on='Gene')

# Dropping irrelevant columns
entropy_genes_df_converted.set_index('ORF', inplace = True)
conservation_genes_df_converted.set_index('ORF', inplace = True)
entropy_genes_df_converted.drop(['Gene', 'file_name'], axis='columns', inplace=True)
conservation_genes_df_converted.drop(['Gene'], axis='columns', inplace=True)


entropy_genes_df_converted.to_csv(os.path.join(RELEVANT_DIR, " ../outputs/entropy features gene names converted.csv"))
conservation_genes_df_converted.to_csv(os.path.join(RELEVANT_DIR, "../outputs/conservation featurs gene names converted.csv"))