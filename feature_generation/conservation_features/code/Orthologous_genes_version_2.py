#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import numpy as np
import glob

INPUT_DIR = r'path_in_privet_server'
OUTPUT_DIR = r'../inputs'

families_hierarchy = ['at4893', 'at4891', 'at4890', 'at4751', 'at2759']

def geometric_mean_of_conservation(AA_cons_tbl, nt_cons_tbl):
    sequence_length = len(AA_cons_tbl)
    sum_AA = np.sum(np.log(AA_cons_tbl))
    curr_AA_geometric_cons = np.exp(sum_AA / sequence_length) # calculate Amino Acids entropy value for all cretures group
    curr_nt_geometric_cons = np.exp(np.log(nt_cons_tbl).values.sum() / (3 * sequence_length)) # calculate nucleotide entropy value for all cretures group
    # return calculate values
    return curr_AA_geometric_cons, curr_nt_geometric_cons
    


# In[4]:


def median_of_conservation (AA_cons_tbl, nt_cons_tbl):
    curr_AA_median_cons = AA_cons_tbl.median()
    curr_nt_median_cons = np.median(nt_cons_tbl.values)
    return curr_AA_median_cons, curr_nt_median_cons


# In[5]:


def find_closest_group (curr_gene):
    for family in families_hierarchy:
        filter_cols = curr_gene.apply(lambda x: x.to_string().find(family) > 0)
        if len(filter_cols.loc[filter_cols]) == 0:
            continue
        closest_col_name = filter_cols.loc[filter_cols].index[0]
        col_index = list(curr_gene.columns).index(closest_col_name) + 1
        nt_closest_col_1 = list(curr_gene.columns)[col_index + 1]
        nt_closest_col_2 = list(curr_gene.columns)[col_index + 2]
        nt_closest_col_3 = list(curr_gene.columns)[col_index + 3]
        AA_closest_col = list(curr_gene.columns)[col_index + 4]
        nt_cons_closest_tbl = curr_gene[[nt_closest_col_1, nt_closest_col_2, nt_closest_col_3]]
        AA_cons_closest = curr_gene[AA_closest_col]

    return AA_cons_closest, nt_cons_closest_tbl


# In[6]:


def find_beginning_and_end_locations (AA_cons, nt_cons_tbl):
    # Taking 50 amino acids/codons from the begining and from the end of each gene
    AA_beg = AA_cons[0:50]
    AA_end = AA_cons[-50:]
    nt_beg = nt_cons_tbl[0:50]
    nt_end = nt_cons_tbl[-50:]
    return AA_beg, AA_end, nt_beg, nt_end


files = glob.glob(os.path.join(INPUT_DIR, "gene_*.csv"))
columns_name = ['AA index', 'n1', 'n2', 'n3', 'AA',
                'first OG', 'location_1','n1_1','n2_1','n3_1','AA_1', 'unknown_1', 
                'second OG', 'location_2','n1_2','n2_2','n3_2','AA_2', 'unknown_2',
                'third OG', 'location_3','n1_3','n2_3','n3_3','AA_3', 'unknown_3',
                'fourth OG', 'location_4','n1_4','n2_4','n3_4','AA_4', 'unknown_4',
                'fifth OG','location_5','n1_5','n2_5','n3_5','AA_5', 'unknown_5'] # names to the columns

columns = ['AA_cons_geometric_all', 'nt_cons_geometric_all', 'AA_cons_median_all',
           'nt_cons_median_all', 'AA_cons_geometric_closest', 'nt_cons_geometric_closest',
          'AA_cons_median_closest', 'nt_cons_median_closest', 'AA_cons_geometric_all_beg', 'nt_cons_geometric_all_beg',
          'AA_cons_median_all_beg', 'nt_cons_median_all_beg', 'AA_cons_geometric_all_end', 'nt_cons_geometric_all_end',
          'AA_cons_median_all_end', 'nt_cons_median_all_end','AA_cons_geometric_closest_beg', 'nt_cons_geometric_closest_beg',
          'AA_cons_median_closest_beg', 'nt_cons_median_closest_beg', 'AA_cons_geometric_closest_end', 'nt_cons_geometric_closest_end',
          'AA_cons_median_closest_end', 'nt_cons_median_closest_end']

results_df = pd.DataFrame(columns=columns)
results_df.index.name = 'Gene'

for file in files:
    file_name = os.path.basename(file)
    file_name_no_extension = os.path.splitext(file_name)[0]
    gene_name = file_name_no_extension.split('_')[1]
    try:
        curr_gene = pd.read_csv(file, header = None, names = columns_name, index_col = 0)
        AA_cons_all = curr_gene['AA'] # Amino acid conservation for all creatures level
        nt_cons_all_tbl = curr_gene[['n1','n2','n3']] # nucleotides conservation for all creatures level
        current_values = []
    
   	# Conservation of all creatures all gene length
        current_values.extend(geometric_mean_of_conservation(AA_cons_all, nt_cons_all_tbl))
        current_values.extend(median_of_conservation(AA_cons_all, nt_cons_all_tbl))
    
   	# Conservation of closest creatures all gene length
        AA_cons_closest, nt_cons_closest_tbl = find_closest_group(curr_gene)
        current_values.extend(geometric_mean_of_conservation(AA_cons_closest, nt_cons_closest_tbl))
        current_values.extend(median_of_conservation(AA_cons_closest, nt_cons_closest_tbl))
    
   	# Conservation of all creatures just on the begining or end of the gene 
        AA_beg_all, AA_end_all, nt_beg_all, nt_end_all = find_beginning_and_end_locations (AA_cons_all, nt_cons_all_tbl)
        current_values.extend(geometric_mean_of_conservation(AA_beg_all, nt_beg_all))
        current_values.extend(median_of_conservation(AA_beg_all, nt_beg_all))
        current_values.extend(geometric_mean_of_conservation(AA_end_all, nt_end_all))
        current_values.extend(median_of_conservation(AA_end_all, nt_end_all))
    
   	# Conservation of closest creatures just on the begining or end of the gene 
        AA_beg_closest, AA_end_closest, nt_beg_closest, nt_end_closest = find_beginning_and_end_locations (AA_cons_closest, nt_cons_closest_tbl)         
        current_values.extend(geometric_mean_of_conservation(AA_beg_closest, nt_beg_closest))
        current_values.extend(median_of_conservation(AA_beg_closest, nt_beg_closest))
        current_values.extend(geometric_mean_of_conservation(AA_end_closest, nt_end_closest))
        current_values.extend(median_of_conservation(AA_end_closest, nt_end_closest))                      

        row = pd.Series(current_values, index=columns, name=gene_name)
        results_df = results_df.append(row)
    except:
        print(gene_name)
results_df.to_csv(os.path.join(OUTPUT_DIR, "conservation_features.csv"))
