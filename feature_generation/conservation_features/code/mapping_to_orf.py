#!/usr/bin/env python
# coding: utf-8

import glob
import os
import pandas as pd


INPUT_DIR = r'../inputs/convert gene name output'
OUTPUT_DIR = r'../inputs'

# Initializing a new DataFrame
mapping_df = pd.DataFrame(columns=['ORF'])
mapping_df.index.name = 'Gene'

files = glob.glob(os.path.join(INPUT_DIR, "AJ*.txt"))

for file in files:
    file_name = os.path.basename(file)
    gene_name = os.path.splitext(file_name)[0]
    try:
        # Reading the output txt files
        with open(file, 'r') as f:
            curr_orf = f.read()
            row = pd.Series(curr_orf, index=['ORF'], name=gene_name)
            mapping_df = mapping_df.append(row)

    except Exception:
        print(gene_name)

mapping_df.to_csv(os.path.join(OUTPUT_DIR, "mapper of genes id.csv"))

# Merging files in the file 'convert_genes_names_in_tables.py'


