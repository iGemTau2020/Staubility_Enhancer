#!/usr/bin/env python
# coding: utf-8


import argparse
import os
import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio import SeqIO
import glob
import pickle

INPUT_DIR = r'path_in_privet_server'
OUTPUT_DIR = r'../inputs'
families_hierarchy = ['at4892', 'at4891', 'at4890'] # 4892 - Orthologous genes from the same order ( they are close to S. cerevisiae but not too much), 4891 - Orthologous genes from the same class
# 4890 Orthologous genes from the same pylum (The most evolutionarily distant amoung the other orthological families (excluding genes from the same kingdom(4751) and super-kingdom(2759)))



def get_df_per_gene(msa):
    cols_num = msa.get_alignment_length()
    chars_cols = list(range(cols_num))
    records_num = len(msa)
    genes_df = pd.DataFrame(columns=['creature', 'gene'] + chars_cols, index=range(records_num))
    genes_df.index.name = 'Gene'
    i = 0
    for record in msa:
        record_seq = str(record.seq) 
        record_id = record.id
        creature_id = record_id.split('_')[0]
        gene_id = record_id.split(':')[1]
        chars = list(record_seq)
        data = [creature_id, gene_id] + chars
        #index = ['creature', 'gene'] + chars_cols
        #seq_series = pd.Series(data = data, index = index)
        genes_df.iloc[i] = data
        i += 1
    return genes_df

def indels_percentage(genes_df):
### get percentage of a specific gene in S. cerevisiae
    number_of_positions = len(genes_df.columns) - 2 # number of columns minus the columns: gene_id and creature_id
    numbers_of_rows = len(genes_df.index) 
    indels_percentage_array = genes_df.apply(lambda x: x.value_counts(normalize=True, sort=False)).loc['-']
    return indels_percentage_array.mean(), indels_percentage_array.median()

def entropy_calc(data_for_entropy_calc, remove_indels=False):
    if remove_indels:
        data_for_entropy_calc = data_for_entropy_calc.copy()
        data_for_entropy_calc = data_for_entropy_calc.replace('-', np.NaN)
    p_df = data_for_entropy_calc.apply(lambda x: x.value_counts(normalize=True, sort=False, dropna=remove_indels))
    entropy_series = p_df.apply(lambda x: -(x * np.log2(x)).sum())

    return entropy_series.values

def prepere_df_for_entropy_func(genes_df):
    # Calculation entropy for all the MSA including indels
    number_of_positions = len(genes_df.columns) - 2 
    all_positions_list = list(range(number_of_positions))
    data_with_indels_for_entropy_calc = genes_df[all_positions_list]
    # Calculation entropy for all the MSA without indels
    return data_with_indels_for_entropy_calc

def feature_extraction(genes_df, columns, remove_indels):
    # Extaction of gene id
    gene_id = genes_df['gene'][genes_df['creature'] == '1294385'].values[0]   
    # Extaction of feaures
    # Indels percentage calculation
    mean_indels_percentage, median_indels_percentage = indels_percentage(genes_df)
    # Entropy calculation
    data_for_entropy_calc = prepere_df_for_entropy_func(genes_df)
    entropy_array = entropy_calc(data_for_entropy_calc, remove_indels)
    mean_entropy = np.mean(entropy_array)
    median_entropy = np.median(entropy_array)
    # Calculation of entropy on the first 100 windows and last 100 windows in window length of 30 amino- acids
    beginning_windows = np.zeros(100)
    end_windows = np.zeros (100)
    entropy_length = len(entropy_array)
    for i in range(100):
        curr_beginning_window = entropy_array[i:i+30]
        curr_end_window = entropy_array[entropy_length-i-30:entropy_length-i]
        beginning_windows[i] = np.mean(curr_beginning_window)
        end_windows[i] = np.mean(curr_end_window)
    # Find the first windows with maximum and minimum value of entropy in the first 100 windows
    first_window_in_the_beginning_with_max_entropy = np.argmax(beginning_windows) + 1
    first_window_in_the_beginning_with_min_entropy = np.argmin(beginning_windows) + 1
    # Find the first windows with maximum and minimum value of entropy in the last 100 windows
    first_window_in_the_end_with_max_entropy = np.argmax(end_windows) + 1
    first_window_in_the_end_with_min_entropy = np.argmin(end_windows) + 1

    # Collecting all the relevant values
    entropy_values = [mean_indels_percentage, median_indels_percentage, mean_entropy, median_entropy]
    windows_with_max_and_min_entropy = [first_window_in_the_beginning_with_max_entropy, first_window_in_the_beginning_with_min_entropy,
                                        first_window_in_the_end_with_max_entropy, first_window_in_the_end_with_min_entropy]
    all_values = entropy_values + list(beginning_windows) + list(end_windows) + windows_with_max_and_min_entropy
    # Saving all the variables that I calculated in a series which is row in statistics_df (the final df)
    row_in_genes_statistics = pd.Series(all_values, index=columns, name=gene_id) # A row in statistics_df
    return row_in_genes_statistics

def run_for_all_files():
    ## Main code and calling to functions
    # Columns names procedure for the final df
    base_columns = ['mean_indels_percentage', 'median_indels_percentage','mean_entropy','median_entropy'] +\
                   ['beginning_window_' + str(x + 1) for x in range(100)] + ['end_window_' + str(x + 1) for x in range(100)] + ['first_window_in_the_beginning_with_max_entropy',
                    'first_window_in_the_beginning_with_min_entropy', 'first_window_in_the_end_with_max_entropy', 'first_window_in_the_end_with_min_entropy']
    columns_with_indels = [x + '_with_indels' for x in base_columns]
    columns_no_indels = [x + '_no_indels' for x in base_columns]
    all_columns = columns_with_indels + columns_no_indels
    # Initializing a df for all the features ( entropy values for each yeast's genes)

    statistics_df = pd.DataFrame(columns=all_columns, index=list(range(7000)))
    statistics_df.index.name = 'Gene'
    i = 0 # number of rows have been populated in statistic_df
    for family in families_hierarchy:
        file_name_pattern = "og_*{}.faa.muscle.clw".format(family)
        files = glob.glob(os.path.join(INPUT_DIR, file_name_pattern))
        for file in files: # Running only on file from the specific family
            print(file)
            entropy_values = []
        #    file_name = os.path.basename(file)
            # Reading each MSA file
            msa = AlignIO.read(file, 'clustal')
            # Extracting functions on each file ( calculate the features first with indels and then without indels)
            genes_df = get_df_per_gene(msa)
            gene_id = genes_df['gene'][genes_df['creature'] == '1294385'].values[0]
            if gene_id in statistics_df.index.values:
                continue
            row_in_genes_statistics_with_indels = feature_extraction(genes_df, columns_with_indels, False)
            row_in_genes_statistics_no_indels = feature_extraction(genes_df, columns_no_indels, True)
            combined_series = pd.concat([row_in_genes_statistics_with_indels, row_in_genes_statistics_no_indels])
            statistics_df.iloc[i] = combined_series
            statistics_df.rename(index={i: gene_id}, inplace=True) # change the index of the static df becaues it is file name defaultly
            i += 1
        statistics_df.to_csv(os.path.join(OUTPUT_DIR, "{}.csv".format(family)))

    rows_to_remove = list(range(i , 7000))
    statistics_df.drop(rows_to_remove, inplace  = True)   #Droping all the empty rows (from i to end)
    statistics_df.to_csv(os.path.join(OUTPUT_DIR, "statistics.csv"))

def run_single_file(file, ouput):
    base_columns = ['mean_indels_percentage', 'median_indels_percentage', 'mean_entropy', 'median_entropy'] + \
                   ['beginning_window_' + str(x + 1) for x in range(100)] + ['end_window_' + str(x + 1) for x in
                                                                             range(100)] + [
                       'first_window_in_the_beginning_with_max_entropy',
                       'first_window_in_the_beginning_with_min_entropy', 'first_window_in_the_end_with_max_entropy',
                       'first_window_in_the_end_with_min_entropy']
    columns_with_indels = [x + '_with_indels' for x in base_columns]
    columns_no_indels = [x + '_no_indels' for x in base_columns]
    all_columns = columns_with_indels + columns_no_indels
    # Initializing a df for all the features ( entropy values for each yeast's genes)


    msa = AlignIO.read(file, 'clustal')
    genes_df = get_df_per_gene(msa)
    gene_id = genes_df['gene'][genes_df['creature'] == '1294385'].values[0]
    row_in_genes_statistics_with_indels = feature_extraction(genes_df, columns_with_indels, False)
    row_in_genes_statistics_no_indels = feature_extraction(genes_df, columns_no_indels, True)
    combined_series = pd.concat([row_in_genes_statistics_with_indels, row_in_genes_statistics_no_indels])
    statistics_df = pd.DataFrame(columns=all_columns, index=[gene_id])
    statistics_df.index.name = 'Gene'
    statistics_df.iloc[0] = combined_series
    statistics_df['file_name'] = os.path.basename(file)
    statistics_df.to_csv(ouput)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process features from msa files')
    parser.add_argument('--run-all-files', dest='run_all_files', action='store_true',
                        help='Run for all files in input directory')
    parser.add_argument('-file', dest='file', type=str,
                        help='path to single file to process')
    parser.add_argument('-output', dest='output', type=str,
                        help='path to output file')
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    if args.run_all_files:
        run_for_all_files()
    else:
        run_single_file(args.file, args.output)

if __name__ == '__main__':
    main()
