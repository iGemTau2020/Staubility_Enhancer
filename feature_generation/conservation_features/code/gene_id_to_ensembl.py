# In this script I convert the gene name from Ortho DB internal IDs to Ensembl IDs
import pandas as pd
import pickle
import os

INPUT_FILE = r"../inputs/filtered.csv"
INPUT_DIR = r'../inputs'
OUTPUT_DIR = r'../inputs'

st = pd.read_csv(INPUT_FILE)
st.set_index('Gene', inplace=True)

# Loading the dicitonaries of mapping between OrthoDB ids and Ensemble ids
with open(os.path.join(INPUT_DIR, 'ensembelsId.pkl'), "rb") as f:
    geneEnsemblIds, ensemblSequenceIdsForGenes = pickle.load(f)
    ensemblSequenceIdsForGenes = ensemblSequenceIdsForGenes[1294385] # Get only genes of S. cerevisiae

odb_gene_names = st.index  # gene names as OrthoDB internal IDs.

def odb_ids_to_ensemble_ids(x):
    tax_id = '1294385_1:' + x   # 1294385_1 is taxonomy id for S. cerevisiae
    xrefid = geneEnsemblIds[tax_id][0]
    ensembl_id = ensemblSequenceIdsForGenes[xrefid]
    return ensembl_id

ensembl_ids = odb_gene_names.map(odb_ids_to_ensemble_ids) # gene names as Ensembl IDs.
st.index = ensembl_ids
st.to_csv(os.path.join(OUTPUT_DIR, "Ensembl_statistics.csv"))
