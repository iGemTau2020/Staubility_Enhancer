# In this script I  took the table extracted in the script "Entropy_version_2.py" and filtered only the orthologous genes from the same class (contain 'at4891' in their  OrthoDB id)
# then, I added only new genes from the same pylum ( contain 'at4890')
import pandas as pd

families_hierarchy = ['at4891', 'at4890'] # 4891 - Orthologous genes from the same class (they are close to S. cerevisiae but not too much),
# 4890 Orthologous genes from the same pylum (The most evolutionarily distant amoung the other orthological families (excluding genes from the same kingdom(4751) and super-kingdom(2759)))

INPUT_FILE = r"../inputs/all_statistics.csv"
st = pd.read_csv(INPUT_FILE)
st.set_index('Gene', inplace=True)
output_df = None
current_position = 0
for family in families_hierarchy:
    family_df = st[st['file_name'].str.contains(family)]
    if output_df is None:
        output_df = family_df
    else:
        f_index = set(family_df.index)
        current_index = set(output_df.index)
        new_genes = f_index.difference(current_index)
        output_df = output_df.append(family_df.loc[new_genes])
output_df.to_csv('../inputs/filtered.csv')