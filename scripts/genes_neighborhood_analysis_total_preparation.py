import pandas as pd
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool, cpu_count

# Load the necessary dataframes
pa_gene_allt_copy2 = pd.read_pickle('/home/omidard/neighbors/pa_gene_allt_copy2.pkl')
locustags_genes_mapping = pd.read_pickle('/home/omidard/neighbors/locustags_genes_mapping.pkl')
pangenome_df = pd.read_pickle('/home/omidard/neighbors/pangenome_df.pkl')

# Step 1: Prepare gene-to-category mapping for fast lookup
gene_to_category = pangenome_df.set_index('Gene')['Category'].to_dict()

# Define a function to process each row in parallel
def process_row(row_tuple):
    index, row = row_tuple  # Unpack the index and row
    genome_name, reaction, locus_tags = row['index'], row['reaction'], row['locustags']

    # Replace .json with .gbk in genome names
    genome_id = genome_name.replace('.json', '.gbk')

    # Handle locus_tags properly based on its type
    if isinstance(locus_tags, float) and pd.isna(locus_tags):
        return pd.DataFrame()  # Skip if NaN

    if isinstance(locus_tags, (list, np.ndarray, pd.Series)):
        if len(locus_tags) == 0:  # Check for empty lists/arrays
            return pd.DataFrame()
    else:
        return pd.DataFrame()  # Skip if it's not iterable

    # Prepare the expanded DataFrame for each locus tag
    rows = []
    for locus_tag in locus_tags:
        gene_name = locustags_genes_mapping.get(locus_tag, None)
        category = gene_to_category.get(gene_name, None)
        rows.append([locus_tag, reaction, genome_id, gene_name, category])

    return pd.DataFrame(rows, columns=['locus_tag', 'reaction', 'genome_id', 'gene_name', 'category'])

# Step 2: Expand `pa_gene_allt_copy2` using melt
def expand_pa_gene_allt_copy2(pa_gene_allt_copy2):
    expanded_df = pa_gene_allt_copy2.reset_index().melt(id_vars='index', var_name='reaction', value_name='locustags')
    return expanded_df

# Step 3: Apply multiprocessing with progress bar
def process_all_reactions(expanded_df):
    with Pool(cpu_count()) as pool:
        results = list(tqdm(pool.imap_unordered(process_row, expanded_df.iterrows(), chunksize=100),
                            total=len(expanded_df), desc="Processing reactions"))
    # Concatenate the results into a final DataFrame
    final_df = pd.concat(results, ignore_index=True)
    return final_df

# Step 4: Main execution
if __name__ == '__main__':
    # Expand the DataFrame first
    expanded_df = expand_pa_gene_allt_copy2(pa_gene_allt_copy2)

    # Process all rows in parallel
    final_df = process_all_reactions(expanded_df)

    # Save the final DataFrame
    final_df.to_csv('/home/omidard/neighbors/final_gene_neighborhood.csv', index=False)
    print('Final dataframe saved as /home/omidard/neighbors/final_gene_neighborhood.csv')

