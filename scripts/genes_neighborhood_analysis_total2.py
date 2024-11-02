import pandas as pd
import os
import multiprocessing as mp
from Bio import SeqIO
from tqdm import tqdm

# Define file paths
final_gene_neighborhood_path = '/home/omidard/neighbors/final_gene_neighborhood.csv'
genomes_dir = '/home/omidard/prokka_genomes/'

# Load filtered DataFrame (final_gene_neighborhood) and drop rows with missing values
filtered_df = pd.read_csv(final_gene_neighborhood_path)

# Drop rows with missing `locus_tag`, `gene_name`, or `category`
filtered_df.dropna(subset=['locus_tag', 'gene_name', 'category'], inplace=True)

# Define a function to find neighbors for a given locus tag in a specific genome file
def find_neighbors(locus_tag, genome_name):
    # Ensure genome_name has .gbk extension
    if not genome_name.endswith('.gbk'):
        genome_name = f"{genome_name}.gbk"
    
    genome_file = os.path.join(genomes_dir, genome_name)
    
    if not os.path.isfile(genome_file):
        print(f"Genome file not found: {genome_file}")
        return None

    neighbors = {
        'Downstream_Gene_1': None,
        'Downstream_Gene_2': None,
        'Target_Gene': None,
        'Upstream_Gene_1': None,
        'Upstream_Gene_2': None
    }

    try:
        # Read the GenBank file
        for record in SeqIO.parse(genome_file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    gene_name = feature.qualifiers.get("locus_tag", ["Unknown"])[0]
                    protein_name = feature.qualifiers.get("product", ["Unknown"])[0]

                    # Check if this is the target gene
                    if gene_name == locus_tag:
                        neighbors['Target_Gene'] = protein_name

                        # Find neighbors in the same record
                        index = record.features.index(feature)
                        if index > 0:
                            downstream_1 = record.features[index - 1]
                            neighbors['Downstream_Gene_1'] = downstream_1.qualifiers.get("product", ["Unknown"])[0] if downstream_1.type == "CDS" else None
                        if index > 1:
                            downstream_2 = record.features[index - 2]
                            neighbors['Downstream_Gene_2'] = downstream_2.qualifiers.get("product", ["Unknown"])[0] if downstream_2.type == "CDS" else None
                        if index < len(record.features) - 1:
                            upstream_1 = record.features[index + 1]
                            neighbors['Upstream_Gene_1'] = upstream_1.qualifiers.get("product", ["Unknown"])[0] if upstream_1.type == "CDS" else None
                        if index < len(record.features) - 2:
                            upstream_2 = record.features[index + 2]
                            neighbors['Upstream_Gene_2'] = upstream_2.qualifiers.get("product", ["Unknown"])[0] if upstream_2.type == "CDS" else None
                        break
    except Exception as e:
        print(f"Error processing genome file {genome_file}: {e}")

    return neighbors

# Define a function to process each row of the filtered dataframe
def process_row(row):
    locus_tag = row['locus_tag']
    genome_name = row['genome_id']
    neighbors = find_neighbors(locus_tag, genome_name)
    
    # Add the reaction column to the output dictionary
    if neighbors is not None:
        neighbors['Reaction'] = row['reaction']
    
    return neighbors

# Define a wrapper function to allow the row to be processed as a tuple
def process_row_wrapper(row_tuple):
    return process_row(row_tuple[1])

# Apply multiprocessing to process each locus tag across all reactions
def process_filtered_df(df):
    with mp.Pool(processes=96) as pool:
        results = list(tqdm(pool.imap(process_row_wrapper, df.iterrows()), total=len(df), desc="Processing Locus Tags"))
    
    # Create DataFrame for results
    results_df = pd.DataFrame(results)
    
    # Merge results with the original filtered DataFrame
    final_df = pd.concat([df.reset_index(drop=True), results_df], axis=1)
    return final_df

# Run the processing
final_df = process_filtered_df(filtered_df)

# Save the final DataFrame with all reactions
final_output_path = '/home/omidard/neighbors/all_reactions_gene_neighborhood.csv'
final_df.to_csv(final_output_path, index=False)
print(f'{final_output_path} saved')
