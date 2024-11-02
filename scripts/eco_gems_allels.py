import os
import pandas as pd
import json
from multiprocessing import Pool, cpu_count, Manager
from functools import partial
from cobra.io import load_json_model
from tqdm import tqdm

# Step 1: Load and modify the header_to_allele dictionary
def load_and_modify_dict(file_path):
    header_to_allele = pd.read_pickle(file_path)
    modified_dict = {}
    for key, value in header_to_allele.items():
        new_value = value.split('A')[0]
        modified_dict[key] = new_value
    return modified_dict

file_path = "/home/omidard/header_to_allele.pickle"

# Load and modify the dictionary once
locustags_genes_mapping = load_and_modify_dict(file_path)

# Step 2: Function to process a single GEM file using COBRApy
def process_gem_file(gem_file, model_folder, locustags_genes_mapping):
    model_path = os.path.join(model_folder, gem_file)
    model = load_json_model(model_path)
    reaction_gprs = {}
    for reaction in model.reactions:
        genes = list(reaction.genes)
        if genes:
            gprs = [locustags_genes_mapping.get(gene.id, gene.id) for gene in genes]
            reaction_gprs[reaction.id] = gprs
    return gem_file, reaction_gprs

# Helper function to update the progress bar (called in the main process)
def update_progress_bar(result, progress_bar, progress_counter, lock):
    with lock:
        progress_counter.value += 1
        progress_bar.update(1)
    return result

model_folder = '/home/omidard/gapfilled3'
gem_files_df = pd.read_csv('complete_gems.csv', index_col=0)
gem_files = gem_files_df['gems'].tolist()

num_cores = min(cpu_count(), 64)

# Processing all GEM files
all_gem_files = gem_files

# Step 3: Use multiprocessing to process GEM files with a progress bar
with Manager() as manager:
    shared_locustags_genes_mapping = manager.dict(locustags_genes_mapping)
    progress_counter = manager.Value('i', 0)
    lock = manager.Lock()
    
    results = []  # Define the results list before the loop

    with tqdm(total=len(all_gem_files), desc="Processing GEM files") as progress_bar:
        partial_process_gem = partial(process_gem_file, model_folder=model_folder, locustags_genes_mapping=shared_locustags_genes_mapping)
        
        with Pool(num_cores) as pool:
            for result in pool.imap_unordered(partial_process_gem, all_gem_files):
                update_progress_bar(result, progress_bar, progress_counter, lock)
                results.append(result)

    # Collect results from multiprocessing
    all_reaction_ids = set()
    gem_reactions = {}
    for gem_file, reaction_gprs in results:
        all_reaction_ids.update(reaction_gprs.keys())
        gem_reactions[gem_file] = reaction_gprs

# Step 4: Create the DataFrame and fill with GPRs
reaction_ids = sorted(all_reaction_ids)
df = pd.DataFrame(index=all_gem_files, columns=reaction_ids, dtype=object)  # Ensure dtype=object to store lists

for gem_file, reaction_gprs in gem_reactions.items():
    for reaction_id, gprs in reaction_gprs.items():
        df.at[gem_file, reaction_id] = str(gprs)  # Convert the list to a string

# Step 5: Save the DataFrame to CSV
df.to_csv('/home/omidard/ecoli_gprs.csv')
