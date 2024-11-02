import os
import cobra
import multiprocessing

# Load the reference model (in MATLAB .mat format)
ref_model = cobra.io.load_matlab_model('/home/omidard/ref_model_dir/marlbr2.mat')

# Filter reactions in the reference model
target_reactions = [reaction for reaction in ref_model.reactions if 's0001' in [gene.id for gene in reaction.genes]]

# Function to process a target model
def process_target_model(target_model_path):
    # Load the target model (in JSON format)
    target_model = cobra.io.load_json_model(target_model_path)

    # Add the target reactions to the target model
    target_model.add_reactions(target_reactions)

    # Save the modified target model (in JSON format)
    target_model_filename = os.path.basename(target_model_path)
    cobra.io.json.save_json_model(target_model, f'/home/omidard/gapfilled/{target_model_filename}')

# List of target model file paths (in JSON format)
target_model_dir = '/home/omidard/gapfilled'
target_model_files = [os.path.join(target_model_dir, filename) for filename in os.listdir(target_model_dir) if filename.endswith('.json')]

# Use multiprocessing to parallelize the process
num_cores = 64
pool = multiprocessing.Pool(processes=num_cores)
pool.map(process_target_model, target_model_files)
pool.close()
pool.join()

