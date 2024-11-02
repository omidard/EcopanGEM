import os
import cobra
import time
from cobra.medium import minimal_medium
from cobra.io.json import save_json_model
import multiprocessing
import warnings
warnings.filterwarnings("ignore")
# Define the M9 medium function
def m9(model):
    for reaction in model.reactions:
        if 'EX_' in  reaction.id:
            reaction.lower_bound=0 
            
    model.reactions.EX_ca2_e.lower_bound=-1000
    model.reactions.EX_cl_e.lower_bound=-1000
    model.reactions.EX_co2_e.lower_bound=-1000
    model.reactions.EX_cobalt2_e.lower_bound=-1000
    model.reactions.EX_cu2_e.lower_bound=-1000
    model.reactions.EX_fe2_e.lower_bound=-1000
    model.reactions.EX_fe3_e.lower_bound=-1000
    model.reactions.EX_h_e.lower_bound=-1000
    model.reactions.EX_h2o_e.lower_bound=-1000
    model.reactions.EX_k_e.lower_bound=-1000
    model.reactions.EX_mg2_e.lower_bound=-1000
    model.reactions.EX_mn2_e.lower_bound=-1000
    model.reactions.EX_mobd_e.lower_bound=-1000
    model.reactions.EX_na1_e.lower_bound=-1000
    model.reactions.EX_tungs_e.lower_bound=-1000
    model.reactions.EX_zn2_e.lower_bound=-1000
    model.reactions.EX_ni2_e.lower_bound=-1000
    model.reactions.EX_sel_e.lower_bound=-1000
    model.reactions.EX_slnt_e.lower_bound=-1000
    model.reactions.EX_glc__D_e.lower_bound=-20
    model.reactions.EX_so4_e.lower_bound=-1000
    model.reactions.EX_nh4_e.lower_bound=-1000
    model.reactions.EX_pi_e.lower_bound=-1000
    model.reactions.EX_cbl1_e.lower_bound=-.01
    model.reactions.EX_o2_e.lower_bound=-20
       
    return model

# Load the reference model (in MATLAB .mat format)
ref_model = cobra.io.load_matlab_model('/home/omidard/ref_model_dir/marlbr2.mat')
"""
# List of target model file paths (in JSON format)
target_model_dir = '/home/omidard/gapfilled'
target_model_files = [os.path.join(target_model_dir, filename) for filename in os.listdir(target_model_dir) if filename.endswith('.json')]
"""

# List of target model file paths (in JSON format) in gapfilled
target_model_dir = '/home/omidard/output_models_dir'
target_model_files_gapfilled = [os.path.join(target_model_dir, filename) for filename in os.listdir(target_model_dir) if filename.endswith('.json')]

# List of target model file paths (in JSON format) in gapfilled3
target_model_dir_gapfilled3 = '/home/omidard/gapfilled3'
target_model_files_gapfilled3 = [os.path.join(target_model_dir_gapfilled3, filename) for filename in os.listdir(target_model_dir_gapfilled3) if filename.endswith('.json')]

# Filter out files that exist in gapfilled but not in gapfilled3
target_model_files = [file for file in target_model_files_gapfilled if os.path.basename(file) not in [os.path.basename(f) for f in target_model_files_gapfilled3]]




def optimize_with_timeout(model, timeout_seconds):
    try:
        return model.optimize()
    except Exception as e:
        print(f"Optimization failed: {e}")
        return None

# Function to process a target model
def process_target_model(target_model_path):
    # Load the target model (in JSON format)
    target_model = cobra.io.load_json_model(target_model_path)

    # Identify missing reactions (reactions in ref_model but not in target_model)
    missing_reactions = [reaction for reaction in ref_model.reactions if reaction.id not in target_model.reactions]
    # Add all missing reactions to the target model with gene names set to 'GAP'
    missing_reactions2=[]
    for reaction in missing_reactions:
        reaction.gene_reaction_rule = '(GAP)'
        missing_reactions2.append(reaction)
    # Add all missing reactions to the target model
    target_model.add_reactions(missing_reactions2)
    m9(target_model)
    
    # List to store essential reactions
    essential_reactions = []

    # Iterate through missing reactions
    for reaction in missing_reactions2:
        # Backup the original model
        original_model = target_model.copy()

        # Temporarily remove the reaction
        target_model.reactions.get_by_id(reaction.id).remove_from_model()
        m9(target_model)

        # Optimize the model with a timeout
        timeout_seconds = 60  # Set your desired timeout in seconds
        result = None
        try:
            result = optimize_with_timeout(target_model, timeout_seconds)
        except multiprocessing.TimeoutError:
            print(f"Optimization timed out for reaction {reaction.id}")

        # Check if the growth rate is zero (essential reaction)
        if result is None or result.objective_value == 0.0:
            essential_reactions.append(reaction)
            print(target_model.id,' = ',reaction.id)

        # Restore the original model
        target_model = original_model

    # Remove non-essential reactions from the target model
    target_model.remove_reactions([reaction for reaction in missing_reactions2 if reaction not in essential_reactions])

    # Save the modified target model (in JSON format)
    target_model_filename = os.path.basename(target_model_path)
    save_json_model(target_model, f'/home/omidard/gapfilled3/{target_model_filename}')

# Use multiprocessing to parallelize the process
num_cores = 64
pool = multiprocessing.Pool(processes=num_cores)
pool.map(process_target_model, target_model_files)
pool.close()
pool.join()
