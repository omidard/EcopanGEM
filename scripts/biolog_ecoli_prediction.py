import os
import cobra
import pandas as pd
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

# Function to apply modifications for M9 media
def m9(model):
    """
    Modify the metabolic model to simulate M9 minimal media conditions.
    """
    for reaction in model.reactions:
        if 'EX_' in reaction.id:
            reaction.lower_bound = 0 

    essential_exchanges = [
        "EX_ca2_e", "EX_cl_e", "EX_co2_e", "EX_cobalt2_e", "EX_cu2_e", "EX_fe2_e", "EX_fe3_e",
        "EX_h_e", "EX_h2o_e", "EX_k_e", "EX_mg2_e", "EX_mn2_e", "EX_mobd_e", "EX_na1_e",
        "EX_tungs_e", "EX_zn2_e", "EX_ni2_e", "EX_sel_e", "EX_slnt_e", "EX_so4_e", "EX_nh4_e",
        "EX_pi_e", "EX_cbl1_e", "EX_o2_e"
    ]

    for ex in essential_exchanges:
        if ex in model.reactions:
            model.reactions.get_by_id(ex).lower_bound = -1000

    model.reactions.get_by_id("EX_cbl1_e").lower_bound = -0.01
    model.reactions.get_by_id("EX_o2_e").lower_bound = -20

    return model

# Function to simulate growth on M9 media with a given carbon source
def simulate_growth_on_m9(model, carbon_source_id):
    """
    Simulate growth on M9 media with a specific carbon source.
    """
    # Ensure model modifications persist
    model = m9(model)

    # Set lower bound of exchange reaction for the carbon source to -1000
    exchange_reaction_id = f'EX_{carbon_source_id}_e'
    if exchange_reaction_id in model.reactions:
        model.reactions.get_by_id(exchange_reaction_id).lower_bound = -1000

    # Run FBA
    solution = model.optimize()

    return solution.objective_value if solution.status == "optimal" else 0

# Function to process a single strain (for parallel execution)
def process_strain(strain, carbon_sources, gem_folder):
    """
    Process a single strain, load the model, and run FBA for each carbon source.
    """
    model_file = os.path.join(gem_folder, f"{strain}")

    # If model does not exist, try GCA/GCF switch
    if not os.path.isfile(model_file):
        if "GCA_" in strain:
            model_file = os.path.join(gem_folder, f"{strain.replace('GCA_', 'GCF_')}")
        elif "GCF_" in strain:
            model_file = os.path.join(gem_folder, f"{strain.replace('GCF_', 'GCA_')}")
        
        if not os.path.isfile(model_file):
            print(f"❌ Skipping (Not Found): {strain}")
            return (strain, "Missing")

    # Load the model
    try:
        model = cobra.io.load_json_model(model_file)
    except Exception as e:
        print(f"⚠️ Error Loading Model: {strain}, {e}")
        return (strain, "Error")

    print(f"✅ Processing Model: {strain}")

    predictions = []
    for carbon_source in carbon_sources:
        if carbon_source == 'No':
            predictions.append(0)
            continue

        # Remove suffixes (_c or _e) from metabolite IDs
        carbon_source_id = carbon_source.rstrip('_ce')

        # Simulate growth on M9 media with the current carbon source
        growth_rate = simulate_growth_on_m9(model, carbon_source_id)

        # Add prediction based on a corrected threshold (0.1)
        predictions.append(1 if growth_rate >= 0.01 else 0)

    return (strain, ';'.join(map(str, predictions)))

# Function to add predictions to Biolog dataframe using multiprocessing
def add_predictions_to_biolog(biolog_df, gem_folder):
    """
    Run the FBA simulations in parallel and store the results.
    """
    # Prepare input for parallel processing
    input_data = [(row['Strain'], row['Met_Id'].split(';'), gem_folder) for _, row in biolog_df.iterrows()]

    # Use multiprocessing pool to process models in parallel
    with Pool(processes=cpu_count()) as pool:
        results = list(tqdm(pool.starmap(process_strain, input_data), total=len(input_data), desc="Processing Strains"))

    # Store results in dataframe
    strain_to_prediction = {strain: prediction for strain, prediction in results}

    # Ensure only strains present in the results are included in the final dataframe
    biolog_df = biolog_df[biolog_df['Strain'].isin(strain_to_prediction.keys())]
    biolog_df['Prediction'] = biolog_df['Strain'].map(strain_to_prediction)

    # Save the updated Biolog dataframe
    output_path = "/home/omidard/biolog_data_with_predictions_panGEM_paper.csv"
    biolog_df.to_csv(output_path, index=False)

    # Summary output
    missing_count = sum(1 for _, pred in results if pred == "Missing")
    error_count = sum(1 for _, pred in results if pred == "Error")
    print(f"\nProcessing complete. Skipped {missing_count} models (not found), {error_count} models (loading errors).")

# Path to the folder containing modified GEMs
gem_folder = '/home/omidard/gapfilled_curated'
# Path to the Biolog dataframe
biolog_df_path = "/home/omidard/biolog_panGEM.csv"

# Load Biolog dataframe
biolog_df = pd.read_csv(biolog_df_path)

# Filter Biolog dataframe to only keep rows with "M9 minimal media no C source" in the 'Media' column
biolog_df = biolog_df[biolog_df['Media'] == 'M9 minimal media no C source']

# Add predictions based on growth simulation results to Biolog dataframe
add_predictions_to_biolog(biolog_df, gem_folder)
