import os
import pandas as pd
from cobra.io import load_json_model
from media import m9_aerobic  # Ensure this is properly defined and accessible
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import signal

# Timeout exception class
class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException

# Register the timeout handler
signal.signal(signal.SIGALRM, timeout_handler)

# Define all GEM names
model_ids = [
    "656416.3", "550691.3", "331112.6", "GCF_000005845.2",
    "656441.3", "656447.3", "550685.3", "656383.3",
    "656433.3", "1355100.3", "656372.3", "656444.3"
]

# Base directory for models
models_folder = '/home/omidard/gapfilled_curated'

# Output directory for results
output_file = '/home/omidard/knock_out_validation_combined_complete.csv'

def calculate_knockout_fitness(args):
    """Calculate fitness for a single reaction knockout."""
    model_id, reaction_id, wildtype_growth_rate = args
    model_path = os.path.join(models_folder, f"{model_id}.json.json")
    model = load_json_model(model_path)
    
    # Apply M9 media
    m9_aerobic(model)

    # Check if reaction exists
    if reaction_id in model.reactions:
        reaction = model.reactions.get_by_id(reaction_id)
        original_lb, original_ub = reaction.lower_bound, reaction.upper_bound

        try:
            signal.alarm(8)  # Set timeout
            reaction.lower_bound = 0
            reaction.upper_bound = 0
            knockout_growth_rate = model.optimize().objective_value
            signal.alarm(0)  # Cancel timeout
        except TimeoutException:
            knockout_growth_rate = None
        finally:
            reaction.lower_bound, reaction.upper_bound = original_lb, original_ub
    else:
        knockout_growth_rate = None

    # Calculate fitness
    fitness = None
    if wildtype_growth_rate and knockout_growth_rate is not None:
        fitness = round((wildtype_growth_rate - knockout_growth_rate) / wildtype_growth_rate * 100, 3)

    # Get GPR rule
    gpr = model.reactions.get_by_id(reaction_id).gene_reaction_rule if reaction_id in model.reactions else None

    return {
        "model_id": model_id,
        "reaction_id": reaction_id,
        "wildtype_growth_rate": round(wildtype_growth_rate, 3),
        "knock_out_growth_rate": round(knockout_growth_rate, 3) if knockout_growth_rate else None,
        "fitness": fitness,
        "gpr": gpr
    }

def get_wildtype_growth_rate(model_id):
    """Calculate the wildtype growth rate for a model."""
    model_path = os.path.join(models_folder, f"{model_id}.json.json")
    model = load_json_model(model_path)
    m9_aerobic(model)
    return model.optimize().objective_value

def run_knockout_analysis():
    """Run knockout analysis for all models and reactions."""
    results = []

    # Iterate over all models
    for model_id in tqdm(model_ids, desc="Processing models"):
        try:
            # Load the model and get reaction IDs
            model_path = os.path.join(models_folder, f"{model_id}.json.json")
            model = load_json_model(model_path)
            reaction_ids = [rxn.id for rxn in model.reactions]

            # Get wildtype growth rate
            wildtype_growth_rate = get_wildtype_growth_rate(model_id)

            # Prepare arguments for multiprocessing
            args = [(model_id, rxn_id, wildtype_growth_rate) for rxn_id in reaction_ids]

            # Use multiprocessing for knockouts
            num_processes = min(cpu_count() - 1, 96)
            with Pool(processes=num_processes) as pool:
                for result in tqdm(pool.imap_unordered(calculate_knockout_fitness, args), total=len(args), desc=f"Knockouts for {model_id}"):
                    results.append(result)
        except Exception as e:
            print(f"Error processing model {model_id}: {e}")

    # Combine results into a DataFrame
    results_df = pd.DataFrame(results)

    # Save the results to a CSV file
    results_df.to_csv(output_file, index=False)
    print(f"Knockout analysis completed. Results saved to {output_file}")

if __name__ == "__main__":
    run_knockout_analysis()
