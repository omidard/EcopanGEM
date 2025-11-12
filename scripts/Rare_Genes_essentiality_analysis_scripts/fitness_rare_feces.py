import os
import pandas as pd
from cobra.io import load_json_model
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from media import feces_media
import signal

# Define a timeout handler
class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException

# Register the signal function handler
signal.signal(signal.SIGALRM, timeout_handler)

def calculate_knockout_fitness(args):
    model_id, reaction_id, media_function, media_name, wildtype_growth_rate = args

    # Append '.json.json' to the model_id
    model_id = f"{model_id}.json.json"

    # Load the original model
    model_file = os.path.join('/home/omidard/gapfilled_curated', model_id)
    original_model = load_json_model(model_file)

    # Apply media condition
    media_function(original_model)

    # Check if the reaction exists in the model
    if reaction_id in original_model.reactions:
        reaction = original_model.reactions.get_by_id(reaction_id)
        original_lower_bound = reaction.lower_bound
        original_upper_bound = reaction.upper_bound
        try:
            # Start the timer for the knockout simulation
            signal.alarm(8)
            # Simulate knock-out by setting bounds to zero
            reaction.lower_bound = 0
            reaction.upper_bound = 0

            # Set the solver time limit
            original_model.solver.configuration.timeout = 10
            
            # Optimize the model to calculate knock-out growth rate
            knock_out_growth_rate = original_model.optimize().objective_value
            signal.alarm(0)  # Disable the alarm
        except TimeoutException:
            print(f"Simulation for reaction {reaction_id} in model {model_id} exceeded time limit")
            knock_out_growth_rate = None
        finally:
            # Reset reaction bounds to original values
            reaction.lower_bound = original_lower_bound
            reaction.upper_bound = original_upper_bound
    else:
        # If the reaction does not exist, set knock_out_growth_rate to None
        knock_out_growth_rate = None

    # Ensure the knockout growth rate is non-negative and handle None cases
    if knock_out_growth_rate is not None:
        if knock_out_growth_rate < 0:
            knock_out_growth_rate = 0

    # Calculate fitness as percentage decrease in growth rate
    if wildtype_growth_rate is not None and wildtype_growth_rate != 0 and knock_out_growth_rate is not None:
        fitness = ((wildtype_growth_rate - knock_out_growth_rate) / wildtype_growth_rate) * 100
        # Cap the fitness score at 100%
        if fitness > 100:
            fitness = 100
        elif fitness < 0:
            fitness = 0
    else:
        fitness = None

    # Round the values to three decimal places
    wildtype_growth_rate = round(wildtype_growth_rate, 3) if wildtype_growth_rate is not None else None
    knock_out_growth_rate = round(knock_out_growth_rate, 3) if knock_out_growth_rate is not None else None
    fitness = round(fitness, 3) if fitness is not None else None

    # Get GPR (Gene-Protein-Reaction) rules
    gpr = reaction.gene_reaction_rule if reaction_id in original_model.reactions else None

    # Return results
    return {
        'model_id': model_id,
        'reaction_id': reaction_id,
        'wildtype_growth_rate': wildtype_growth_rate,
        'knock_out_growth_rate': knock_out_growth_rate,
        'fitness': fitness,
        'media_name': media_name,
        'gpr': gpr
    }

def get_wildtype_growth_rate(model_id, media_function, media_name):
    model_id_str = f"{model_id}.json.json"
    model_file = os.path.join('/home/omidard/gapfilled_curated', model_id_str)
    model = load_json_model(model_file)
    media_function(model)
    wildtype_growth_rate = model.optimize().objective_value
    return model_id, media_function, media_name, wildtype_growth_rate

def run_simulation():
    # Load fitness_rare DataFrame
    fitness_rare = pd.read_csv('/home/omidard/fitness_rare.csv', index_col=0, dtype=str)

    # Extract models and reactions
    model_ids = fitness_rare.index[921:].tolist()
    reaction_ids = fitness_rare.columns.tolist()

    # Define the m9_aerobic media condition
    media_conditions = [
        (feces_media, 'feces_media')
    ]

    # Process each model independently
    for model_id in tqdm(model_ids, desc="Processing models"):
        wildtype_growth_rates = []

        # Get the wild-type growth rate for m9_aerobic media
        for media_function, media_name in media_conditions:
            wildtype_growth_rates.append(get_wildtype_growth_rate(model_id, media_function, media_name))

        # Create the arguments for starmap
        args = [
            (model_id, reaction_id, media_function, media_name, wildtype_growth_rate)
            for model_id, media_function, media_name, wildtype_growth_rate in wildtype_growth_rates
            for reaction_id in reaction_ids
        ]

        # Define the number of processes to use
        num_processes = min(96, cpu_count() - 1)  # Use up to 96 available CPU cores

        # Initialize an empty list to store the processed data for the current model
        processed_data = []

        # Process fitness_rare using multiprocessing Pool for each media condition
        with Pool(processes=num_processes) as pool:
            for result in tqdm(pool.imap_unordered(calculate_knockout_fitness, args), total=len(args), desc=f"Processing {model_id}"):
                processed_data.append(result)

        # Convert the processed data to a DataFrame
        print(f"Converting processed data for model {model_id} to DataFrame...")
        processed_df = pd.DataFrame(processed_data)

        # Save the DataFrame to a CSV file
        output_file = f'/home/omidard/fitness_results/feces_fitness_{model_id}.csv'
        print(f"Saving DataFrame to {output_file}...")
        try:
            processed_df.to_csv(output_file, index=False)
            print(f"Saved model {model_id} successfully to {output_file}.")
        except Exception as e:
            print(f"Error saving model {model_id}: {e}")

    print("Processing complete and data saved.")

if __name__ == "__main__":
    # Ensure the output directory exists
    os.makedirs('/home/omidard/fitness_results', exist_ok=True)
    run_simulation()
