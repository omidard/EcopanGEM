import cobra
from cobra.io import load_json_model, save_json_model
from cobra.flux_analysis import flux_variability_analysis
import os
import pandas as pd
from tqdm import tqdm
import pathos.multiprocessing as mp
from media import urine_media, feces_media, serum_media

def load_and_apply_media(args):
    model_id, media_function, media_name, output_dir = args
    try:
        model_file = os.path.join('/home/omidard/gapfilled_curated', f"{model_id}.json.json")
        if not os.path.exists(model_file):
            raise FileNotFoundError(f"Model file {model_file} not found.")
        
        model = load_json_model(model_file)
        media_function(model)
        fva_results = flux_variability_analysis(model)
        fva_df = pd.DataFrame(fva_results)
        fva_df.columns = [f"{model_id}_min_{media_name}", f"{model_id}_max_{media_name}"]
        output_file = os.path.join(output_dir, f"{model_id}_{media_name}_fva_results.csv")
        fva_df.to_csv(output_file)
        return fva_df
    except Exception as e:
        print(f"Error processing model {model_id} for {media_name}: {e}")
        return pd.DataFrame()

def process_models(model_ids, media_function, media_name, output_dir):
    args = [(model_id, media_function, media_name, output_dir) for model_id in model_ids]
    results = []
    with mp.Pool(processes=mp.cpu_count()) as pool:
        for result in tqdm(pool.imap_unordered(load_and_apply_media, args), total=len(model_ids)):
            results.append(result)
    return results

def concatenate_results(results):
    concatenated_df = pd.concat(results, axis=1)
    concatenated_df.fillna(0, inplace=True)
    return concatenated_df

def main():
    # Load model IDs
    fitness_rare = pd.read_csv('/home/omidard/fitness_rare.csv', index_col=0, dtype=str)
    model_ids = fitness_rare.index.tolist()
   
    # Define output directory
    output_dir = '/home/omidard/eco_fva_results'
    os.makedirs(output_dir, exist_ok=True)

    # Process models for each media
    media_functions = [(feces_media, "feces"), (serum_media, "serum"), (urine_media, "urine")]
    for media_function, media_name in media_functions:
        print(f"Processing {media_name} media...")
        results = process_models(model_ids, media_function, media_name, output_dir)
        concatenated_df = concatenate_results(results)
        output_file = os.path.join(output_dir, f"all_models_{media_name}_fva_results.csv")
        concatenated_df.to_csv(output_file)
        print(f"Completed processing for {media_name} media. Results saved to {output_file}")

if __name__ == "__main__":
    main()

