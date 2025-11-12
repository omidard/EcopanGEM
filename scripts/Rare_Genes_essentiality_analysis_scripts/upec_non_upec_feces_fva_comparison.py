import cobra
from cobra.io import load_json_model, save_json_model
from cobra.flux_analysis import flux_variability_analysis
import os
import pandas as pd
from media import urine_media, feces_media, serum_media

# Define model file paths
upec_model_file = '/home/omidard/gapfilled_curated/1355101.3.json.json'
non_upec_model_file = '/home/omidard/gapfilled_curated/1328859.4.json.json'

# Load models
upec_model = load_json_model(upec_model_file)
non_upec_model = load_json_model(non_upec_model_file)

# Apply urine media constraints
feces_media(upec_model)
feces_media(non_upec_model)

# Perform FVA
upec_fva = flux_variability_analysis(upec_model)
non_upec_fva = flux_variability_analysis(non_upec_model)

# Convert FVA results to DataFrames
upec_fva_df = pd.DataFrame(upec_fva)
non_upec_fva_df = pd.DataFrame(non_upec_fva)

# Define output directory
output_dir = '/home/omidard/eco_fva_results_upec'
os.makedirs(output_dir, exist_ok=True)

# Save models after applying urine media constraints
upec_model_output = os.path.join(output_dir, 'upec_model_feces_media.json')
non_upec_model_output = os.path.join(output_dir, 'non_upec_model_feces_media.json')
save_json_model(upec_model, upec_model_output)
save_json_model(non_upec_model, non_upec_model_output)

# Save FVA results
upec_fva_output = os.path.join(output_dir, 'upec_feces_fva_results.csv')
non_upec_fva_output = os.path.join(output_dir, 'non_upec_feces_fva_results.csv')
upec_fva_df.to_csv(upec_fva_output)
non_upec_fva_df.to_csv(non_upec_fva_output)

print("FVA analysis completed and results saved.")
