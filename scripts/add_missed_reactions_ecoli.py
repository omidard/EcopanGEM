import cobra
import pandas as pd
from multiprocessing import Pool
import os
import re

def process_model(model_id, df, input_path, output_path):
    """
    Process a single metabolic model: Add/modify reactions and metabolites.
    """
    model_file_path = f"{input_path}/{model_id}.json"
    print(model_file_path)    
    # Check if the model file exists
    if not os.path.exists(model_file_path):
        print(f"Model file {model_file_path} does not exist. Skipping...")
        return

    # Load the model
    model = cobra.io.load_json_model(model_file_path)
    
    # Filter the dataframe for reactions related to this model
    model_df = df[df['model_id'] == model_id]
    
    for _, row in model_df.iterrows():
        reaction_id = row['KEGG_ID_y']
        # Check if the reaction is in the model
        if reaction_id in model.reactions:
            reaction = model.reactions.get_by_id(reaction_id)
            # Modify GPR
            new_gpr = f"({reaction.gene_reaction_rule}) and ({row['missing_locus']})"
            reaction.gene_reaction_rule = new_gpr
        else:
            # Add the reaction
            reaction = cobra.Reaction(reaction_id)
            model.add_reactions([reaction])
            reaction.name = row['Sysname']
            reaction.gene_reaction_rule = row['missing_locus']
            # Assuming the reaction string in 'Curated_Reaction' is in a format cobra understands
            reaction.build_reaction_from_string(row['Curated_Reaction'])

    # Save the modified model
    cobra.io.save_json_model(model, f"{output_path}/{model_id}.json")

def parallelize_models(df, input_path, output_path):
    """
    Parallelize the processing of metabolic models.
    """
    model_ids = df['model_id'].unique()
    with Pool(64) as p:
        p.starmap(process_model, [(model_id, df, input_path, output_path) for model_id in model_ids])

def correct_reaction(reaction_str):
    # Remove standalone _c patterns that represent stoichiometry values
    reaction_str = re.sub(r'(\d+)_c', r'\1', reaction_str)
    reaction_str = re.sub(r' n_c', '', reaction_str)
    return reaction_str

if __name__ == "__main__":
    # Load your dataframe
    df = pd.read_csv("/home/omidard/merged_df.csv")

    # Drop rows where 'KEGG_ID_y' or 'Curated_Reaction' is NaN
    df = df.dropna(subset=['KEGG_ID_y', 'Curated_Reaction'])

    # Apply corrections to the Curated_Reaction column
    df['Curated_Reaction'] = df['Curated_Reaction'].apply(correct_reaction)

    # Ensure 'KEGG_ID_y' is of string type
    df['KEGG_ID_y'] = df['KEGG_ID_y'].astype(str)

    # Process models in parallel
    parallelize_models(df, "/home/omidard/gapfilled", "/home/omidard/gapfilled3")
