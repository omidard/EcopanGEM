"""
import pandas as pd

# Specify the path to the CSV file
file_path = "/home/omidard/miss_lucci_df.csv"

# Read the file into a DataFrame
df = pd.read_csv(file_path)
# Create a DataFrame containing only rows with 'hypothetical protein' in the 'Gene_Description' column
df_hypothetical = df[df['Gene_Description'].str.contains('hypothetical protein', case=False, na=False)]
# Create a DataFrame excluding the rows with 'hypothetical protein' in the 'Gene_Description' column
df_without_hypothetical = df[~df['Gene_Description'].str.contains('hypothetical protein', case=False, na=False)]
df_hypothetical.to_csv("/home/omidard/ecoli_hypothetical_proteins.csv", index=False)


def extract_info(s):
    parts = str(s).strip().split(']')
    data_dict = {}
    for part in parts:
        if '=' in part:
            key, value = part.replace('[', '').split('=', 1)
            data_dict[key.strip()] = value.strip()
    return pd.Series(data_dict)

# Apply the function to each row in the dataframe
df_info = df_without_hypothetical['Gene_Description'].apply(extract_info)

# Joining the new columns to the original dataframe
df_tidy = pd.concat([df_without_hypothetical, df_info], axis=1)

# Drop the original 'Gene_Description' column (optional)
df_tidy = df_tidy.drop('Gene_Description', axis=1)
df_tidy.to_csv("/Users/omidard/Desktop/tidy_missing_locci_ecoli.csv", index=False)



import pandas as pd
import concurrent.futures

# Specify the path to the CSV file
file_path = "/home/omidard/miss_lucci_df.csv"

# Read the file into a DataFrame
df = pd.read_csv(file_path)

# Filtering operations
df_hypothetical = df[df['Gene_Description'].str.contains('hypothetical protein', case=False, na=False)]
df_without_hypothetical = df[~df['Gene_Description'].str.contains('hypothetical protein', case=False, na=False)]
df_hypothetical.to_csv("/home/omidard/ecoli_hypothetical_proteins.csv", index=False)

def extract_info(s):
    parts = str(s).strip().split(']')
    data_dict = {}
    for part in parts:
        if '=' in part:
            key, value = part.replace('[', '').split('=', 1)
            data_dict[key.strip()] = value.strip()
    return data_dict

# Apply the function in parallel using concurrent.futures
with concurrent.futures.ProcessPoolExecutor() as executor:
    data_dicts = list(executor.map(extract_info, df_without_hypothetical['Gene_Description']))

df_info = pd.DataFrame(data_dicts)

# Joining the new columns to the original dataframe
df_tidy = pd.concat([df_without_hypothetical.reset_index(drop=True), df_info], axis=1)

# Drop the original 'Gene_Description' column (optional)
df_tidy = df_tidy.drop('Gene_Description', axis=1)
df_tidy.to_csv("/Users/omidard/Desktop/tidy_missing_locci_ecoli.csv", index=False)
"""
import pandas as pd
import concurrent.futures

# Specify the path to the CSV file
file_path = "/home/omidard/ecoli_missing_locus_tags.csv"

# Read the file into a DataFrame
df = pd.read_csv(file_path)

# Drop rows containing 'transposase' in the 'gene_product' column
df = df[~df['gene_product'].str.contains('transposase', case=False, na=False)]

# Filtering operations
df_hypothetical = df[df['gene_product'].str.contains('hypothetical protein', case=False, na=False)]
df_without_hypothetical = df[~df['gene_product'].str.contains('hypothetical protein', case=False, na=False)]
df_hypothetical.to_csv("/home/omidard/ecoli_hypothetical_proteins.csv", index=False)

def extract_info(s):
    parts = str(s).strip().split(']')
    data_dict = {}
    for part in parts:
        if '=' in part:
            key, value = part.replace('[', '').split('=', 1)
            data_dict[key.strip()] = value.strip()
    return data_dict

# Apply the function in parallel using concurrent.futures
with concurrent.futures.ProcessPoolExecutor() as executor:
    data_dicts = list(executor.map(extract_info, df_without_hypothetical['gene_product']))

df_info = pd.DataFrame(data_dicts)

# Joining the new columns to the original dataframe
df_tidy = pd.concat([df_without_hypothetical.reset_index(drop=True), df_info], axis=1)

# Drop the original 'gene_product' column (optional)
df_tidy = df_tidy.drop('gene_product', axis=1)
df_tidy.to_csv("/home/omidard/tidy_missing_locci_ecoli.csv", index=False)

