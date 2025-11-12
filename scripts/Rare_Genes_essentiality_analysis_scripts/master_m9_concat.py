import pandas as pd
import os

# Define the directory and the file pattern
directory = '/home/omidard/fitness_results'
file_pattern = 'm9_fitness_'

# List all files in the directory
all_files = os.listdir(directory)

# Filter files that start with the pattern
filtered_files = [f for f in all_files if f.startswith(file_pattern) and f.endswith('.csv')]

# Create a list to hold dataframes
df_list = []

# Read each file and append to the list
for file in filtered_files:
    file_path = os.path.join(directory, file)
    df = pd.read_csv(file_path)
    df_list.append(df)

# Concatenate dataframes column-wise
master_df = pd.concat(df_list, axis=0)

# Reindex the rows
master_df.reset_index(drop=True, inplace=True)

# Save the concatenated dataframe to a new CSV file
output_path = os.path.join(directory, 'master_m9_knockout.csv')
master_df.to_csv(output_path, index=False)

print(f'Master dataframe saved to {output_path}')
