import os
import cobra
#from cobra.io import load_json_model
from Bio import SeqIO
import pandas as pd
from multiprocessing import Pool
"""
def process_model(model_file):
    # Lists to store results
    missing_loci = []

    # Load model
    model_path = os.path.join(model_dir, model_file)
    model = cobra.io.load_json_model(model_path)

    # Get model gene IDs
    model_genes = {gene.id for gene in model.genes}

    # Read the corresponding genome file
    genome_file = model_file.replace('.json.json', '.gbk')
    genome_path = os.path.join(genome_dir, genome_file)
    genome = SeqIO.parse(genome_path, "genbank")

    # Check genes in genome not in model
    for feature in genome.features:
        if feature.type == "CDS":
            locus_tag = feature.qualifiers.get("locus_tag")[0]
            if locus_tag not in model_genes:
                product = feature.qualifiers.get("product", ["Unknown"])[0]
                missing_loci.append((model.id, locus_tag, product))

    return missing_loci
"""

def process_model(model_file):
    # Lists to store results
    missing_loci = []

    # Load model
    model_path = os.path.join(model_dir, model_file)
    model = cobra.io.load_json_model(model_path)

    # Get model gene IDs
    model_genes = {gene.id for gene in model.genes}

    # Read the corresponding genome file
    genome_file = model_file.replace('.json.json', '.gbk')
    genome_path = os.path.join(genome_dir, genome_file)
    genome_iterator = SeqIO.parse(genome_path, "genbank")
    genome_record = next(genome_iterator)  # Get the first record from the iterator

    # Check genes in genome not in model
    for feature in genome_record.features:
        if feature.type == "CDS":
            locus_tag = feature.qualifiers.get("locus_tag")[0]
            if locus_tag not in model_genes:
                product = feature.qualifiers.get("product", ["Unknown"])[0]
                missing_loci.append((model.id, locus_tag, product))

    return missing_loci



genome_dir = "/home/omidard/prokka_genomes"
model_dir = "/home/omidard/output_models_dir"

# Initialize the pool with 64 workers
with Pool(64) as p:
    results = p.map(process_model, os.listdir(model_dir))  # [:2] for testing, remove to process all

# Flatten the list of lists
flat_results = [item for sublist in results for item in sublist]

# Create DataFrame
df = pd.DataFrame(flat_results, columns=["model_id", "missing_locus", "gene_product"])

# Save DataFrame to CSV
df.to_csv("/home/omidard/ecoli_missing_locus_tags.csv", index=False)
