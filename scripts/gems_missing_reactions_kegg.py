import pandas as pd
import pandas as pd
from bioservices import KEGG




file_path = "/home/omidard/tidy_missing_locci_ecoli.csv"
# Read the file into a DataFrame
df = pd.read_csv(file_path)

# Metabolic keywords
metabolic_keywords = [
    "dehydrogenase", "kinase", "phosphatase", "carboxylase", "synthase", "hydratase", 
    "isomerase", "lyase", "ligase", "oxidase", "reductase", "transferase", "translocase", 
    "mutase", "decarboxylase", "catalase", "amylase", "esterase", "hydrolase", "cyclase",
    "oxidoreductase", "cleavage", "coenzyme", "transporter", "symporter", "antiporter", 
    "channel", "permease", "porin", "pump", "carrier", "biosynthetic", "catabolic", 
    "anabolic", "metabolism", "conversion", "glycolytic", "TCA", "electron transfer", 
    "substrate", "cofactor", "bioconversion", "metabolic pathway"
]

# Non-metabolic keywords
non_metabolic_keywords = [
    "transcriptional", "RNA", "DNA", "ribosomal", "polymerase", "repressor", "activator", 
    "nucleotide", "histone", "transposase", "integrase", "helicase", "replicase", "protein",
    "chaperone", "protease", "ubiquitin", "peptidase", "phosphorylation", "methylation",
    "acetylation", "lipoprotein", "cytoskeletal", "membrane", "flagellar", "cell wall",
    "capsular", "filament", "binding", "recognition", "sensory", "response regulator", 
    "kinase receptor", "adhesion", "secreted", "sensor", "signal","phage"
]

def contains_metabolic_keyword(protein_desc):
    if isinstance(protein_desc, str):  # Only proceed if protein_desc is a string
        for keyword in metabolic_keywords:
            if keyword in protein_desc:
                return True
    return False

def contains_non_metabolic_keyword(protein_desc):
    if isinstance(protein_desc, str):  # Only proceed if protein_desc is a string
        for keyword in non_metabolic_keywords:
            if keyword in protein_desc:
                return True
    return False


# Applying filters
df_metabolic = df[df['protein'].apply(contains_metabolic_keyword) & ~df['protein'].apply(contains_non_metabolic_keyword)]
df_non_metabolic = df[~df['protein'].apply(contains_metabolic_keyword) | df['protein'].apply(contains_non_metabolic_keyword)]


kegg = KEGG()

# Sample list of unique proteins, replace this with df_metabolic.protein.unique()
unique_proteins = df_metabolic['protein'].unique()
# Initialize lists to store results
protein_names = []
kegg_ids = []
ec_numbers = []

def fetch_details_from_kegg(protein_name):
    # Query KEGG for the protein name in the 'genes' database
    result = kegg.find("genes", protein_name)
    
    # If result is not a string, return None for both kegg_id and ec_number
    if not isinstance(result, str):
        print(f"Error fetching details for {protein_name}: Bad response from KEGG")
        return None, None
    
    # Split the result into individual entries
    entries = result.split("\n")

    # Retrieve the KEGG gene identifier for the first entry
    kegg_id = entries[0].split("\t")[0].strip() if entries and len(entries) > 0 else None
    
    ec_number = None
    if kegg_id:  # If KEGG ID is found, then get the EC number.
        try:
            gene_details = kegg.get(kegg_id)  # Fetch details for the gene using its KEGG identifier
            # Check if the gene details contain an EC number
            for line in gene_details.split("\n"):
                if "EC:" in line:
                    ec_number = line.split("EC:")[-1].strip()  # Extract the EC number
        except Exception as e:
            print(f"Error fetching details for {kegg_id}: {e}")
    
    return kegg_id, ec_number


for protein in unique_proteins:
    print(f"Fetching data for: {protein}")
    kegg_id, ec_number = fetch_details_from_kegg(protein)
    
    # Append results to the lists
    protein_names.append(protein)
    kegg_ids.append(kegg_id)
    ec_numbers.append(ec_number)

# Create a DataFrame from the results
df_results = pd.DataFrame({
    'Protein_Name': protein_names,
    'KEGG_ID': kegg_ids,
    'EC_Number': ec_numbers
})

df_results

# Filter rows where EC_Number is not None and not an empty string
filtered_ec = df_results[df_results['EC_Number'].notnull() & (df_results['EC_Number'] != '')]


def fetch_details_from_ec(ec):
    if not ec:
        return "", "", "", "", []
    
    kegg = KEGG()

    enzyme_details = kegg.get(ec)
    
    # Extract relevant information
    parsed_data = kegg.parse(enzyme_details)
    class_name = parsed_data.get("CLASS", "")
    sysname = parsed_data.get("SYSNAME", "")
    reaction = parsed_data.get("REACTION", "")
    all_reac = parsed_data.get("ALL_REAC", "")
    
    # Extract all reference PMIDs
    references = parsed_data.get("REFERENCE", [])
    pmids = [ref.get("REFERENCE").split("[PMID:")[1].split("]")[0] for ref in references if "[PMID:" in ref.get("REFERENCE", "")]
    
    return class_name, sysname, reaction, all_reac, pmids

# Initialize lists to store results
classes, sysnames, reactions, all_reacs, pmid_lists = [], [], [], [], []

for ec in filtered_ec['EC_Number']:
    if ec:
        ec = "ec:" + ec.replace(']','') if not ec.startswith("ec:") else ec  # Ensure EC number has correct format
    class_name, sysname, reaction, all_reac, pmids = fetch_details_from_ec(ec)
    
    classes.append(class_name)
    sysnames.append(sysname)
    reactions.append(reaction)
    all_reacs.append(all_reac)
    pmid_lists.append(", ".join(pmids))

# Update df_results with the fetched information
filtered_ec['Class'] = classes
filtered_ec['Sysname'] = sysnames
filtered_ec['Reaction'] = reactions
filtered_ec['All_Reac'] = all_reacs
filtered_ec['PMIDs'] = pmid_lists

filtered_ec.to_csv("/home/omidard/ecoli_missed_reactions_refined.csv", index=False)


