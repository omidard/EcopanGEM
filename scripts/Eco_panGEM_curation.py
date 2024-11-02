import os
import cobra
from cobra import Model, Reaction, Metabolite
from multiprocessing import Pool, Manager, Lock, Value
from tqdm import tqdm
import logging

# Directory paths
input_dir = '/home/omidard/gapfilled3'
output_dir = '/home/omidard/gapfilled_curated'
special_model_path = '/home/omidard/iB21_1397.json'

# Dictionary for missing reactions
missing_reactions = {}

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# List of files to process
model_files = [f for f in os.listdir(input_dir) if f.endswith('.json.json')]

# Placeholder for curation instructions
curation_instructions = {
    'R00004': ['PPA'],
    'R00017': {'new_id': 'CCP', 'lower_bound': 0},
    'R00036': ['PPBNGS'],
    'R00084': ['HMBS'],
    'R00104': ['NADK'],
    'R00138': ['PPA2'],
    'R00161': ['FMNAT'],
    'R00178': ['ADMDC'],
    'R00191': ['PDE1'],
    'R00197': ['LDH_D'],
    'R00259': ['ACGS'],
    'R00275': ['SPODM'],
    'R00286': ['UDPGD'],
    'R00289': ['GALUi'],
    'R00339': ['TARTD'],
    'R00351': ['CS'],
    'R00352': {'lower_bound': 0},
    'R00371': ['GLYAT'],
    'R00451': ['DAPDC'],
    'R00476': ['GLYCTO1', 'GLYCTO2', 'GLYCTO3'],
    'R00479': ['ICL'],
    'R00480': ['ASPK'],
    'R00489': ['GLUDC'],
    'R00566': ['ARGDC'],
    'R00577': {'lower_bound': 0},
    'R00586': ['SERAT'],
    'R00601': {'lower_bound': 0},
    'R00619': {'new_id': 'TMDPK', 'lower_bound': 0},
    'R00670': ['ORNDC'],
    'R00693': {'lower_bound': 0},
    'R00709': ['ICDHyr'],
    'R00768': ['GF6PTA'],
    'R00787': ['NTRIR2x'],
    'R00804': {
        'add_gprs': ['F1PP', 'R5PP', 'G1PP'],
        'new_reactions': [
            {'id': 'E4PP', 'reaction': 'e4p_c + h2o_c -> erthrs_c + pi_c', 'name': 'D-Erythrose 4-phosphate phosphatase', 'lower_bound': 0, 'upper_bound': 1000},
            {'id': 'SBT6PP', 'reaction': 'sbt6p_c + h2o_c -> sbt__D_c + pi_c', 'name': 'D-Sorbitol 6-phosphate phosphatase', 'lower_bound': 0, 'upper_bound': 1000},
            {'id': 'MAN1PP', 'reaction': 'man1p_c + h2o_c -> man_c + pi_c', 'name': 'D-Sorbitol 6-phosphate phosphatase', 'lower_bound': 0, 'upper_bound': 1000}
        ]
    },
    'URIt2pp_copy1': {'new_id': 'URIt2pp', 'merge_gpr': 'URIt2pp_copy2'},
    'URIt2pp_copy2': {'new_id': 'URIt2pp', 'merge_gpr': 'URIt2pp_copy1'},
    'R00833': ['MMM'],
    'R00842': {'lower_bound': 0},
    'R00885': ['MAN1PT2'],
    'R00888': ['GMAND'],
    'R00948': ['GLGC'],
    'R00955': ['UGLT'],
    'R00965': ['OMPDC'],
    'R01009': {'new_id': 'DOLPMT', 'reaction': 'gdpmann_c + dolp_c ⇌ gdp_c + dolmanp_c', 'name': 'Dolichyl-phosphate beta-D-mannosyltransferase'},
    'R01015': ['TPI'],
    'R01025': ['CHOLD'],
    'R01049': ['PRPPS'],
    'R01061': {'add_gprs': ['GAPD', 'E4PD']},
    'R01071': ['ATPPRT'],
    'R01072': ['GLUPRT'],
    'R01086': ['ARGSL'],
    'R01092': ['GALKr'],
    'R01123': ['IPDDI'],
    'R01134': ['GMPR'],
    'R01201': ['ACGAMK'],
    'R01213': ['IPPS'],
    'R01224': ['MTHFR2'],
    'R01291': ['RHCCE'],
    'R01394': ['HPYRI'],
    'R01397': ['ASPCT'],
    'R01398': ['OCBT'],
    'R01465': ['THRD'],
    'R01512': ['PGK'],
    'R01518': ['PGM'],
    'R01526': {'new_id': 'RBK_Dr'},
    'R01530': ['A5PISO'],
    'R01540': ['ALTRH'],
    'R01542': {'lower_bound': 0},
    'R01639': ['XYLK'],
    'R01663': {'new_id': 'DCMPDA'},
    'R01708': {'lower_bound': 0},
    'R01714': ['CHORS'],
    'R01717': ['ICHORS'],
    'R01724': ['NAMNPP'],
    'R01737': ['GNK'],
    'R01797': {'lower_bound': 0},
    'R01799': ['DASYN120', 'DASYN140', 'DASYN141', 'DASYN160', 'DASYN161', 'DASYN180', 'DASYN181'],
    'R01801': ['PGSA120', 'PGSA140', 'PGSA141', 'PGSA160', 'PGSA161', 'PGSA180', 'PGSA181'],
    'R01818': ['PMANM'],
    'R01825': ['E4PD'],
    'R01826': ['DDPA'],
    'R01868': ['DHORD2', 'DHORD5'],
    'R01931': ['CYANST'],
    'R02030': {'lower_bound': 0},
    'R02055': ['PSD120', 'PSD140', 'PSD141', 'PSD160', 'PSD161', 'PSD180', 'PSD181'],
    'R02058': ['G1PACT'],
    'R02101': ['TMDS'],
    'R02240': ['DAGK120', 'DAGK140', 'DAGK141', 'DAGK160', 'DAGK161', 'DAGK180', 'DAGK181'],
    'R02291': ['ASAD'],
    'R02301': ['FOMETRi'],
    'R02412': ['SHKK'],
    'R02437': ['RMI'],
    'R02454': ['MANAO'],
    'R02472': ['DPR'],
    'R02555': ['TAGURr'],
    'R02783': ['UAMAGS'],
    'R02921': {'lower_bound': 0},
    'R03033': ['GALCTND'],
    'R03035': ['PTPATi'],
    'R03165': ['UPP3S'],
    'R03231': ['AMAOTr'],
    'R03254': ['KDOPS'],
    'R03313': ['G5SD'],
    'R03317': ['UACMAMO'],
    'R03348': ['NNDPR'],
    'R03350': ['KDOPP'],
    'R03351': ['KDOCT2'],
    'R03387': ['DDGALK'],
    'R03411': ['CFAS160E', 'CFAS160G', 'CFAS180E', 'CFAS180G'],
    'R03425': ['GLYCDx'],
    'R03443': ['AGPR'],
    'R03449': {'lower_bound': 0},
    'R03472': ['AMPMS2'],
    'R03576': ['ALLK'],
    'R03595': ['SELNPS'],
    'R03774': ['RHMND'],
    'R04189': ['SADH'],
    'R04208': ['PRAIS'],
    'R04273': {'lower_bound': 0},
    'R04292': ['QULNS'],
    'R04383': ['KDUI'],
    'R04424': ['MCITD'],
    'R04448': ['HETZK'],
    'R04477': ['ACMAMUT'],
    'R04606': ['LPADSS'],
    'R04638': ['DNTPPA', 'NTPP5'],
    'R04657': ['TDSK'],
    'R05182': {'lower_bound': 0},
    'R05553': ['ADCL'],
    'R05554': ['UGCIAMH'],
    'R05571': ['GLTPD'],
    'R05606': ['MNNH'],
    'R05627': ['UDCPDP'],
    'R05634': ['CDPMEK'],
    'R05636': ['DXPS'],
    'R05637': ['MECDPS'],
    'R05683': ['IDOND', 'IDOND2'],
    'R05688': ['DXPRIi'],
    'R05712': ['NTRIR3pp', 'NTRIR4pp'],
    'R05747': ['ASR'],
    'R05790': {'lower_bound': 0},
    'R06264': {'lower_bound': 0},
    'R06513': ['TDPGDH'],
    'R06836': ['R15BPK'],
    'R06895': ['CPPPGO2'],
    'R07125': ['KG6PDC'],
    'R07229': ['SELR'],
    'R07262': ['DHNCOAT'],
    'R07359': ['NADH9', 'NADH10', 'NADH5'],
    'R07411': ['HEMEOS'],
    'R07414': ['GGPTRCS', 'GLNS'],
    'R07607': ['METSOXR2'],
    'R07644': ['ENTCS'],
    'R07659': ['UDPKAAT'],
    'R08057': ['DGUNC'],
    'R08224': {'lower_bound': 0},
    'R08572': ['GLYCK2'],
    'R08878': ['ALR2', 'DKGLCNR1'],
    'R08991': {'lower_bound': 0},
    'R09248': ['OCTDPS'],
    'R09394': ['CPMPS'],
    'R09395': ['MPTS'],
    'R09497': ['NO3R1bpp'],
    'R09501': ['DMSOR1'],
    'R09541': {'lower_bound': 0},
    'R09641': ['UDPGPT'],
    'R09675': ['TPRDCOAS'],
    'R09837': ['OXDHCOAT', 'REPHACCOAI'],
    'R09838': ['PACCOAE'],
    'R09958': {'lower_bound': 0},
    'R09978': ['CCGS'],
    'R10002': {'lower_bound': 0},
    'R10093': {'lower_bound': 0},
    'R10142': ['TDPADGAT'],
    'R10185': ['ARMEPNS'],
    'R10186': ['RPNTPH'],
    'R10204': ['CPL'],
    'R10205': ['PRCPD'],
    'R10247': ['TYRL', 'THZPSN'],
    'R10303': ['AADDGT'],
    'R10463': {'lower_bound': 0},
    'R10648': {'lower_bound': 0},
    'R10706': ['CITL'],
    'R10901': ['UM3PL', 'UM4PL'],
    'R10939': ['PAI2I'],
    'R10970': ['6D6SFK'],
    'R11183': ['AI2K'],
    'R11186': ['COLIPAKpp'],
    'R11201': {'lower_bound': 0},
    'R11202': {'lower_bound': 0},
    'R11225': ['OPHBDC'],
    'R11335': ['CYTBO3_4pp'],
    'R11479': ['AMPNTAT'],
    'R11581': ['BMOGDS1', 'BMOGDS2', 'BWCOGDS1', 'BWCOGDS2', 'MOGDS'],
    'R11582': ['MOCDS'],
    'R11591': ['LDH_D', 'LDH_D2'],
    'R11593': {'lower_bound': 0},
    'R11600': ['CRNCAL2', 'CRNDCAL2', 'CTBTCAL2'],
    'R11659': ['THZPSN'],
    'R11773': ['PPCSCT'],
    'R11861': {'new_id': 'GGCLUT2'},
    'R11945': ['NADH16pp', 'NADH17pp', 'NADH18pp'],
    'R12035': ['ENLIPIDAt2ex'],
    'R12146': ['GALCTLO'],
    'R12500': ['CDGR'],
    'R12570': ['GTHPi', 'THIORDXi'],
    'R12579': ['NADHPO'],
    'R12635': {'lower_bound': 0},
    'R12925': ['CU1Opp', 'FEROpp'],
    'R13005': ['HEPT3']}


# Special reaction IDs to be prioritized
special_reaction_ids = {'ENLIPIDAt2ex', 'THZPSN', 'ICHORS', 'GLYCTO1'}

# Extract the list of needed reaction IDs from curation instructions
reaction_ids = {reaction for value in curation_instructions.values() if isinstance(value, list) for reaction in value}
reaction_ids -= special_reaction_ids  # Remove special reaction IDs from the main set

# Function to load model from JSON file
def load_model(file_path):
    return cobra.io.load_json_model(file_path)

# Function to save model to JSON file
def save_model(model, file_path):
    cobra.io.save_json_model(model, file_path)

# Function to create new reactions
def create_reaction(reaction_id, reaction_str, reaction_name, lower_bound, upper_bound, gpr, metabolites):
    reaction = Reaction(reaction_id)
    reaction.name = reaction_name
    reaction.lower_bound = lower_bound
    reaction.upper_bound = upper_bound
    
    # Parse the reaction string and add metabolites
    mets = {}
    for part in reaction_str.split(' + '):
        if ':' in part:
            met_id, stoich = part.split(':')
            if met_id in metabolites:
                mets[metabolites[met_id]] = float(stoich)
            else:
                logging.warning(f"Skipping malformed part '{part}' in reaction string for reaction '{reaction_id}'")
    reaction.add_metabolites(mets)
    
    reaction.gene_reaction_rule = gpr
    return reaction

# Function to clean up GPRs
def clean_gpr(gpr):
    # Remove leading "or " and any excess whitespace
    cleaned_gpr = gpr.strip()
    if cleaned_gpr.startswith("or "):
        cleaned_gpr = cleaned_gpr[3:].strip()
    return cleaned_gpr

# Function to extract required reactions from a model
def extract_required_reactions(file, needed_reactions, all_reactions, lock, progress, total):
    if not needed_reactions:  # Stop if all reactions are found
        return
    model = load_model(os.path.join(input_dir, file))
    with lock:
        for reaction in model.reactions:
            if reaction is None or reaction.id is None:
                continue
            if reaction.id in needed_reactions and reaction.id not in all_reactions:
                reaction.gene_reaction_rule = clean_gpr(reaction.gene_reaction_rule)
                all_reactions[reaction.id] = reaction.copy()
                needed_reactions.remove(reaction.id)
                if not needed_reactions:  # Stop if all reactions are found
                    break
        progress.value += 1
        print(f'Extracting reactions: {progress.value}/{total}', end='\r')

# Function to curate model based on instructions
def curate_model(model, all_reactions, instructions, missing_reactions, metabolites):
    for reaction_id, action in instructions.items():
        if reaction_id not in model.reactions:
            continue  # Skip if the target reaction does not exist in the model
        if isinstance(action, dict):
            # Special cases
            if 'new_id' in action:
                # Change reaction ID
                if reaction_id in model.reactions:
                    new_id = action['new_id']
                    if new_id in model.reactions:
                        model.reactions.remove(new_id)
                    reaction = model.reactions.get_by_id(reaction_id)
                    reaction.lower_bound = 0
                    reaction.id = new_id
                else:
                    missing_reactions[reaction_id] = action
            if 'lower_bound' in action:
                if reaction_id in model.reactions:
                    reaction = model.reactions.get_by_id(reaction_id)
                    reaction.lower_bound = action['lower_bound']
            if 'add_gprs' in action:
                if reaction_id in model.reactions:
                    gpr = model.reactions.get_by_id(reaction_id).gene_reaction_rule
                    for add_gpr in action['add_gprs']:
                        if add_gpr in model.reactions:
                            model.reactions.get_by_id(add_gpr).gene_reaction_rule += ' or ' + gpr
                        elif add_gpr in all_reactions:
                            new_reaction = all_reactions[add_gpr].copy()
                            new_reaction.gene_reaction_rule = ''
                            model.add_reactions([new_reaction])
                            model.reactions.get_by_id(add_gpr).gene_reaction_rule += ' or ' + gpr
                        else:
                            missing_reactions[reaction_id] = add_gpr
            if 'new_reactions' in action:
                for new_reaction_info in action['new_reactions']:
                    if reaction_id in model.reactions:
                        gpr = model.reactions.get_by_id(reaction_id).gene_reaction_rule
                    else:
                        gpr = ''
                    if new_reaction_info['id'] in model.reactions:
                        existing_reaction = model.reactions.get_by_id(new_reaction_info['id'])
                        existing_reaction.gene_reaction_rule += ' or ' + clean_gpr(gpr)
                        logging.info(f"Updated GPR for existing reaction '{new_reaction_info['id']}'")
                    else:
                        new_reaction = create_reaction(new_reaction_info['id'], new_reaction_info['reaction'], new_reaction_info['name'],
                                                       new_reaction_info['lower_bound'], new_reaction_info['upper_bound'], clean_gpr(gpr), metabolites)
                        model.add_reactions([new_reaction])
                if reaction_id in model.reactions:
                    model.reactions.remove(reaction_id)
        else:
            # Regular cases
            for target_id in action:
                if target_id in model.reactions:
                    target_reaction = model.reactions.get_by_id(target_id)
                    if reaction_id in model.reactions:
                        reaction = model.reactions.get_by_id(reaction_id)
                        gpr = clean_gpr(reaction.gene_reaction_rule)
                        if 'and' in target_reaction.gene_reaction_rule:
                            target_reaction.gene_reaction_rule += ' and ' + gpr
                        if 'or' in target_reaction.gene_reaction_rule:
                            target_reaction.gene_reaction_rule += ' or ' + gpr
                        else:
                            target_reaction.gene_reaction_rule = gpr
                        model.reactions.remove(reaction_id)
                    else:
                        if reaction_id not in missing_reactions:
                            missing_reactions[reaction_id] = target_id
                elif target_id in all_reactions:
                    new_reaction = all_reactions[target_id].copy()
                    new_reaction.gene_reaction_rule = ''
                    model.add_reactions([new_reaction])
                    if reaction_id in model.reactions:
                        reaction = model.reactions.get_by_id(reaction_id)
                        gpr = clean_gpr(reaction.gene_reaction_rule)
                        if 'and' in new_reaction.gene_reaction_rule:
                            new_reaction.gene_reaction_rule += ' and ' + gpr
                        if 'or' in new_reaction.gene_reaction_rule:
                            new_reaction.gene_reaction_rule += ' or ' + gpr
                        else:
                            new_reaction.gene_reaction_rule = gpr
                        model.reactions.remove(reaction_id)
                    else:
                        if reaction_id not in missing_reactions:
                            missing_reactions[reaction_id] = target_id
                else:
                    if reaction_id not in missing_reactions:
                        missing_reactions[reaction_id] = target_id
    return model

# Function to process a single model file
def process_model_file(args):
    file, all_reactions, instructions, missing_reactions, lock, progress, total, metabolites = args
    try:
        model_path = os.path.join(input_dir, file)
        model = load_model(model_path)
        curated_model = curate_model(model, all_reactions, instructions, missing_reactions, metabolites)
        save_model(curated_model, os.path.join(output_dir, file))
    except Exception as e:
        logging.warning(f"Skipping {file} due to error: {e}")
    finally:
        with lock:
            progress.value += 1
            print(f'Curating models: {progress.value}/{total}', end='\r')

# Function to extract metabolites from a model
def extract_metabolites(model):
    metabolites = {}
    for metabolite in model.metabolites:
        metabolites[metabolite.id] = metabolite
    return metabolites

# Use multiprocessing to extract required reactions
def extract_required_reactions_parallel(model_files, special_model_path):
    with Manager() as manager:
        all_reactions = manager.dict()
        all_metabolites = manager.dict()
        lock = manager.Lock()
        progress = manager.Value('i', 0)
        total = len(model_files)
        needed_reactions = manager.list(reaction_ids)

        # Extract special reactions and metabolites from the special model first
        special_model = load_model(special_model_path)
        with lock:
            for reaction in special_model.reactions:
                if reaction.id in special_reaction_ids:
                    all_reactions[reaction.id] = reaction.copy()
                    special_reaction_ids.remove(reaction.id)
            all_metabolites.update(extract_metabolites(special_model))
        
        with Pool() as pool:
            pool.starmap(extract_required_reactions, [(file, needed_reactions, all_reactions, lock, progress, total) for file in model_files])
        
        # Convert needed_reactions to a regular list to check if it's empty
        needed_reactions_list = list(needed_reactions)
        return dict(all_reactions), dict(all_metabolites), needed_reactions_list

# Extract required reactions
all_reactions, all_metabolites, remaining_needed_reactions = extract_required_reactions_parallel(model_files, special_model_path)

# Use multiprocessing to process all models with a progress bar
if __name__ == '__main__':
    with Manager() as manager:
        lock = manager.Lock()
        progress = manager.Value('i', 0)
        total = len(model_files)
        with Pool() as pool:
            pool.map(process_model_file, [(file, all_reactions, curation_instructions, missing_reactions, lock, progress, total, all_metabolites) for file in model_files])

# Print missing reactions
print("\nMissing Reactions:")
for key, value in missing_reactions.items():
    print(f"{key}: {value}")

# Print any remaining reactions that couldn't be found
if remaining_needed_reactions:
    print("\nRemaining Needed Reactions that couldn't be found:")
    for reaction_id in remaining_needed_reactions:
        print(reaction_id)
