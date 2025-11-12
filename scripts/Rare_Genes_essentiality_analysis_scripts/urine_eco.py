import os
import pandas as pd
from cobra.io import load_json_model
from cobra.flux_analysis import flux_variability_analysis
from concurrent.futures import ProcessPoolExecutor, TimeoutError, as_completed
from tqdm import tqdm

def urine_media(model, exchanges):
    """Modify the specified model's medium to simulate urine-based medium."""
    biomass_reaction_id = 'BIOMASS_Ec_iML1515_core_75p37M'
    model.objective = biomass_reaction_id

    # Set the lower bound of all exchange reactions to zero
    for reaction in model.reactions:
        if reaction.id.startswith('EX_'):
            reaction.lower_bound = 0

    # Set the lower bound for specified exchange reactions in the urine medium
    for ex in exchanges:
        if ex in model.reactions:
            model.reactions.get_by_id(ex).lower_bound = -0.5

    # Set lower bounds for EX_h_e and EX_h2o_e
    if 'EX_h_e' in model.reactions:
        model.reactions.EX_h_e.lower_bound = -1000
    if 'EX_h2o_e' in model.reactions:
        model.reactions.EX_h2o_e.lower_bound = -1000

    return model

def prepare_and_run_fva(model_path, exchange_reactions):
    """Load a model, apply the urine medium modifications, and perform FVA on exchange reactions."""
    model = load_json_model(model_path)
    urine_media(model, exchange_reactions)
    # Get only exchange reactions for FVA
    exchange_rxns = [rxn for rxn in model.reactions if rxn.id in exchange_reactions]
    fva_result = flux_variability_analysis(model, exchange_rxns, processes=1)
    fva_result.columns = [model.id + '_urine_min', model.id + '_urine_max']
    return fva_result

def process_model_with_timeout(model_path, exchange_reactions, timeout=60):
    """Process a single model with a timeout."""
    with ProcessPoolExecutor(max_workers=1) as executor:
        future = executor.submit(prepare_and_run_fva, model_path, exchange_reactions)
        try:
            result = future.result(timeout=timeout)
        except TimeoutError:
            print(f"Model {model_path} exceeded the timeout of {timeout} seconds and was skipped.")
            result = None
    return result

def process_models(model_folder, exchanges, num_models=None, timeout=20):
    """Process a specified number of GEMs using multiprocessing with a timeout for each model."""
    model_files = [os.path.join(model_folder, f) for f in os.listdir(model_folder)]
    if num_models is not None:
        model_files = model_files[:num_models]
    results = []
    with ProcessPoolExecutor(max_workers=64) as executor:
        future_to_model = {executor.submit(process_model_with_timeout, path, exchanges, timeout): path for path in model_files}
        with tqdm(total=len(future_to_model), desc="Processing models") as pbar:
            for future in as_completed(future_to_model):
                res = future.result()
                if res is not None:
                    results.append(res)
                pbar.update(1)
    return pd.concat(results, axis=1)

# List of exchange reactions for urine media simulation
exchange_reactions = ['EX_duri_e','EX_dcyt_e','EX_ade_e','EX_adn_e','EX_cellb_e','EX_ala_B_e','EX_din_e','EX_56dura_e','EX_cgly_e',
 'EX_dgsn_e','EX_g3pc_e','EX_cytd_e','EX_chol_e','EX_xyl__D_e','EX_cyst__L_e','EX_dad_2_e','EX_galt_e',
 'EX_etoh_e','EX_4abut_e','EX_gly_e','EX_glyc_e','EX_gua_e','EX_gsn_e','EX_gal_e','EX_hxan_e',
 'EX_tyr__L_e','EX_phe__L_e','EX_ala__L_e','EX_pro__L_e','EX_thr__L_e','EX_asn__L_e','EX_ile__L_e',
 'EX_his__L_e','EX_ser__L_e','EX_cysi__L_e','EX_ins_e','EX_inost_e','EX_acgam_e','EX_pydxn_e','EX_sbt__D_e',
 'EX_taur_e','EX_sucr_e','EX_thym_e','EX_sarcs_e','EX_thymd_e','EX_rib__D_e','EX_urate_e','EX_cbl1_e',
 'EX_xan_e','EX_urea_e','EX_uri_e','EX_xtsn_e','EX_ura_e','EX_alltn_e','EX_ca2_e','EX_pi_e',
 'EX_cl_e','EX_rbt_e','EX_mg2_e','EX_abt__D_e','EX_cys__L_e','EX_k_e','EX_na1_e','EX_cobalt2_e',
 'EX_csn_e','EX_gln__L_e','EX_cu2_e','EX_leu__L_e','EX_fe2_e','EX_met__L_e','EX_hom__L_e','EX_indole_e','EX_mnl_e',
 'EX_rnam_e','EX_val__L_e','EX_citr__L_e','EX_tmao_e','EX_trp__L_e','EX_ahcys_e','EX_tre_e','EX_acald_e','EX_5aop_e',
 'EX_gbbtn_e','EX_mthgxl_e','EX_5mta_e','EX_amet_e','EX_ala__D_e','EX_mn2_e','EX_26dap_LL_e','EX_ncam_e',
 'EX_fald_e','EX_tcynt_e','EX_pydx_e','EX_4crsol_e','EX_meoh_e','EX_12ppd__R_e','EX_dha_e','EX_co2_e','EX_dad_5_e',
 'EX_metsox_S__L_e','EX_23dappa_e','EX_dms_e','EX_hqn_e','EX_ni2_e','EX_ch4_e','EX_no3_e','EX_acser_e','EX_h2o2_e',
 'EX_5aptn_e','EX_acorn_e','EX_ser__D_e','EX_4gudbutn_e','EX_hg2_e','EX_cd2_e','EX_athr__L_e','EX_nh4_e',
 'EX_5hoxindact_e','EX_acetol_e','EX_progly_e','EX_2ameph_e','EX_cm_e','EX_ttrcyc_e','EX_sel_e','EX_slnt_e',
 'EX_mincyc_e','EX_fe3_e','EX_zn2_e','EX_iad_e','EX_4hphac_e','EX_6apa_e','EX_mobd_e','EX_tungs_e','EX_so4_e','EX_glc__D_e']

# Directory where models are stored
model_folder = '/home/omidard/gapfilled3'

# Run the processing on all models
final_df = process_models(model_folder, exchange_reactions, timeout=20)
final_df.to_csv('/home/omidard/Urine_eco_fva.csv')
