import os
import pandas as pd
from cobra.io import load_json_model
from cobra.flux_analysis import flux_variability_analysis
from concurrent.futures import ProcessPoolExecutor, TimeoutError, as_completed
from tqdm import tqdm

def serum_media(model, exchanges):
    """Modify the specified model's medium to simulate serum-based medium."""
    biomass_reaction_id = 'BIOMASS_Ec_iML1515_core_75p37M'
    model.objective = biomass_reaction_id

    # Set the lower bound of all exchange reactions to zero
    for reaction in model.reactions:
        if reaction.id.startswith('EX_'):
            reaction.lower_bound = 0

    # Set the lower bound for specified exchange reactions in the serum medium
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
    """Load a model, apply the serum medium modifications, and perform FVA on exchange reactions."""
    model = load_json_model(model_path)
    serum_media(model, exchange_reactions)
    # Get only exchange reactions for FVA
    exchange_rxns = [rxn for rxn in model.reactions if rxn.id in exchange_reactions]
    fva_result = flux_variability_analysis(model, exchange_rxns, processes=1)
    fva_result.columns = [model.id + '_min_serum', model.id + '_max_serum']
    return fva_result

def process_model_with_timeout(model_path, exchange_reactions, timeout=20):
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

# List of exchange reactions for serum media simulation
exchange_reactions = ['EX_duri_e', 'EX_dcyt_e', 'EX_cala_e', 'EX_ade_e', 'EX_adn_e', 'EX_nh4_e', 'EX_ala_B_e',
 'EX_acac_e', 'EX_din_e', 'EX_acon_C_e', 'EX_56dura_e', 'EX_cgly_e', 'EX_56dthm_e', 'EX_dgsn_e', 'EX_g3pc_e',
 'EX_cytd_e', 'EX_chol_e', 'EX_xyl__D_e', 'EX_cyst__L_e', 'EX_dad_2_e', 'EX_galt_e', 'EX_etoh_e', 'EX_4abut_e',
 'EX_glc__D_e', 'EX_gly_e', 'EX_glyc_e', 'EX_gua_e', 'EX_gsn_e', 'EX_for_e', 'EX_gal_e', 'EX_glu__L_e',
 'EX_hxan_e', 'EX_tyr__L_e', 'EX_phe__L_e', 'EX_ala__L_e', 'EX_pro__L_e', 'EX_thr__L_e', 'EX_asn__L_e',
 'EX_man_e', 'EX_ile__L_e', 'EX_fuc__L_e', 'EX_his__L_e', 'EX_ser__L_e', 'EX_cysi__L_e', 'EX_ins_e',
 'EX_ind3ac_e', 'EX_ocdcea_e', 'EX_pnto__R_e', 'EX_inost_e', 'EX_acgam_e', 'EX_nadp_e', 'EX_nadph_e',
 'EX_nmn_e', 'EX_quln_e', 'EX_pydxn_e', 'EX_so3_e', 'EX_ppp9_e', 'EX_pyr_e', 'EX_ribflv_e', 'EX_ppbng_e',
 'EX_sbt__D_e', 'EX_ppi_e', 'EX_taur_e', 'EX_succ_e', 'EX_tsul_e', 'EX_sucr_e', 'EX_srtn_e', 'EX_thym_e',
 'EX_pep_e', 'EX_sarcs_e', 'EX_pser__L_e', 'EX_thymd_e', 'EX_rib__D_e', 'EX_udpg_e', 'EX_urate_e', 'EX_xan_e',
 'EX_urea_e', 'EX_uri_e', 'EX_xtsn_e', 'EX_ura_e', 'EX_trypta_e', 'EX_3hpppn_e', 'EX_23dhb_e', 'EX_3c3hmp_e',
 'EX_3hoxpac_e', 'EX_alltn_e', 'EX_ca2_e', 'EX_3mop_e', 'EX_cl_e', 'EX_rbt_e', 'EX_hxa_e', 'EX_mg2_e',
 'EX_galctn__D_e', 'EX_abt__D_e', 'EX_cys__L_e', 'EX_k_e', 'EX_na1_e', 'EX_s_e', 'EX_cobalt2_e', 'EX_ru5p__D_e',
 'EX_glcn_e', 'EX_csn_e', 'EX_gln__L_e', 'EX_gal1p_e', 'EX_arab__L_e', 'EX_cu2_e', 'EX_glutar_e', 'EX_f_e',
 'EX_2hyoxplac_e', 'EX_pa_EC_e', 'EX_leu__L_e', 'EX_fe2_e', 'EX_S2hglut_e', 'EX_met__L_e', 'EX_hom__L_e',
 'EX_scys__L_e', 'EX_indole_e', 'EX_hcys__L_e', 'EX_xylu__L_e', 'EX_glycogen_e', 'EX_pppn_e', 'EX_mnl_e',
 'EX_psuri_e', 'EX_orot5p_e', 'EX_ppoh_e', 'EX_val__L_e', 'EX_nad_e', 'EX_citr__L_e', 'EX_tma_e', 'EX_tmao_e',
 'EX_trp__L_e', 'EX_ahcys_e', 'EX_tartr__L_e', 'EX_acon_T_e', 'EX_acald_e', 'EX_3sala_e', 'EX_q8h2_e',
 'EX_s7p_e', 'EX_g3p_e', 'EX_acglu_e', 'EX_5aop_e', 'EX_gbbtn_e', 'EX_mthgxl_e', 'EX_5mta_e', 'EX_amet_e',
 'EX_ap4a_e', 'EX_13dpg_e', 'EX_N1aspmd_e', 'EX_1315507', 'EX_maltttr_e', 'EX_ala__D_e', 'EX_mn2_e', 'EX_26dap_LL_e',
 'EX_o2_e', 'EX_4abz_e', 'EX_ncam_e', 'EX_fald_e', 'EX_pi_e', 'EX_tcynt_e', 'EX_pydx_e', 'EX_r5p_e', 'EX_5fthf_e',
 'EX_xylu__D_e', 'EX_34257', 'EX_3hcinnm_e', 'EX_4crsol_e', 'EX_meoh_e', 'EX_12ppd__R_e', 'EX_salc_e', 'EX_co2_e',
 'EX_dad_5_e', 'EX_metsox_S__L_e', 'EX_23dappa_e', 'EX_3uib_e', 'EX_cynt_e', 'EX_h2o_e', 'EX_n8aspmd_e', 'EX_cur_e',
 'EX_dms_e', 'EX_15dap_e', 'EX_hqn_e', 'EX_ni2_e', 'EX_ag_e', 'EX_ch4_e', 'EX_no2_e', 'EX_no3_e', 'EX_41475',
 'EX_42776', 'EX_acser_e', 'EX_skm_e', 'EX_quin_e', 'EX_h2o2_e', 'EX_raffin_e', 'EX_h2s_e', 'EX_gthox_e', 'EX_glu__D_e',
 'EX_5aptn_e', 'EX_acorn_e', 'EX_no_e', 'EX_ser__D_e', 'EX_4gudbutn_e', 'EX_minohp_e', 'EX_hg2_e', 'EX_cd2_e',
 'EX_5hoxindact_e', 'EX_thcur_e', 'EX_tyrp_e', 'EX_pa_EC_e', 'EX_pa_EC_e', 'EX_pa_EC_e', 'EX_pa_EC_e', 'EX_pa_EC_e',
 'EX_rml_e', 'EX_slnt_e', 'EX_progly_e', 'EX_thrp_e', 'EX_peamn_e', 'EX_ps180_e', 'EX_cm_e', 'EX_ttrcyc_e', 'EX_doxrbcn_e',
 'EX_mincyc_e', 'EX_rfamp_e', 'EX_novbcn_e', 'EX_peng_e', 'EX_fe3_e', 'EX_zn2_e', 'EX_fusa_e', 'EX_iad_e', 'EX_frulys_e',
 'EX_nh4_e', 'EX_6apa_e', 'EX_35307', 'EX_mepn_e', 'EX_dhna_e', 'EX_dhnpt_e', 'EX_4h3npacald_e', 'EX_5fthf_e', 'EX_4ahmmp_e',
 'EX_allphn_e', 'EX_aacald_e', 'EX_chols_e', 'EX_ctbt_e', 'EX_pa140_e', 'EX_lipoate_e', 'EX_ethso3_e', 'EX_frulys_e',
 'EX_g3ps_e', 'EX_galctn__L_e', 'EX_mqn8_e', 'EX_suchms_e', 'EX_s7p_e', 'EX_sulfac_e', 'EX_tet_e']

# Directory where models are stored
model_folder = '/home/omidard/gapfilled3'

# Run the processing on all models
final_df = process_models(model_folder, exchange_reactions, timeout=20)
final_df.to_csv('/home/omidard/serum_eco_fva.csv')
