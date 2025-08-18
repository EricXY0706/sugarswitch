user_configs = {
    "protein_chain_id": "A",
    "functional_hotspots": [1], # e.g. [1,2,3,"4-10"]
}
# URIC: [58,160,177,229]
# FTL: [155,156,158,159,161,164,165]
# GH: [1,4,8,116,120,41,45,64,172,175,178], 26
# IL2: [42,44,45,65,72,61,62,35,38,43,20,16,19,84,88,91,126,15,123,129,12,19], 20
# NB1: [27,29,30,32,53,54,57,74,100,103,104,106]
# NB2: ["51-58","97-105","26-34"]
# NB3: [29,37,98,"101-105"]

pipeline_configs = {
    "Cb_interaction_threshold": 6.0,
    "num_neighbors_to_shield": 3,
    "evc_min_sequence_distance": 6,
    "evc_theta": 0.8,
    "evc_num_iterations": 100,
    "evc_lambda_h": 0.01,
    "evc_lambda_J": 0.01,
    "evc_num_cpu": 10,
    "conservation_threshold": {"loop": 0.7, "ss": 0.6}, # uricase: (0.65, 0.55)
    "evc_coupling_threshold": 0.5,
    "sasa_cutoff": 0.5,
    "bond_length_C1_ND2": 1.43,
    "angle_C1_ND2_CG": 120.0,
    "angle_C2_C1_ND2": 109.5,
    "angle_O5_C1_ND2": 109.5,
    "dihedral_C1_ND2_CG_CB": {"loop": -120.0, "helix": -180.0, "sheet": 180.0},
    "dihedral_C2_C1_ND2_CG": {"loop": 120.0, "helix": 140.0, "sheet": 100.0},
    "dihedral_O5_C1_ND2_CG": {"loop": -120.0, "helix": -100.0, "sheet": -140.0},
}

basic_configs = {**user_configs, **pipeline_configs}
ranker_configs = {
    "conservation_weight": 1.0,
    "coupling_weight": 0.5,
    "sasa_weight": 1.0,
    "sasa_next1_weight": 0.7,
    "sasa_next2_weight": 0.5,
    "sasa_around_weight": 1.0,
    "sasa_next_weight": 0.6,
    "ddG_weight": 0.3,
    "dTm_weight": 0.3,
    "ddG_S_weight": 0.2,
    "dTm_S_weight": 0.2,
    "ddG_T_weight": 0.2,
    "dTm_T_weight": 0.2,
    "mut_score_weight": 0.3,
    "mut_score_S_weight": 0.2,
    "mut_score_T_weight": 0.2,
}