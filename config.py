user_configs = {
    "name": "protein_name",
    "protein_chain_id": "A",
    "functional_hotspots": [1,2,3,"4-10"],
}

pipeline_configs = {
    "Cb_interaction_threshold": 6.0,
    "num_neighbors_to_shield": 3,
    "evc_min_sequence_distance": 6,
    "evc_theta": 0.8,
    "evc_num_iterations": 100,
    "evc_lambda_h": 0.01,
    "evc_lambda_J": 0.01,
    "evc_num_cpu": 10,
    "conservation_threshold": {"loop": 0.5, "ss": 0.5},
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