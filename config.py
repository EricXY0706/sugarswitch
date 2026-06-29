pipeline_configs = {
    "evc_min_sequence_distance": 6,
    "evc_theta": 0.8,
    "evc_num_iterations": 100,
    "evc_lambda_h": 0.01,
    "evc_lambda_J": 0.01,
    "evc_num_cpu": 10,
    "conservation_threshold": 0.8,
    "evc_coupling_threshold": 0.5,
    "rsasa_threshold": 0.3,
    "bond_length_C1_ND2": 1.43,
    "angle_C1_ND2_CG": 120.0,
    "angle_C2_C1_ND2": 109.5,
    "angle_O5_C1_ND2": 109.5,
    "dihedral_C1_ND2_CG_CB": {"loop": 178.5, "helix": 178.5, "sheet": 178.5},
    "dihedral_C2_C1_ND2_CG": {"loop": 90.0, "helix": 90.0, "sheet": 90.0},
    "dihedral_O5_C1_ND2_CG": {"loop": -95.0, "helix": -95.0, "sheet": -95.0},
}

basic_configs = {**pipeline_configs}
ranker_configs = {
    "rSASA_weight": 1.0,
    "ddG_weight": 0.3,
    "dTm_weight": 0.3,
    "ddG_NXST_weight": 0.2,
    "dTm_NXST_weight": 0.2,
    "MutScore_weight": 0.3,
    "MutScore_NXST_weight": 0.2,
}