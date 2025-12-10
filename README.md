# SugarSwitch🪄
**SugarSwitch** is a pipeline for protein N-glycosylation modification.
> 📄 **Paper**: [This is the title for the article](https://scholar.google.com/)  
> 🌐 **Webserver**: [This is the name for the SugarSwitch webserver](https://scholar.google.com/)

---
## Quick start
### Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/EricXY0706/sugarswitch.git
   cd sugarswitch
   ```
2. **Run the automated setup procedure**
   ```bash
   chmod +x setup.sh
   ./setup.sh
   ```
> ⚠️ **Notice**: [Foldseek](https://github.com/steineggerlab/foldseek) binary need to be downloaded from [here](https://drive.google.com/file/d/1B_9t3n_nlj8Y3Kpc_mMjtMdY0OPYa7Re/view) and place it in the `./SaProt/foldseek` folder.  
>
> ✅ **This procedure will automatically**:
> - Create and initiate **sugarswitch** conda environment with Python 3.10
> - Install all the required dependencies
> - Download SaProt and SPIRED weights
> - Download and install PyRosetta and DSSP

### Usage
1. **Update the configrations in `config.py`**
   
   - **Must-ToDos**:  
      ➡️ Manually spefify the protein chain ID to be modified with `protein_chain_id` and the important sites which you think the modification pipeline should **ignore** with `functional_hotspots` in **`user_configs`**
      ```Python
      # An example
      
      user_configs = {
       "protein_chain_id": "A", # Specify the chain ID to be modified and glycosylated.
       "functional_hotspots": [1,2,3], # e.g. [1,2,3,"4-10"], if an interval is input, both the lower and the upper bond are included.
      }
      ```
   - **Optional-ToDos**:
      - **pipeline configurations**:  
        ⚠️ We recommand not adjusting the **EVCouplings** related configs (EVC, unless you are familiar with the parameters) and empirically pre-defined **glycan chain topological parameters** (bond length, angles, and dihedrals).  
        ✅ `Conservation_threshold`, `evc_coupling_threshold`, and `sasa_cutoff` can be appropriately adjusted.
         ```Python
         # An example
         
         pipeline_configs = {
          "Cb_interaction_threshold": 6.0,
          "num_neighbors_to_shield": 3,
          "evc_min_sequence_distance": 6,
          "evc_theta": 0.8,
          "evc_num_iterations": 100,
          "evc_lambda_h": 0.01,
          "evc_lambda_J": 0.01,
          "evc_num_cpu": 10,
          "conservation_threshold": {"loop": 0.5, "ss": 0.5}, # The residues on loop region/other secondary structures with conservation value above the given threshold will be discarded from the prefilter. The higher the thresholds, the less strict the filtering.
          "evc_coupling_threshold": 0.5, # The residues with co-evolving strength above the given threshold will be discarded from the prefilter. The higher the threshold, the less strict the filtering.
          "sasa_cutoff": 0.5, # The residues with SASA values falling in the last n% will be discarded from the prefilter. The higher the threshold, the more strict the filtering.
          "bond_length_C1_ND2": 1.43,
          "angle_C1_ND2_CG": 120.0,
          "angle_C2_C1_ND2": 109.5,
          "angle_O5_C1_ND2": 109.5,
          "dihedral_C1_ND2_CG_CB": {"loop": -120.0, "helix": -180.0, "sheet": 180.0},
          "dihedral_C2_C1_ND2_CG": {"loop": 120.0, "helix": 140.0, "sheet": 100.0},
          "dihedral_O5_C1_ND2_CG": {"loop": -120.0, "helix": -100.0, "sheet": -140.0},
         }
         ```
      - **ranker configurations**:  
        ➡️ Manually set the weights for all the factors for glycosylation sites ranking. The higher the value, the more important the factor is.
        ```Python
        # An example

        ranker_configs = {
          "conservation_weight": 1.0, # Conservation level. The next 2 values are computed by EVCouplings.
          "coupling_weight": 0.5, # Co-evolving strength.
          "sasa_weight": 1.0, # SASA of the current site. The next 5 values are computed by EVCouplings.
          "sasa_next1_weight": 0.7, # SASA of the next site to the current site.
          "sasa_next2_weight": 0.5, # SASA of the next 2nd site to the current site.
          "sasa_around_weight": 1.0, # Mean SASA value of the triad sequon around the current site.
          "sasa_next_weight": 0.6, # Mean SASA value of the triad sequon following the current site.
          "ddG_weight": 0.3, # ΔΔG of the protein mutant (NXX) and WT. The next 6 values are computed by SaProt.
          "dTm_weight": 0.3, # ΔTm of the protein mutant (NXX) and WT
          "ddG_S_weight": 0.2, # ΔΔG of the protein mutant (NXS) and WT
          "dTm_S_weight": 0.2, # ΔTm of the protein mutant (NXS) and WT
          "ddG_T_weight": 0.2, # ΔΔG of the protein mutant (NXT) and WT
          "dTm_T_weight": 0.2, # ΔTm of the protein mutant (NXT) and WT
          "mut_score_weight": 0.3,
          "mut_score_S_weight": 0.2,
          "mut_score_T_weight": 0.2,
         }
        ```
