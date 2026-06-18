# SugarSwitch🪄
**SugarSwitch** is an AI-driven pipeline for protein N-glycoengineering, integrating glycosylation site screening, prioritization, and multi-site N-glycosylated sequence design.
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
   # Running the setup process with AI agent assistance to minimize software and hardware compatibility issues is recommended.
   chmod +x setup.sh
   ./setup.sh
   ```
> ✅ **This procedure will automatically**:
> - Create and initiate **sugarswitch** conda environment with Python 3.10
> - Install all the required dependencies
> - Download Boltz2, SaProt, SPIRED, and ESM-LoRA-Gly weights
> - Download and install PyRosetta and DSSP

### Usage
1. Simple Usage
   - **Prefilter**:  
      ➡️ Run the Prefilter to scan and get the detailed information of the potential N-glycosylation sites
      ```Python
      # An example

      python run_sugarswitch.py prefilter \
         --input "input_fasta_file" \
         --input_structure_file "optinal_input_structure_file" \
         --out_dir "output_dir" \
         --chain_id "CHAIN_ID" \
         --functional_hotspots "[1,2,3,'4-10']" \ # This has to be a string with quotation marks
         --gpu_id 0
      ```
      ⚠️ The input file in FASTA format should contain at least one sequence.
      ```txt
      # An example: In this case, test_seq_1 will be chain A-C, and test_seq_2 will be chain D

      >test_seq_1 copies: 3
      MSSQIRQNYSTDVEAAVNSLVNLYLQASYTYLSLGFYFDR...
      >test_seq_2 copies: 1
      WNPPTFSPALLVVTEGDNATFTCSFSNTSESFVLNWYRMS...
      ```
   - **Designer**:  
      ➡️ Run the Designer to design the N-glycosylated sequence with mutations and possible indels given the WT sequence
      ```Python
      # An example: Design 3 sequences with 5 N-glycosylation sites each

      python run_sugarswitch.py designer \
         --input "input_fasta_file" \
         --out_dir "output_dir" \
         --chain_id "CHAIN_ID" \
         --num_designs 3 \
         --num_gly_sites 5 \
         --num_steps 100 \ # Number of hallucination steps
         --add_pll_loss True \ # This makes the sequence more natural but significantly increases computational overhead
         --gpu_id 0
      ```
   
   - **Pipeline**:  
      ➡️ Run the pipeline with both **Prefilter** and **Designer**
      ```Python
      # An example

      python run_sugarswitch.py pipeline \
         --input "input_fasta_file" \
         --input_structure_file "optinal_input_structure_file" \
         --out_dir "output_dir" \
         --chain_id "CHAIN_ID" \
         --functional_hotspots "[1,2,3,'4-10']" \ # This has to be a string with quotation marks
         --num_designs 3 \
         --num_gly_sites 5 \
         --num_steps 100 \ # Number of hallucination steps
         --add_pll_loss True \ # This makes the sequence more natural but significantly increases computational overhead
         --gpu_id 0
      ```
2. Optionally update the configrations for **Prefilter** in `config.py`
   
   - **pipeline configurations**:  
      ⚠️ We recommand not adjusting the **EVCouplings** related configs (EVC, unless you are familiar with the parameters) and empirically pre-defined **glycan chain topological parameters** (bond length, angles, and dihedrals).  
      ✅ `conservation_threshold`, `evc_coupling_threshold`, and `sasa_cutoff` can be appropriately adjusted.
      ```Python
      # An example
      
      pipeline_configs = {
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
         "dihedral_C1_ND2_CG_CB": {"loop": 178.5, "helix": 178.5, "sheet": 178.5},
         "dihedral_C2_C1_ND2_CG": {"loop": 90.0, "helix": 90.0, "sheet": 90.0},
         "dihedral_O5_C1_ND2_CG": {"loop": -95.0, "helix": -95.0, "sheet": -95.0},
      }
      ```
   - **ranker configurations**:  
      ➡️ Manually set the weights for all the factors for glycosylation sites ranking. The higher the value, the more important the factor is.
      ```Python
      # An example

      ranker_configs = {
         "sasa_weight": 1.0, # SASA of the current site. The next 5 values are computed by PyRosetta.
         "gvp_unbind_score_weight": 0.5, # Unbinding level. This is computed by GVP-Bind
         "ddG_weight": 0.3, # ΔΔG of the protein mutant (NXX) and WT. The next 6 values are computed by SPIRED.
         "dTm_weight": 0.3, # ΔTm of the protein mutant (NXX) and WT
         "ddG_S_weight": 0.2, # ΔΔG of the protein mutant (NXS) and WT
         "dTm_S_weight": 0.2, # ΔTm of the protein mutant (NXS) and WT
         "ddG_T_weight": 0.2, # ΔΔG of the protein mutant (NXT) and WT
         "dTm_T_weight": 0.2, # ΔTm of the protein mutant (NXT) and WT
         "mut_score_weight": 0.3, # Mutation score of the protein mutant (NXX). The next 3 values are computed by SaProt.
         "mut_score_S_weight": 0.2, # Mutation score of the protein mutant (NXS).
         "mut_score_T_weight": 0.2, # Mutation score of the protein mutant (NXT).
      }
      ```

## Acknowledgements
1. Hopf, Thomas A., et al. "The EVcouplings Python framework for coevolutionary sequence analysis." Bioinformatics 35.9 (2019): 1582-1584.
2. Chaudhury, Sidhartha, Sergey Lyskov, and Jeffrey J. Gray. "PyRosetta: a script-based interface for implementing molecular modeling algorithms using Rosetta." Bioinformatics 26.5 (2010): 689-691.
3. Chen, Yinghui, et al. "An end-to-end framework for the prediction of protein structure and fitness from single sequence." Nature Communications 15.1 (2024): 7400.
4. Su, Jin, et al. "Saprot: Protein language modeling with structure-aware vocabulary." BioRxiv (2023): 2023-10.
5. Feng, Zhiyong, et al. "ESM-LoRA-Gly: Improved prediction of N-and O-linked glycosylation sites by tuning protein language models with low-rank adaptation (LoRA)." bioRxiv (2025): 2025-08.

## Licence & Citation
**License**: MIT License - See LICENSE file for details  
**Citation**: If you use SugarSwitch in your research, please cite:
```
@article{xxx,
  title={title},
  author={authors},
  journal={journal},
  pages={pages},
  year={2025},
  publisher={publisher}
}
```