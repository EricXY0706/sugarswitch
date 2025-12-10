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
> - Download PyRosetta

### Usage
1. **Update the configrations in `config.py`**  
1.1 **Must-ToDos**: Manually spefify the protein chain ID to be modified with `protein_chain_id` and the import sites which the modification pipeline should ignore with `functional_hotspots` in **`user_configs`**
   ```Python
   # An Example
   
   user_configs = {
    "protein_chain_id": "A", # Specify the chain ID to be modified and glycosylated
    "functional_hotspots": [1,2,3], # e.g. [1,2,3,"4-10"], if an interval is input, both the lower and the upper bond are included
   }
   ```
