#!/bin/bash
set -euo pipefail

echo "=== [Sugarswitch Setup] Initializing environment ==="

# -------------------------------
# 1. Conda environment preparation
# -------------------------------
if ! command -v conda &> /dev/null; then
    echo "[ERROR] Conda not found. Please install Miniconda/Anaconda first."
    exit 1
fi

echo "[INFO] Initializing conda..."
source "$(conda info --base)/etc/profile.d/conda.sh"

ENV_NAME="sugarswitch"
PY_VERSION="3.10"

if conda env list | grep -q "^$ENV_NAME "; then
    echo "[WARN] Conda environment '$ENV_NAME' already exists. Skipping creation."
else
    echo "[INFO] Creating conda environment '$ENV_NAME'..."
    conda create -n "$ENV_NAME" python="$PY_VERSION" -y
fi

conda activate "$ENV_NAME"

# -------------------------------
# 2. Boltz2 installation
# -------------------------------
echo "[INFO] Installing Boltz2..."
export HF_ENDPOINT="https://hf-mirror.com"

pip install -U "boltz[cuda]"
pip install -U huggingface_hub

mkdir -p boltz_ckpt
huggingface-cli download \
    --resume-download \
    --repo-type model \
    --local-dir boltz_ckpt \
    --local-dir-use-symlinks False \
    boltz-community/boltz-2

if [[ -f boltz_ckpt/mols.tar ]]; then
    mkdir -p boltz_ckpt/mols
    tar -xvf boltz_ckpt/mols.tar -C boltz_ckpt/mols
else
    echo "[WARN] mols.tar not found in boltz_ckpt/"
fi

# -------------------------------
# 3. Install EVCouplings + plmc
# -------------------------------
echo "[INFO] Installing EVCouplings..."
pip install evcouplings

echo "[INFO] Cloning and building plmc..."
if [[ ! -d plmc ]]; then
    git clone https://github.com/debbiemarkslab/plmc.git
fi

(
    cd plmc
    make all-openmp32
    cd ..
)

# -------------------------------
# 4. SaProt + Spired checkpoints
# -------------------------------
echo "[INFO] Downloading SaProt checkpoints..."
mkdir -p SaProt/weights/PLMs

huggingface-cli download \
    --resume-download \
    --repo-type model \
    --local-dir SaProt/weights/PLMs \
    --local-dir-use-symlinks False \
    westlake-repl/SaProt_650M_PDB

echo "[INFO] Downloading Spired model..."
mkdir -p Spired/model

wget -c "https://zenodo.org/records/10589086/files/model.zip?download=1" \
    -O Spired/model/model.zip

if ! command -v unzip &> /dev/null; then
    echo "[ERROR] unzip not found. Please install unzip."
    exit 1
fi

unzip -o Spired/model/model.zip -d Spired/model/

echo "[WARN] Please manually download Foldseek binary (see SaProt/foldseek/README.md)."

# -------------------------------
# 5. PyRosetta + DSSP
# -------------------------------
echo "[INFO] Installing PyRosetta..."
conda config --add channels https://conda.graylab.jhu.edu
conda install pyrosetta=2024.39+release.59628fb -y

echo "[INFO] Installing DSSP..."
conda install -c ostrokach dssp -y

# -------------------------------
# 5. PyRosetta + DSSP
# -------------------------------
echo "[INFO] Installing LigandMPNN..."
git clone https://github.com/dauparas/LigandMPNN.git
cd LigandMPNN
mkdir -p ./model_params
wget -q https://files.ipd.uw.edu/pub/ligandmpnn/solublempnn_v_48_002.pt -O ./model_params/solublempnn_v_48_002.pt
wget -q https://files.ipd.uw.edu/pub/ligandmpnn/solublempnn_v_48_010.pt -O ./model_params/solublempnn_v_48_010.pt
wget -q https://files.ipd.uw.edu/pub/ligandmpnn/solublempnn_v_48_020.pt -O ./model_params/solublempnn_v_48_020.pt
wget -q https://files.ipd.uw.edu/pub/ligandmpnn/solublempnn_v_48_030.pt -O ./model_params/solublempnn_v_48_030.pt
cd ..

# -------------------------------
# 6. Python dependencies
# -------------------------------
echo "[INFO] Installing Python dependencies..."

pip install \
    biopython==1.84 \
    einops==0.8.0 \
    matplotlib==3.10.5 \
    numpy==1.26.4 \
    pandas==2.3.1 \
    rdkit==2025.3.5 \
    scikit-learn==1.6.1 \
    seaborn==0.13.2 \
    torch==2.8.0 \
    prody==2.4.1 \
    ml_collections==0.1.1 \
    --no-cache-dir

echo ""
echo "=== [Sugarswitch Setup] Completed Successfully ==="
echo "Activate the environment using:"
echo "    conda activate sugarswitch"