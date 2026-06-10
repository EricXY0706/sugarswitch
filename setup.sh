#!/bin/bash
###############################################################################
# setup.sh — SugarSwitch environment setup
# Running the setup process with AI agent assistance to minimize software and hardware compatibility issues is recommended.

set -uo pipefail

cd "$(dirname "$(readlink -f "$0")")"

echo "=== [SugarSwitch Setup] Initializing environment ==="

# ------------------------------- Configuration ------------------------------
ENV_NAME="sugarswitch"
PY_VERSION="3.10"
export HF_ENDPOINT="${HF_ENDPOINT:-https://hf-mirror.com}"
GH_PROXY="${GH_PROXY:-https://ghfast.top/}"

# ------------------------- 1. Conda environment -----------------------------
if ! command -v conda &>/dev/null; then
    echo "[ERROR] Conda not found. Please install Miniconda/Anaconda first."
    exit 1
fi

echo "[INFO] Initializing conda..."
source "$(conda info --base)/etc/profile.d/conda.sh"

if conda env list | grep -qE "/${ENV_NAME}$|^${ENV_NAME}[[:space:]]"; then
    echo "[WARN] Conda environment '${ENV_NAME}' already exists. Skipping creation."
else
    echo "[INFO] Creating conda environment '${ENV_NAME}' (python ${PY_VERSION})..."
    conda create -n "${ENV_NAME}" python="${PY_VERSION}" -y
fi

conda activate "${ENV_NAME}"

if ! command -v aria2c &>/dev/null; then
    echo "[INFO] Installing aria2 (multi-connection downloader)..."
    conda install -c conda-forge aria2 -y
fi

# ------------------------------ 2. Boltz2 -----------------------------------
echo "[INFO] Installing Boltz2..."
pip install -U "boltz[cuda]"
pip install -U huggingface_hub

mkdir -p boltz_ckpt
hf download --repo-type model boltz-community/boltz-2 --local-dir boltz_ckpt

if [[ -f boltz_ckpt/mols.tar ]]; then
    echo "[INFO] Extracting mols.tar -> boltz_ckpt/mols/ ..."
    tar -xf boltz_ckpt/mols.tar -C boltz_ckpt
else
    echo "[WARN] mols.tar not found in boltz_ckpt/"
fi

# ----------------------- 3. EVCouplings + plmc ------------------------------
echo "[INFO] Installing EVCouplings..."
pip install evcouplings
pip install "setuptools<81"

echo "[INFO] Cloning and building plmc..."
if [[ ! -d plmc ]]; then
    git clone "${GH_PROXY}https://github.com/debbiemarkslab/plmc.git"
fi
( cd plmc && make all-openmp32 )

# -------------------- 4. SaProt + SPIRED + Foldseek -------------------------
echo "[INFO] Downloading SaProt checkpoints..."
mkdir -p SaProt/weights/PLMs
hf download --repo-type model westlake-repl/SaProt_650M_PDB --local-dir SaProt/weights/PLMs

echo "[INFO] Downloading SPIRED model..."
mkdir -p Spired/model
aria2c -x16 -s16 -k1M --file-allocation=none --console-log-level=warn \
    -d Spired/model -o model.zip \
    "https://zenodo.org/records/10589086/files/model.zip?download=1"

if ! command -v unzip &>/dev/null; then
    echo "[ERROR] unzip not found. Please install unzip."
    exit 1
fi
unzip -o Spired/model/model.zip -d Spired/model/

echo "[INFO] Downloading Foldseek binary (official mmseqs.com build)..."
if grep -qm1 avx2 /proc/cpuinfo; then FS_VARIANT="avx2"; else FS_VARIANT="sse41"; fi
mkdir -p SaProt/foldseek
aria2c -x8 -s8 --file-allocation=none --console-log-level=warn \
    -d /tmp -o "foldseek-${FS_VARIANT}.tar.gz" \
    "https://mmseqs.com/foldseek/foldseek-linux-${FS_VARIANT}.tar.gz"
tar -xzf "/tmp/foldseek-${FS_VARIANT}.tar.gz" -C /tmp/
cp /tmp/foldseek/bin/foldseek SaProt/foldseek/foldseek
chmod -R 777 SaProt/foldseek
echo "[INFO] Foldseek installed: $(./SaProt/foldseek/foldseek version 2>/dev/null)"

# -------------------------- 5. PyRosetta + DSSP -----------------------------
echo "[INFO] Installing PyRosetta..."
conda install -c https://conda.graylab.jhu.edu pyrosetta=2024.39+release.59628fb -y

echo "[INFO] Installing DSSP..."
conda install -c ostrokach dssp -y
MKDSSP="$(command -v mkdssp || true)"
if [[ -n "${MKDSSP}" ]]; then
    ln -sf "${MKDSSP}" "$(dirname "${MKDSSP}")/dssp"
    echo "[INFO] Linked dssp -> ${MKDSSP}"
fi

# ------------------------ 6. Python dependencies ----------------------------
echo "[INFO] Installing Python dependencies..."
pip install \
    biopython==1.84 \
    einops==0.8.0 \
    matplotlib==3.10.5 \
    pandas==2.3.1 \
    rdkit==2024.3.2 \
    scikit-learn==1.6.1 \
    seaborn==0.13.2 \
    torch==2.6.0 \
    ml_collections==0.1.1 \
    --no-cache-dir

pip install prody==2.4.1 --no-deps --no-cache-dir
pip install "transformers<5" peft fair-esm --no-cache-dir

echo "[INFO] Pinning numpy==1.26.4 ..."
pip install numpy==1.26.4 --no-cache-dir
NPV="$(python -c 'import numpy; print(numpy.__version__)' 2>/dev/null || echo none)"
if [[ "${NPV}" != "1.26.4" ]]; then
    echo "[WARN] numpy is ${NPV} after pin; forcing a clean reinstall..."
    pip uninstall -y numpy; pip uninstall -y numpy
    pip install numpy==1.26.4 --no-cache-dir --no-deps --force-reinstall
fi
python -c "import numpy; print('[OK] numpy', numpy.__version__)"

echo ""
echo "=== [SugarSwitch Setup] Completed Successfully ==="
echo "Activate the environment using:"
echo "    conda activate ${ENV_NAME}"