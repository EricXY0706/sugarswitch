# GVP-Bind — protein binding-site prediction

Predicts, **per residue**, the probability that it lies on a protein–protein
**binding interface**, from a *single* (apo) chain — **no MSA, no complex, no
docking, no folding**.

The model is a GVP-GNN (geometric vector perceptron graph net) over the backbone,
fused with two extra signals:

- a **SaProt** structure-aware language-model embedding (foldseek 3Di + sequence), and
- a **heavy-atom graph** of the query chain,

trained with a decoupled-contrastive (DCL) auxiliary loss. One checkpoint is
bundled — `checkpoints/gvpbind_saprot_atg_dcl.ckpt` — and is the default; no flag
needed.

## How good is it?

The headline is the *capability*: competitive interface prediction **from a
single apo chain, with no MSA and no complex folding**.

The numbers below are from our **internal evaluation** (apo input, single chain;
per-chain median ROC-AUC). They are recorded here for orientation, not as a
published result — reproduce them before citing. On the shared sets the baselines
reproduce their own papers' values, which is our sanity check that the harness is
fair.

| benchmark | GVP-Bind | PeSTo | ScanNet (noMSA) | note |
|---|---|---|---|---|
| post-2022 temporal holdout (536 chains) | 0.855 | 0.828 | 0.790 | clean / leakage-controlled; Wilcoxon p≈1e-5 vs PeSTo |
| ScanNet∩PeSTo test set (670 chains) | 0.947 | 0.927 | 0.881 | both baselines reproduce their papers |
| MaSIF-site transient set (~50) | 0.858 | 0.850 | — | p≈0.03 vs PeSTo; head-to-head fair, but our PeSTo regen is below its paper's 0.92 (input-provenance caveat) |
| AlphaFold-Multimer dimers (24) | 0.951 | 0.963 | — | statistical tie, p≈0.57 — on par with AFM, no MSA/folding |

Read this as: **on our leakage-controlled / external benchmarks it matches or
exceeds PeSTo and ScanNet, and is on par with AlphaFold-Multimer on dimer
interfaces.** Important caveat: on the older **PPBS** splits (chains ≤2015, which
likely sit inside PeSTo's training set) it only **ties** PeSTo on this metric —
so we do *not* headline PPBS, and the clean benchmarks above are the honest
comparison.

> **Use it for "target prediction":** high probability = likely interaction site;
> **low probability = interaction-quiet surface** (e.g. candidate regions for
> modifications that shouldn't disrupt a partner interface). Combine with solvent
> accessibility for that use — the model tells you about *interfaces*, not burial.

---

## Install

```bash
conda create -n gvpbind python=3.10 && conda activate gvpbind
cd GVP-Bind
pip install -e '.[saprot]'      # core + the SaProt embedding backend (transformers 4.28)
```

Core deps: `torch`, `pytorch-lightning`, `torchmetrics`, `numpy`, `biopython`,
`gemmi`. CUDA is used if present (SaProt 650M runs much faster on GPU; CPU works
for a single small chain).

### External assets (not bundled — required by the SaProt embedding)

The bundled model needs a SaProt embedding, which needs two things that are **too
large / not pip-installable** and must be obtained once:

1. **foldseek binary** — https://github.com/steineggerlab/foldseek (single static
   binary).
2. **SaProt_650M_AF2 weights** — https://github.com/westlake-repl/SaProt
   (a HuggingFace-format model directory, ~2.5 GB).

Point the code at them with environment variables (or the `--saprot-model-dir` /
`--foldseek-bin` flags):

```bash
export SAPROT_MODEL_DIR=/path/to/SaProt_650M_AF2     # dir with config.json + pytorch_model.bin
export FOLDSEEK_BIN=/path/to/foldseek
```

---

## Quickstart

Input is `path/to/structure.pdb_<CHAIN>` — the structure and the chain to score.

```bash
gvpbind-predict examples/ferritin_2ffx_chainA.pdb_A --predictions_folder out/
```

Produces in `out/`:

| file | contents |
|------|----------|
| `predictions_*.csv` | per-residue `Binding site probability` |
| `annotated_*.pdb`   | same probability written into the **B-factor** column |
| `annotated_*.cxc`   | ChimeraX script: open it to colour the structure by probability |

View it: `chimerax out/annotated_*.cxc` (or open the PDB in PyMOL and `spectrum b`).

---

## How the model treats the input

The model is **apo / query-chain-only**: it scores the residues of the chain you
name, using *only* that chain's structure (it was trained with the partner chain
dropped). If you pass a multi-chain complex it simply slices to the query chain
first — passing the lone chain gives the identical result.

This means it does **not** need the partner to be present and its probabilities
are meaningful on a bare monomer. Two read-outs:

- default — calibrated per-residue interface probability;
- `--rank-normalize` — within-chain percentile in `[0,1]`, handy when you only
  care about the *ranking* of candidate sites.

---

## Options

- `--checkpoint PATH`     model checkpoint (default: bundled `gvpbind_saprot_atg_dcl.ckpt`).
- `--saprot-model-dir D`  SaProt weights dir (or `$SAPROT_MODEL_DIR`).
- `--foldseek-bin PATH`   foldseek binary (or `$FOLDSEEK_BIN`).
- `--temperature T`       soften over-confident probabilities: `sigmoid(logit / T)`, `T>1`.
- `--rank-normalize`      emit within-chain rank in `[0,1]` instead of probability.
- `--device cpu|cuda`     (auto-detects).
- `--name RUN`            output filename tag.
- `--predictions_folder`  output directory.

## What's inside (for reference)

```
src/gvpbind/
├── model/   gvp_layers.py, atom_encoder.py, scannetpp.py  # GVP-GNN + atom encoder + per-residue head
├── task/    multitask.py, ddg_regression.py               # LightningModule wrappers (checkpoint loading)
├── infer/   saprot_embed.py, atomgraph_infer.py           # inline SaProt embedding + heavy-atom graph
├── data/    dockground.py (PDB→arrays), dataset.py (k-NN graph)
└── cli/     predict.py, _annotate.py
checkpoints/gvpbind_saprot_atg_dcl.ckpt                    # the bundled model (~10 MB; SaProt weights are external)
```

This is the **prediction-only** subset of the GVP-Bind research project; training,
labeling, and benchmarking code are not included. The checkpoint is small because
the heavy SaProt weights live outside it (loaded at run time, see above).
