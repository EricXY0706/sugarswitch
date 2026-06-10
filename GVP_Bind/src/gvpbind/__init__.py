"""GVP-Bind: per-residue protein binding-site (interface) prediction.

Production model = GVP-GNN over the apo backbone, fused with a SaProt structure-
aware language-model embedding and a heavy-atom graph (trained with a decoupled
contrastive auxiliary loss). Predicts, per residue of a single query chain, the
probability that it lies on a protein-protein interface.
"""
__version__ = "0.2.0"
