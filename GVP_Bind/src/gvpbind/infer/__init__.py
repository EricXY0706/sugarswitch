"""On-the-fly feature computation for inference on *novel* structures.

The training pipeline reads SaProt PLM embeddings (``--esm-cache-dir``) and the
heavy-atom graph (``--atomgraph-cache-dir``) from pre-built caches. Those caches
do not exist for structures generated at inference time (e.g. Boltz-2 refolds in
the proseek loop). This subpackage recomputes both features *inline* from a bare
PDB, matching the training recipe, so a SaProt+atomgraph checkpoint can run with
no cache.

Public entry points:
    saprot_embed.compute_saprot_embedding(pdb, chain, ...) -> [L, 1280] float
    atomgraph_infer.compute_atomgraph(pdb, chain, pdb_resnum_order) -> dict
"""
