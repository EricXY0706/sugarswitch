"""Write ScanNet-style structure annotations: a PDB whose B-factor column holds
the per-residue binding-site probability, plus a ChimeraX command file to colour
by it. Matches the `annotated_<name>.pdb` / `.cxc` files ScanNet emits.

The B-factor convention matches ScanNet: the raw probability in [0, 1] is written
into the temperature-factor field (PDB columns 61-66) for every atom of a scored
residue; unscored residues (non-query / hetero) get `default`.
"""
from __future__ import annotations

from pathlib import Path


def write_annotated_pdb(src_pdb, out_pdb, prob_by_res: dict, default: float = 0.0) -> None:
    """Copy `src_pdb` to `out_pdb`, overwriting the B-factor with per-residue prob.

    prob_by_res maps (chain_letter, pdb_resnum) -> probability. Residues absent
    from the map get `default`. PDB residue numbering is matched exactly, so the
    output overlays cleanly on the original coordinates.
    """
    src_pdb, out_pdb = Path(src_pdb), Path(out_pdb)
    with open(src_pdb) as f, open(out_pdb, "w") as o:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
                chain = line[21]
                try:
                    resnum = int(line[22:26])
                except ValueError:
                    o.write(line)
                    continue
                p = prob_by_res.get((chain, resnum), default)
                line = line[:60] + ("%6.2f" % p) + line[66:]
            o.write(line)


def write_chimerax_script(cxc_path, pdb_name: str) -> None:
    """Minimal ChimeraX command file: open the annotated PDB and colour by
    B-factor (= binding probability) on a blue-white-red 0..1 scale."""
    cxc_path = Path(cxc_path)
    cxc_path.write_text(
        "\n".join([
            f"open {pdb_name}",
            "color byattribute bfactor palette #4575b4:#ffffff:#d73027 range 0,1 key true",
            "show surface",
            "",
        ])
    )
