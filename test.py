from util import StructureFileEditor, StructureLoader
from Bio.PDB import PDBIO
from io import StringIO
from pyrosetta import init, Pose
import pyrosetta.rosetta.core.import_pose as ip
import pandas as pd


prot = "IL2"
structure_file = f"/sdata2/WORK/Xuyi/{prot}/{prot}.pdb"
scores = pd.read_csv(f"/sdata2/WORK/Xuyi/{prot}/evc/evc_aa_freq.csv")["conservation"].tolist()
data = []
for i, score in enumerate(scores):
    data.append((1, i+1, score*100))

init("-mute all")
structure, _ = StructureLoader.load_structure(structure_file)
io = PDBIO()
io.set_structure(structure)
buf = StringIO()
io.save(buf)
pdb_str = buf.getvalue()
pose = Pose()
ip.pose_from_pdbstring(pose, pdb_str)

StructureFileEditor.write_bfactor(pose, data, f"/sdata2/WORK/Xuyi/{prot}.pdb")
