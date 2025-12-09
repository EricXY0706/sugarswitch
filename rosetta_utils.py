import pyrosetta
from pyrosetta.rosetta.core.scoring.sasa import SasaCalc
from pyrosetta import init, Pose
import pyrosetta.rosetta.core.import_pose as ip
from pyrosetta.toolbox.mutants import mutate_residue

from util import StructureLoader, StructureFileEditor
from Bio.PDB import PDBIO
from io import StringIO

class Rosetta_funcs:

    def __init__(self) -> None:
        self.init = self.initiation()
        self.chains = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

    def initiation(self) -> None:
        init("-mute all")
    
    @staticmethod
    def get_pose(
        structure_file: str
    ) -> Pose:
        '''
        Get the pose object from the structure file.
        
        structure_file: Path to the structure file.
        '''

        structure, _ = StructureLoader.load_structure(structure_file)
        io = PDBIO()
        io.set_structure(structure)
        buf = StringIO()
        io.save(buf)
        pdb_str = buf.getvalue()
        pose = Pose()
        ip.pose_from_pdbstring(pose, pdb_str)
        
        return pose

    def get_SASA(
        self, 
        structure_file: str, 
        cutoff: float = 0.5,
        chain: str = "A",
    ):
        '''
        Get the SASA values of each residue in the structure.
        
        structure_file: Path to the structure file.
        cutoff: SASA cutoff value.
        chain: Chain ID to be calculated.
        '''
        sasa_calc = SasaCalc()
        
        # Get pose
        pose = self.get_pose(structure_file)
        
        # SASA calculation
        sasa_calc.calculate(pose)
        residue_sasa_list = sasa_calc.get_residue_sasa()
        
        # Store the info in sasa_result: [chain_id, res_id, res_name, sasa_value] style
        sasa_result, sasa_values = [], []
        low_sasa_regions = set()
        sasa_index_dict = {}
        for i, sasa in enumerate(residue_sasa_list):
            chain_id = pose.chain(i+1)
            res_id = pose.pdb_info().number(i+1)
            sasa_result.append([chain_id, res_id, sasa])
            sasa_values.append(sasa)
            if chain_id == (self.chains.index(chain) + 1):
                sasa_index_dict[res_id] = sasa
        
        # Retrieve the SASA threshold
        sasa_values.sort(reverse=False)
        sasa_cutoff = sasa_values[round(len(sasa_values) * cutoff)]
        
        # Filter the high SASA regions out
        
        filtered_sasa_regions = [item for item in sasa_result if item[-1] <= sasa_cutoff and item[0] == (self.chains.index(chain) + 1)]
        for item in filtered_sasa_regions:
            low_sasa_regions.add(item[1])
            
        # Save to file
        StructureFileEditor.write_bfactor(pose, sasa_result, structure_file)

        return sasa_cutoff, low_sasa_regions, sasa_index_dict

    def mutate(
        self,
        structure_file: str,
        output_file: str,
        chain_id: str,
        mutate_position: int, 
        mutation: str
    ) -> None:
        '''
        Mutate the structure file.
        
        structure_file: Path to the structure file.
        output_file: Path to the output file.
        chain_id: Chain ID to be mutated.
        mutate_position: Position to be mutated, starting from 1.
        mutation: AA after mutation
        '''

        pose = self.get_pose(structure_file)
        target_pose_index = pose.pdb_info().pdb2pose(chain_id, mutate_position)
        mutate_residue(pack_or_pose=pose, mutant_position=target_pose_index, mutant_aa=mutation, pack_radius=8.0)
        pose.dump_pdb(output_file)