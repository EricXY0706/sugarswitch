import logging
import os
import shutil
import subprocess
import re
import time
from typing import Set, List, Tuple
from tqdm import tqdm
import tarfile
import yaml

from Bio import SeqIO
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select, Structure, DSSP
import pandas as pd
import numpy as np
from math import degrees
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler
import requests
from requests.auth import HTTPBasicAuth

from config import basic_configs

SS_TAG = {
    "-": ["loop", "loop"], 
    "S": ["bend", "loop"],
    "T": ["turn", "loop"],
    "B": ["bridge", "loop"],
    "H": ["a-helix", "helix"],
    "G": ["310-helix", "helix"],
    "I": ["pi-helix", "helix"],
    "P": ["k-helix", "helix"],
    "E": ["strand", "sheet"],
    }

AA_INTERACTIONS = {
    "H_bond": ["R", "K", "D", "E", "S", "T", "N", "Q", "Y"],
    "Hydrophobic": ["A", "V", "I", "L", "M", "F", "W", "Y"],

}

CHAIN_IDS = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

logger = logging.getLogger(__name__)

class FlowList(list):
    pass

def flow_list_representer(dumper, data):
    return dumper.represent_sequence("tag:yaml.org,2002:seq", data, flow_style=True)

yaml.add_representer(FlowList, flow_list_representer, Dumper=yaml.SafeDumper)

class FastaLoader:

    @staticmethod
    def get_sequence(sequence_file, chain_id):

        chain_id = CHAIN_IDS.index(chain_id) + 1
        seqs = {rec.id: (str(rec.seq), int(rec.description.split(" ")[-1].split(":")[-1])) for rec in SeqIO.parse(sequence_file, "fasta")}
        s = None
        total_num = 0
        for k, (seq, num) in seqs.items():
            total_num += num
            if total_num < chain_id:
                continue
            else:
                s = seq
                break
        return s
        
class StructureLoader:
    
    @staticmethod
    def cif_to_pdb(structure_file):
        
        # Load cif structure file
        filename = structure_file.split(".")[0]
        parser = MMCIFParser()
        structure = parser.get_structure("structure", structure_file)
        
        # Write the structure into pdb structure file
        pdb_io = PDBIO()
        pdb_io.set_structure(structure)
        pdb_io.save(f"{filename}_temp.pdb")
        
        return f"{filename}_temp.pdb"
    
    @staticmethod
    def load_structure(structure_file: str):
        
        if structure_file.endswith("cif"):
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure("structure", structure_file)

        elif structure_file.endswith("pdb"):
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("structure", structure_file)

        return structure, structure_file
    
    def get_sequences(structure_file: str):

        structure, _ = StructureLoader.load_structure(structure_file)
        chain_seq_dict = {}
        for model in structure:
            for chain in model:
                seq = "".join([residue.resname for residue in chain])
                chain_seq_dict[chain.id] = seq
        
        return chain_seq_dict
    
    @staticmethod
    def get_ca_coords(structure_file: str, chain_id: str):

        structure, filename = StructureLoader.load_structure(structure_file)
        ca_coords = []
        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    for residue in chain:
                        if "CA" in residue:
                            ca_coords.append(residue["CA"].coord)
        
        return np.array(ca_coords)
    
    @staticmethod
    def get_chainids(structure_file: str):
        
        structure, filename = StructureLoader.load_structure(structure_file)
        chain_ids = set()
        for model in structure:
            for chain in model:
                chain_ids.add(chain.id)
        
        return chain_ids
    
    @staticmethod
    def get_mean_bfactor(structure_file, chain_ids_list):
        
        # Load structure
        structure, filename = StructureLoader.load_structure(structure_file)
        
        # Enumerate the chains and get mean bfactors
        structure_bfactors = []
        for chain_id in chain_ids_list:
            chain_bfactors = []
            for res in structure[0][chain_id]:
                for atom in res:
                    chain_bfactors.append(atom.bfactor)
            if chain_bfactors:
                chain_mean_bfactor = sum(chain_bfactors) / len(chain_bfactors)
            else:
                chain_mean_bfactor = None
            structure_bfactors.append(chain_mean_bfactor)
        
        if structure_bfactors:
            structure_mean_bfactor = sum(structure_bfactors) / len(structure_bfactors)
        else:
            structure_mean_bfactor = None
        
        return structure_mean_bfactor
    
    @staticmethod
    def remove_temp_file(structure_file):
        
        # Get file name
        filename = structure_file.split(".")[0]
        
        # Delete temp file
        if structure_file.endswith("cif"):
            os.system(f"rm -f {filename}_temp.pdb")
    
    @staticmethod
    def get_secondary_structure(structure_file, chain_id):

        """
        Run DSSP and return a dictionary mapping (chain_id, residue_id) to secondary structure type.
        
        Args:
            pdb_file (str): Path to the PDB file.

        Returns:
            Dict[(chain_id, res_id), ss]: A dictionary where ss is one of 'H', 'E', 'C', etc.
        """
        structure, filename = StructureLoader.load_structure(structure_file)
        model = structure[0]

        # Run DSSP
        dssp = DSSP(model, structure_file)

        # Build the dictionary
        ss_dict = {}
        for key in dssp.keys():
            chain, res_id = key
            if chain == chain_id:
                ss_dict[(chain, res_id[1])] = dssp[key][2]

        return ss_dict

class StructureFileEditor:
    
    @staticmethod
    def write_bfactor(pose, valuelist, filename):
        
        # valuelist: [chain_id, res_id, bfactor]
        value_map = {(chain_id, res_id): bfactor for chain_id, res_id, bfactor in valuelist}
        
        # Enumerate the pose and change the bfactors
        for i in range(1, pose.size()+1):
            chain_id = pose.chain(i)
            res_id = pose.pdb_info().number(i)
            bfactor = value_map.get((chain_id, res_id), 0.00)
            for j in range(1, pose.residue(i).natoms()+1):
                pose.pdb_info().bfactor(i, j, bfactor)
        
        # Save to file
        pose.dump_pdb(filename)
    
    @staticmethod
    def write_score_as_bfactor(pose, structure_file, chain_id, df_file):
        
        df = pd.read_csv(df_file)
        sites = df["Site"].tolist()
        chain_seq_dict = StructureLoader.get_sequences(structure_file)
        data_pose = []
        for i, cs in enumerate(chain_seq_dict.items()):
            chain, seq = cs
            data = [(i+1, s+1, df.loc[df["Site"] == s+1, "Borda_score"].values[0] * 100 if chain == chain_id and (s+1) in sites else 0.) for s in range(len(seq))]
        data_pose.extend(data)
        StructureFileEditor.write_bfactor(pose, data_pose, structure_file)

class MsaFileGenerator:
    
    def __init__(self, input_fasta_file) -> None:
        """
        MSA file generator using ColabFold API.
        Adopted from ColabFold: https://github.com/sokrypton/ColabFold
        """
        self.username = "example_user"
        self.password = "example_password"
        self.user_agent = "colabfold/1.5.5"
        self.host_url = "https://api.colabfold.com"
        # self.host_url = "https://protenix-server.com/api/msa"
        self.email = ""
        self.tqdm_bar_format = "{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} estimate remaining: {remaining}]"
        self.seq_id_pairs = {rec.id: (str(rec.seq), int(rec.description.split(" ")[-1].split(":")[-1])) for rec in SeqIO.parse(input_fasta_file, "fasta")}

    def _submit_to_mmseqs2_service(
        self,
        seq: str,
        mode: str,
        use_pairing: bool = False,
    ):
        headers = {}
        headers["User-Agent"] = self.user_agent
        submission_endpoint = "ticket/pair" if use_pairing else "ticket/msa"
        error_count = 0
        while True:
            try:
                res = requests.post(
                    url=f"{self.host_url}/{submission_endpoint}",
                    data={
                        "q": seq,
                        "mode": mode,
                        "email": self.email,
                    },
                    timeout=6.02,
                    headers=headers,
                    auth=HTTPBasicAuth(self.username, self.password)
                )
            except requests.exceptions.Timeout:
                continue
            except Exception as e:
                error_count += 1
                logger.warning(f"Error while fetching result from MSA server. Retrying ... ({error_count}/5)")
                time.sleep(5)
                if error_count > 5:
                    raise e
                continue
            break
        try:
            out = res.json()
        except ValueError:
            logger.warning(f"Server didn't reply with json: {res.text}")
            out = {"status": "ERROR"}
        return out
    
    def _status(
        self,
        ID: str,
    ):
        headers = {}
        headers["User-Agent"] = self.user_agent
        error_count = 0
        while True:
            try:
                res = requests.get(
                    url=f"{self.host_url}/ticket/{ID}",
                    timeout=6.02,
                    headers=headers,
                    auth=HTTPBasicAuth(self.username, self.password)
                )
            except requests.exceptions.Timeout:
                continue
            except Exception as e:
                error_count += 1
                logger.warning(f"Error while fetching result from MSA server. Retrying ... ({error_count}/5)")
                time.sleep(5)
                if error_count > 5:
                    raise e
                continue
            break
        try:
            out = res.json()
        except ValueError:
            logger.warning(f"Server didn't reply with json: {res.text}")
            out = {"status": "ERROR"}
        return out
    
    def _download_from_mmseqs2_service(
        self,
        ID: str,
        path: str,
    ):
        headers = {}
        headers["User-Agent"] = self.user_agent
        error_count = 0
        while True:
            try:
                res = requests.get(
                    url=f"{self.host_url}/result/download/{ID}",
                    timeout=6.02,
                    headers=headers,
                    auth=HTTPBasicAuth(self.username, self.password)
                )
            except requests.exceptions.Timeout:
                continue
            except Exception as e:
                error_count += 1
                logger.warning(f"Error while fetching result from MSA server. Retrying ... ({error_count}/5)")
                time.sleep(5)
                if error_count > 5:
                    raise e
                continue
            break
        with open(path, "wb") as f:
            f.write(res.content)
        
    def run_mmseqs2(
        self,
        x,
        prefix,
        use_env = True,
        use_filter = True,
        use_template = False,
        filter = None,
        use_pairing = False,
        pairing_strategy = "greedy",
    ):
        seqs = [x] if isinstance(x, str) else x

        if filter is not None:
            use_filter = filter
        if use_filter:
            mode = "env" if use_env else "all"
        else:
            mode = "env-nofilter" if use_env else "nofilter"
        if use_pairing:
            use_template = False
            use_env = False
            mode = ""
            if pairing_strategy == "greedy":
                mode = "pairgreedy"
            elif pairing_strategy == "complete":
                mode = "paircomplete"
        os.makedirs(prefix, exist_ok=True)
        tar_gz_file = f"{prefix}/out.tar.gz"
        REDO = True
        if not os.path.exists(tar_gz_file):
            TIME_ESTIMATE = 100
            with tqdm(total=TIME_ESTIMATE, bar_format=self.tqdm_bar_format) as pbar:
                while REDO:
                    pbar.set_description("SUBMIT")
                    out = self._submit_to_mmseqs2_service(
                        seq=seqs,
                        mode=mode,
                    )
                    while out["status"] in ["UNKNOWN", "RATELIMIT"]:
                        sleep_time = 60
                        logger.error(f"MMSeqs service is not available. Retrying in {sleep_time} seconds ...")
                        time.sleep(sleep_time)
                        out = self._submit_to_mmseqs2_service(
                            seq=seqs,
                            mode=mode,
                        )
                    if out["status"] == "ERROR":
                        raise Exception(f"MMSeqs API is giving errors. Please confirm your input is a valid protein sequence. If error persists, please try again an hour later.")
                    if out["status"] == "MAINTENANCE":
                        raise Exception(f"MMSeqs API is undergoing maintenance. Please try again later.")
                    
                    ID, TIME = out["id"], 0
                    pbar.set_description(out["status"])
                    while out["status"] in ["UNKNOWN", "RUNNING", "PENDING"]:
                        t = 60
                        logger.error(f"Sleeping for {t}s. Reason: {out['status']}")
                        time.sleep(t)
                        out = self._status(ID)
                        pbar.set_description(out["status"])
                        if out["status"] == "RUNNING":
                            TIME += t
                        pbar.n = min(99, int(100 * TIME / (30.0 * 60)))
                        pbar.refresh()
                    if out["status"] == "COMPLETE":
                        pbar.n = 100
                        pbar.refresh()
                        REDO = False
                    if out["status"] == "ERROR":
                        REDO = False
                        raise Exception(f"MMSeqs API is giving errors. Please confirm your input is a valid protein sequence. If error persists, please try again an hour later.")
                self._download_from_mmseqs2_service(
                    ID=ID,
                    path=tar_gz_file,
                )
                with tarfile.open(tar_gz_file, "r:gz") as tar:
                    tar.extractall(os.path.dirname(tar_gz_file))
                lines = open(f"{prefix}/uniref.a3m", "r").readlines()[:-1]
                pats = [re.compile(rf".*{k}") for k, _ in self.seq_id_pairs.items()]
                idx = [[i for i, line in enumerate(lines) if pats[j].match(line.strip())] for j in range(len(pats))]
                idx.append([len(lines)])
                for i, _ in enumerate(range(len(idx) - 1)):
                    start_idx = idx[i][0]
                    end_idx = idx[i + 1][0]
                    first_line = lines[start_idx][1:] if start_idx != 0 else lines[start_idx]
                    with open(f"{prefix}/uniref_{i+1}.a3m", "w") as f:
                        f.write(first_line)
                        f.writelines(lines[start_idx + 1: end_idx])

class MsaFileEditor:
    
    @staticmethod
    def load_msa_seqs(msa_file):
        
        # Load the MSA sequences from a FASTA file into a dictionary
        msa_seqs = {}
        for seq_rec in SeqIO.parse(msa_file, "fasta"):
            msa_seqs[seq_rec.description] = str(seq_rec.seq)
        return msa_seqs
    
    @staticmethod
    def write_msa_record(description, sequence, f):
        
        # Write a sequence record to the output FASTA file
        f.write(f'>{description}\n')
        f.write(f'{sequence}\n')
    
    @staticmethod
    def find_seq_change(old_seq, new_seq):
        
        # Align the two seqs
        aligned_seq_1 = old_seq + "-"*(len(new_seq)-len(old_seq))
        aligned_seq_2 = new_seq
        
        # Find the gaps for insertion
        gaps_index = next((i for i, (a, b) in enumerate(zip(aligned_seq_1, aligned_seq_2)) if a != b), len(old_seq))
        
        return gaps_index, len(new_seq)-len(old_seq)
    
    @staticmethod
    def update_msa_file_single(msa_file_in, msa_file_out, old_seq, new_seq):
        
        # Load the original MSA seqs
        msa_seqs = MsaFileEditor.load_msa_seqs(msa_file_in)
        
        # Find the gaps indices
        gaps_index, gaps_length = MsaFileEditor.find_seq_change(old_seq, new_seq)
        assert gaps_length >= 0
        
        # Update the MSA file
        open(msa_file_out, "w")
        with open(msa_file_out, "a") as f:
            for desc, seq in msa_seqs.items():
                
                # Query seq
                if desc == "101":
                    MsaFileEditor.write_msa_record(desc, new_seq, f)
                
                # MSA seqs
                else:
                    seq = "".join([c for c in list(seq) if c.isupper() or c == "-"])
                    seq_ = seq[:gaps_index] + ("-" * max(0, gaps_length)) + seq[gaps_index:]
                    assert len(seq_) == len(new_seq), f"Error in {desc} alignment"
                    MsaFileEditor.write_msa_record(desc, seq_, f)
                    
    @staticmethod
    def update_msa_file_multiple(msa_file_in, msa_file_out, new_seq, mutations):
        
        # Update the MSA file with a new sequence and adjust gaps only for insertions.
        # mutations = {index: n_ins}
        msa_seqs = MsaFileEditor.load_msa_seqs(msa_file_in)

        with open(msa_file_out, "w") as f:
            for desc, seq in msa_seqs.items():
                if desc == "101":
                    MsaFileEditor.write_msa_record(desc, new_seq, f)
                else:
                    seq = "".join([c for c in seq if c.isupper() or c == "-"])
                    new_aligned_seq = list(seq)     
                    for index, n_ins in sorted(mutations.items(), reverse=True):
                        if index < len(new_aligned_seq):
                            if n_ins >0:
                                new_aligned_seq.insert(index, "-"*n_ins)
                        else:
                            new_aligned_seq.append("-")
                    
                    new_aligned_seq = "".join(new_aligned_seq)

                    assert len(new_aligned_seq) == len(new_seq), f"Error in {desc} alignment"
                    MsaFileEditor.write_msa_record(desc, new_aligned_seq, f)

def update_infer(
        input_fasta_file: str,
        output_dir: str,
        suffix: str = "",
) -> str:
    """
    Update the infer result dir with the input fasta file.
    """
    if not os.path.exists(input_fasta_file):
        raise FileNotFoundError(f"Input fasta file `{input_fasta_file}` not found.")
    filename = f"{basic_configs['name']}{suffix}"
    msa = MsaFileGenerator(input_fasta_file=input_fasta_file)
    with open(input_fasta_file, "r") as f:
        query = f.read()
    msa.run_mmseqs2(x=query, prefix=f"{output_dir}/msa{suffix}")
    
    # Boltz prediction
    seqs = {rec.id: (str(rec.seq), int(rec.description.split(" ")[-1].split(":")[-1])) for rec in SeqIO.parse(input_fasta_file, "fasta")}
    chain_ptr = 0
    out = []
    for name, (seq, num) in seqs.items():
        ids = FlowList(CHAIN_IDS[chain_ptr : chain_ptr + num])
        out.append({
            "protein": {
                "id": ids,
                "sequence": seq,
                "msa": f"{output_dir}/msa{suffix}/uniref_{len(out)+1}.a3m",
            }
        })
        chain_ptr += num

    with open(f"{output_dir}/{filename}.yaml", "w") as f:
        yaml.safe_dump({"version": 1, "sequences": out}, f, sort_keys=False, indent=2)
    
    cmd = [
        "boltz", "predict", f"{output_dir}/{filename}.yaml",
        "--cache", f"./boltz_ckpt",
        "--output_format", "pdb",
        "--out_dir", output_dir,
    ]
    subprocess.run(cmd, check=True)
    subprocess.run(f"mv {output_dir}/boltz_results_{filename}/predictions/{filename}/{filename}_model_0.pdb {output_dir}/{filename}.pdb", shell=True)
    subprocess.run(f"rm -r {output_dir}/boltz_results_{filename}", shell=True)
    subprocess.run(f"rm -f {output_dir}/{filename}.yaml", shell=True)
    shutil.rmtree(f"{output_dir}/msa{suffix}")
    
    return f"{output_dir}/{filename}.pdb"
    
class GlycanMover:

    def __init__(self,
        bond_length: float = 1.43,
        angle_C1: float = 120.0,
        angle_C2: float = 109.5,
        angle_O5: float = 109.5,
        dihedral_C1: float = 178.5,
        dihedral_C2: float = 90.0,
        dihedral_O5: float = -95.0
        ) -> None:
        """
        Initialize default bond geometry values for glycan placement.
        Atoms series: NAG@C2/O5 -> NAG@C1 -> Asn@ND2 -> Asn@CG -> Asn@CB

        Sets:
            self.bond_length (float): Ideal bond length in Ångstroms.
                NAG@C1-Asn@ND2, default is 1.43Å
            self.angle (float): Ideal bond angle (degrees).
                C1: NAG@C1-Asn@ND2-Asn@CG, default is 120.0°
                C2: NAG@C2-NAG@C1-Asn@ND2, default is 109.5°
                O5: NAG@O5-NAG@C1-Asn@ND2, default is 109.5°
            self.dihedral (float): Ideal dihedral angle (degrees).
                C1: NAG@C1-Asn@ND2-Asn@CG-Asn@CB, default is -120.0° (loop) / 180° (helix)
                C2: NAG@C2-NAG@C1-Asn@ND2-Asn@CG, default is 120.0°
                O5: NAG@O5-NAG@C1-Asn@ND2-Asn@CG, default is -120.0°
        """
        
        self.bond_length = bond_length
        self.angle_c1 = angle_C1
        self.angle_c2 = angle_C2
        self.angle_o5 = angle_O5
        self.dihedral_c1 = dihedral_C1
        self.dihedral_c2 = dihedral_C2
        self.dihedral_o5 = dihedral_O5

    def _get_atom_coords(
        self,
        structure,
        chain_id: str,
        res_id: int,
        atom_names: List[str]
    ) -> np.ndarray:
        """
        Retrieve coordinates of specific atoms from a residue.

        Args:
            structure: Bio.PDB structure object.
            chain_id (str): Chain ID (e.g., "A").
            res_id (int): Residue ID.
            atom_names (list[str]): List of atom names to extract.

        Returns:
            np.ndarray: Array of shape (N, 3) with atom coordinates.
        """
        model = structure[0]
        residue = model[chain_id][res_id]
        coords = []
        for atom_name in atom_names:
            coords.append(residue[atom_name].coord)
        return np.array(coords)

    def _place_atom(
        self,
        origin: np.ndarray,
        bond_length: float,
        angle_deg: float,
        dihedral_deg: float,
        prev1: np.ndarray,
        prev2: np.ndarray,
        prev3: np.ndarray
    ) -> np.ndarray:
        """
        Place a new atom in 3D space given bond length, angle, and dihedral.

        Args:
            origin (np.ndarray): Origin point (not used here but kept for consistency).
            bond_length (float): Desired bond length to prev1.
            angle_deg (float): Bond angle between new_atom-prev1-prev2.
            dihedral_deg (float): Dihedral angle defined by new_atom-prev1-prev2-prev3.
            prev1 (np.ndarray): First atom (closest to the new atom).
            prev2 (np.ndarray): Second atom (middle point).
            prev3 (np.ndarray): Third atom (farthest).

        Returns:
            np.ndarray: The 3D coordinates of the new atom.
        """

        theta = np.deg2rad(angle_deg)
        phi = np.deg2rad(dihedral_deg)

        b1 = prev1 - prev2
        b2 = prev2 - prev3
        b1 /= np.linalg.norm(b1)
        b2 /= np.linalg.norm(b2)

        n1 = np.cross(b1, b2)
        n1 /= np.linalg.norm(n1)
        n2 = np.cross(n1, b1)
        n2 /= np.linalg.norm(n2)

        d = bond_length * (np.cos(theta)*(-b1) + np.sin(theta)*(np.cos(phi)*n2 + np.sin(phi)*n1))
        new_atom = prev1 + d
        return new_atom
    
    def _build_frame(
        self,
        origin: np.ndarray,
        pt1: np.ndarray,
        pt2: np.ndarray
    ) -> np.ndarray:
        """
        Build a local orthonormal frame from three points.

        Args:
            origin (np.ndarray): Reference point.
            pt1 (np.ndarray): Defines x-axis direction from origin.
            pt2 (np.ndarray): Helps define plane and z-axis.

        Returns:
            np.ndarray: 3x3 rotation matrix representing the frame (column-wise axes).
        """
        x = pt1 - origin
        x /= np.linalg.norm(x)
        z = np.cross(x, pt2 - origin)
        z /= np.linalg.norm(z)
        y = np.cross(z, x)
        return np.vstack([x, y, z]).T

    def _dihedral_angle(
        self,
        p0: np.ndarray,
        p1: np.ndarray,
        p2: np.ndarray,
        p3: np.ndarray
    ) -> float:
        """
        Compute the dihedral angle (in degrees) defined by four points.

        Args:
            p0, p1, p2, p3 (np.ndarray): Four sequential atoms.

        Returns:
            float: Dihedral angle in degrees in range [-180, 180].
        """
        b0 = p1 - p0
        b1 = p2 - p1
        b2 = p3 - p2
        b1 /= np.linalg.norm(b1)
        v = np.cross(b0, b1)
        w = np.cross(b2, b1)
        x = np.dot(v, w)
        y = np.dot(np.cross(v, w), b1)
        angle = degrees(np.arctan2(y, x))

        if angle > 180:
            angle -= 360
        if angle < -180:
            angle += 360
        return angle
    
    def _rotation_matrix(
        self,
        axis: np.ndarray,
        theta_deg: float
    ) -> np.ndarray:
        """
        Compute 3D rotation matrix about an axis using Rodrigues' formula.

        Args:
            axis (np.ndarray): Rotation axis (3D vector).
            theta_deg (float): Rotation angle in degrees.

        Returns:
            np.ndarray: 3x3 rotation matrix.
        """
        theta = np.deg2rad(theta_deg)
        axis = axis / np.linalg.norm(axis)
        a = np.cos(theta / 2)
        b, c, d = -axis * np.sin(theta / 2)
        return np.array([
            [a*a + b*b - c*c - d*d, 2*(b*c - a*d),     2*(b*d + a*c)],
            [2*(b*c + a*d),     a*a + c*c - b*b - d*d, 2*(c*d - a*b)],
            [2*(b*d - a*c),     2*(c*d + a*b),     a*a + d*d - b*b - c*c]
        ])
    
    def _change_glycan_serial_id(
        self,
        glycan_lines: List[str],
        glycan_idx: int,
        glycan_chain: str = "G",
    ) -> List[str]:
        """
        Re-label a glycan block with a unique glycan ID.

                NOTE: PDB chain ID supports only 1 character. We keep chain ID as "G"
                and store "Gn" in segID field (cols 73-76).

        Args:
            glycan_lines (List[str]): PDB lines of one glycan block.
            glycan_idx (int): 1-based index of glycan insertion.
            glycan_chain (str): One-letter chain ID for glycans (default "G").

        Returns:
            List[str]: Updated PDB lines for the glycan block.
        """
        tag = f"G{glycan_idx}"[-4:]
        new_lines = []

        for line in glycan_lines:
            if not (line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("TER")):
                new_lines.append(line)
                continue

            padded = line.rstrip("\n").ljust(80)
            chars = list(padded)
            chars[21:23] = list(tag)
            chars[72:76] = list("    ")
            new_lines.append("".join(chars) + "\n")

        return new_lines

    def _merge_structures(
        self,
        protein_structure,
        glycan_structure,
        glycan_chain: str,
        output_path: str,
        glycan_idx: int,
    ) -> None:
        """
        Merge glycan chain into protein structure and save to file.

        Args:
            protein_structure: Bio.PDB structure object (protein).
            glycan_structure: Bio.PDB structure object (glycan).
            glycan_chain (str): Chain ID for the glycan to merge.
            output_path (str): Path to write merged PDB file.
            glycan_idx (int): 1-based index of glycan insertion.
        """

        class GlycanSelect(Select):
            def accept_chain(self, chain):
                return chain.get_id() == glycan_chain
        io = PDBIO()
        io.set_structure(glycan_structure)
        io.save("temp_glycan.pdb", GlycanSelect())
        with open("./temp_glycan.pdb") as f:
            glycan_lines = f.readlines()
        glycan_lines = self._change_glycan_serial_id(
            glycan_lines=glycan_lines,
            glycan_idx=glycan_idx,
        )
        if os.path.exists(output_path):
            lines = open(output_path).readlines()
        else:
            with open(output_path, "w") as f:
                io.set_structure(protein_structure)
                io.save(f)
            lines = open(output_path).readlines()

        while lines and lines[-1].strip() == "END":
            lines = lines[:-1]

        lines.extend(glycan_lines)
        with open(output_path, "w") as f:
            for line in lines:
                if line.startswith("ATOM") and ("ROH" in line and "SYST" in line):
                    continue
                else:
                    f.write(line)
            f.write("END")
        os.system("rm -f ./temp_glycan.pdb")

    def move(
        self,
        protein_structure_file: str,
        glycan_structure_file: str,
        output_pdb: str,
        glycan_positions: dict = None,
        glycan_chain: str = "X",
        nag_res_id: int = 2,
    ) -> None:
        """
        Attach a glycan to an Asn residue in a protein structure with specified geometry.

        Args:
            protein_structure_file (str): Path to PDB file of the protein.
            glycan_structure_file (str): Path to PDB file of the glycan.
            output_pdb (str): Path to save the merged structure.
            glycan_positions (dict): Mapping of protein chain ID to Asn residue ID for glycan attachment (e.g., {"A": [2,5]}).
            glycan_chain (str, optional): Chain ID of the glycan. Default is "X".
            nag_res_id (int, optional): Residue ID of the first glycan residue. Default is 2.
        """

        prot, _ = StructureLoader.load_structure(structure_file=protein_structure_file)
        glyc, _ = StructureLoader.load_structure(structure_file=glycan_structure_file)
        io = PDBIO()
        io.set_structure(prot)
        io.save(output_pdb, Select())

        if glycan_positions is not None:
            glycan_idx = 1
            for protein_chain, asn_res_ids in glycan_positions.items():
                for asn_res_id in asn_res_ids:
                    ND2, CG, CB = self._get_atom_coords(prot, protein_chain, asn_res_id, ["ND2", "CG", "CB"])
                    C1, C2, O5 = self._get_atom_coords(glyc, glycan_chain, nag_res_id, ["C1", "C2", "O5"])

                    C1_target = self._place_atom(ND2, self.bond_length, self.angle_c1, self.dihedral_c1 + 180, ND2, CG, CB)
                    C2_target = self._place_atom(C1_target, self.bond_length, self.angle_c2, self.dihedral_c2, C1_target, ND2, CG)
                    O5_target = self._place_atom(C1_target, self.bond_length, self.angle_o5, self.dihedral_o5, C1_target, ND2, CG)

                    R_local = self._build_frame(C1, C2, O5)
                    R_target = self._build_frame(C1_target, C2_target, O5_target)

                    R = R_target @ R_local.T
                    t = C1_target - R @ C1

                    for res in glyc[0][glycan_chain]:
                        for atom in res:
                            atom.coord = R @ atom.coord + t

                    C1_new, O5_new = self._get_atom_coords(glyc, glycan_chain, nag_res_id, ["C1", "O5"])
                    ND2_new = ND2
                    CG_new = CG

                    current_dihedral = self._dihedral_angle(O5_new, C1_new, ND2_new, CG_new)
                    dihedral_diff = self.dihedral_o5 + 180 - current_dihedral
                    if dihedral_diff > 180:
                        dihedral_diff -= 360
                    elif dihedral_diff < -180:
                        dihedral_diff += 360

                    axis = ND2_new - C1_new
                    axis /= np.linalg.norm(axis)
                    R_adjust = self._rotation_matrix(axis, dihedral_diff)

                    for res in glyc[0][glycan_chain]:
                        for atom in res:
                            atom.coord = C1_new + R_adjust @ (atom.coord - C1_new)

                    self._merge_structures(prot, glyc, glycan_chain, output_pdb, glycan_idx)
                    glycan_idx += 1
    
class InteractionCheck:
    def __init__(self) -> None:
        """
        Initialize the InteractionCheck object.
        """
        pass

    def _extract_cb(
        self,
        structure: Structure,
        chain_id: str = None,
    ) -> Tuple[np.ndarray, List[int]]:
        """
        Extract Cβ (or Cα for Gly) coordinates and their residue IDs.

        Args:
            structure (Structure): Biopython Structure object parsed from a PDB/mmCIF file.

        Returns:
            coords (np.ndarray): Cβ/Cα coordinates with shape (N, 3).
            res_ids (List[int]): Corresponding residue sequence numbers from the structure.
        """
        coords = []
        res_ids = []
        for model in structure:
            for chain in model:
                if chain_id is None or chain.get_id() == chain_id:
                    for residue in chain:
                        res = residue.resname
                        for atom in residue:
                            name = atom.get_fullname().strip(" ")
                            if (res == "GLY" and name == "CA") or (res != "GLY" and name == "CB"):
                                coords.append(atom.coord)
                                res_ids.append(f"{chain.get_id()}_{residue.get_id()[1]}")
                                break
            break
        return np.array(coords), res_ids

    def _compute_contact_map(
        self,
        coords: np.ndarray
    ) -> np.ndarray:
        """
        Compute pairwise Euclidean distances between coordinates.

        Args:
            coords (np.ndarray): An array of atom coordinates with shape (N, 3).

        Returns:
            np.ndarray: A distance matrix of shape (N, N), where each entry (i, j) is the distance between coords[i] and coords[j].
        """
        diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
        dists = np.linalg.norm(diff, axis=-1)
        return dists
    
    def get_inter_interaction_aa(
        self,
        structure_file: str,
        chain_id: str,
        dist_threshold: float = 6.0,
    ) -> Set[int]:
        """
        Get inter-chain interacting residues on one specific chain.

        Args:
            structure_file (str): Path to the input structure file (e.g., mmCIF or PDB).
            chain_id (str): Chain ID to extract Cβ/Cα coordinates.

        Returns:
            Set[int]: A set of residue sequence numbers that are in contact with the target residues.
        """
        structure, _ = StructureLoader.load_structure(structure_file)
        coords_cb, res_ids = self._extract_cb(structure=structure)
        contact = self._compute_contact_map(coords=coords_cb)
        mask = np.full(contact.shape, 1.0)
        groups = defaultdict(list)
        for i, res_id in enumerate(res_ids):
            chain, _ = res_id.split("_")
            groups[chain].append(i)
        
        for _, idx_list in groups.items():
            idx = np.array(idx_list)
            mask[np.ix_(idx, idx)] = 1e8
        
        contact = contact * mask
        np.fill_diagonal(contact, 1e8)

        interacting_idxs = np.unique(np.where(contact <= dist_threshold)[0])
        inter_interactions_aa = [int(res_ids[i].split("_")[1]) for i in interacting_idxs if res_ids[i].split("_")[0] == chain_id]
        inter_interactions_aa = set(inter_interactions_aa)

        return inter_interactions_aa

    def get_intra_interaction_aa(
        self,
        structure_file: str,
        chain_id: str,
        positions: List[int],
        dist_threshold: float = 6.0,
        is_self_included: bool = True,
        num_neighbors: int = 3
    ) -> Set[int]:
        """
        Extract possible interactions in the structure, given the target residues.

        Args:
            structure_file (str): Path to the input structure file (e.g., mmCIF or PDB).
            positions (List[int]): List of target residue sequence numbers (starting from 1).
            dist_threshold (float): CB distance threshold to define contact (in Å). Default is 6.0.
            is_self_included (bool): Whether to include themselves in the result, for detecting interacting residues around non-editable residues. Default is True.
            num_neighbors (int): Number of neighbor residues (in sequence) to include around each target residue.

        Returns:
            Set[int]: A set of residue sequence numbers that are in contact with the target residues.
        """
        structure, _ = StructureLoader.load_structure(structure_file)
        coords_cb, res_ids = self._extract_cb(structure=structure, chain_id=chain_id)
        contact = self._compute_contact_map(coords=coords_cb)
        
        pts = []
        for pt in positions:
            if isinstance(pt, int):
                pts.append(pt)
            elif isinstance(pt, str) and "-" in pt:
                start, end = int(pt.split("-")[0]), int(pt.split("-")[1])
                for p in range(start, end+1):
                    pts.append(p)

        interacting_aas = []
        for p in pts:
            if p not in res_ids:
                continue
            p_idx = res_ids.index(p)
            contact_p = contact[p_idx, :]
            interacting_idxs = np.where(contact_p <= dist_threshold)[0]
            interacting_res_ids = [res_ids[i] for i in interacting_idxs]
            interacting_aas.extend(interacting_res_ids)

        if is_self_included:
            return set(interacting_aas) | set(
                [p + i for p in pts for i in range(1, num_neighbors + 1) if (p + i) <= int(res_ids[-1].split("_")[1])]
            ) | set(
                [p - i for p in pts for i in range(1, num_neighbors + 1) if (p - i) > 0]
            )
        else:
            return set(interacting_aas) - set(pts) - set(
                [p + i for p in pts for i in range(1, num_neighbors + 1) if (p + i) <= int(res_ids[-1].split("_")[1])]
            ) - set(
                [p - i for p in pts for i in range(1, num_neighbors + 1) if (p - i) > 0]
            )

class ClashCheck:
    def __init__(self) -> None:
        """
        Initialize the ClashCheck class with van der Waals radii for common elements.
        These values are used to compute steric clash thresholds between atoms.
        """
        self.vdw = {"C": 1.82, "N": 1.64, "O": 1.44, "S": 1.77}

    def _extract_all_coords(
        self, 
        structure: Structure
    ) -> Tuple[np.ndarray, List[str], List[str], List[str]]:
        """
        Extract coordinates, element types, and chain IDs of all non-hydrogen atoms in the structure.

        Args:
            structure (Structure): Bio.PDB Structure object representing the macromolecule.

        Returns:
            Tuple containing:
                - coords (np.ndarray): Atom coordinates of shape (N, 3).
                - elems (List[str]): Element symbols corresponding to each atom.
                - chain_ids (List[str]): Chain ID for each atom.
                - res_ids (List[str]): Residue ID for each atom.
        """
        coords = []
        elems = []
        chain_ids = []
        res_ids = []

        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        elem = atom.element
                        if elem != "H":
                            coords.append(atom.get_coord())
                            elems.append(elem)
                            chain_ids.append(chain.id)
                            res_ids.append(f"{chain.id}_{residue.id[1]}")
            break

        return np.array(coords), elems, chain_ids, res_ids

    def _compute_contact_map(
        self, 
        coords: np.ndarray, 
        elems: List[str], 
        chain_ids: List[str],
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute pairwise Euclidean distances between atoms and corresponding steric cutoff thresholds.

        Args:
            coords (np.ndarray): Atom coordinates of shape (N, 3).
            elems (List[str]): Element types of each atom.
            chain_ids (List[str]): Chain ID of each atom.

        Returns:
            Tuple containing:
                - dists (np.ndarray): Pairwise distance matrix (N, N) with intra-chain distances set to inf.
                - steric (np.ndarray): Steric cutoff matrix (N, N) based on van der Waals radii.
        """
        diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
        dists = np.linalg.norm(diff, axis=-1)

        n = len(elems)
        steric = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                steric_val = self.vdw.get(elems[i], 1.5) + self.vdw.get(elems[j], 1.5) - 0.4
                steric[i][j] = steric_val
                steric[j][i] = steric_val

        for i in range(n):
            for j in range(i, n):
                if chain_ids[i] == chain_ids[j]:
                    dists[i][j] = np.inf
                    dists[j][i] = np.inf

        return dists, steric

    def has_clash(
        self,
        chain_id: str,
        secondary_structure: dict,
        structure_file: str,
    ) -> Tuple:
        """
        Determine whether there are any steric clashes between atoms across different chains.

        Args:
            structure_file (str): Path to the input structure file (e.g., mmCIF or PDB).

        Returns:
            bool: True if no inter-chain clashes are found, False otherwise.
        """
        structure, _ = StructureLoader.load_structure(structure_file)

        coords, elems, chain_ids, res_ids = self._extract_all_coords(structure)
        contact_map, steric_cutoff = self._compute_contact_map(coords, elems, chain_ids)

        clash = (contact_map < steric_cutoff).astype(np.uint8)
        clash_res_index = [item for item in [[res_ids[clash_pair[0]], res_ids[clash_pair[1]]] for clash_pair in np.argwhere(clash == 1) if clash_pair[0] < clash_pair[1]]]

        if clash_res_index:
            clash_res_index = [item for item in clash_res_index if item[1] != "X_2" and not item[1].endswith("_1") and SS_TAG[secondary_structure[(chain_id, int(item[0].split('_')[-1]))]][-1] != "loop"]

            return clash_res_index if clash_res_index else None

class BordaCount:
    def __init__(
        self,
        conservation_weight: float = 1.0,
        coupling_weight: float = 0.2,
        sasa_weight: float = 1.0,
        sasa_next1_weight: float = 0.7,
        sasa_next2_weight: float = 0.5,
        sasa_around_weight: float = 1.0,
        sasa_next_weight: float = 0.6,
        ddG_weight: float = 0.5,
        dTm_weight: float = 0.5,
        ddG_S_weight: float = 0.3,
        dTm_S_weight: float = 0.3,
        ddG_T_weight: float = 0.3,
        dTm_T_weight: float = 0.3,
        mut_score_weight: float = 0.3,
        mut_score_S_weight: float = 0.2,
        mut_score_T_weight: float = 0.2,
    ) -> None:
        """
        Initialize the BordaCount class with weights for different features.
        """

        self.weights = {
            "ConservationScore": conservation_weight,
            "CouplingScore": coupling_weight,
            "SASA_i": sasa_weight,
            "SASA_i+1": sasa_next1_weight,
            "SASA_i+2": sasa_next2_weight,
            "SASA_(i-1:i+1)": sasa_around_weight,
            "SASA_(i:i+2)": sasa_next_weight,
            "ddG": ddG_weight,
            "dTm": dTm_weight,
            "ddG_NXS": ddG_S_weight,
            "dTm_NXS": dTm_S_weight,
            "ddG_NXT": ddG_T_weight,
            "dTm_NXT": dTm_T_weight,
            "MutScore": mut_score_weight,
            "MutScore_NXS": mut_score_S_weight,
            "MutScore_NXT": mut_score_T_weight,
        }
        self.higher_is_better = {
            "ConservationScore": False,
            "CouplingScore": False,
            "SASA_i": True,
            "SASA_i+1": True,
            "SASA_i+2": True,
            "SASA_(i-1:i+1)": True,
            "SASA_(i:i+2)": True,
            "ddG": True,
            "dTm": True,
            "ddG_NXS": True,
            "dTm_NXS": True,
            "ddG_NXT": True,
            "dTm_NXT": True,
            "MutScore": True,
            "MutScore_NXS": True,
            "MutScore_NXT": True,
        }

    def compute_score(self, dataframe: pd.DataFrame) -> pd.DataFrame:
        """
        Compute the Borda score for each row in the dataframe.

        Args:
            dataframe (pd.DataFrame): Input dataframe with features.

        Returns:
            pd.DataFrame: Dataframe with added Borda score column.
        """
        df = dataframe.copy()
        scores = np.zeros(len(df))

        for col, weight in self.weights.items():
            if col not in df.columns:
                continue

            ascending = not self.higher_is_better[col]
            ranks = df[col].rank(method="min", ascending=ascending)
            max_rank = ranks.max()
            borda_scores = max_rank - ranks + 1
            scores += borda_scores * weight

        min_score = scores.min()
        max_score = scores.max()
        norm_scores = (scores - min_score) / (max_score - min_score + 1e-8)

        df.insert(0, "Borda_score", norm_scores)
        df["Borda_score"] = df["Borda_score"].round(3)
        df.sort_values(by="Borda_score", ascending=False, inplace=True)
        return df

class prefilter_report:
    def __init__(self, input_fasta_file: str, query_sequence: str, editable_regions: set, df_file: str, output_html: str):
        self.input_fasta_file = input_fasta_file
        self.query_sequence = query_sequence
        self.editable_regions = editable_regions
        self.df_file = df_file
        self.output_html = output_html
    
    def plot_heatmap(
        self,
        out_file: str,
    ):

        df = pd.read_csv(self.df_file)
        higher_is_better = {
            "ConservationScore": False,
            "CouplingScore": False,
            "SASA_i": True,
            "SASA_i+1": True,
            "SASA_i+2": True,
            "SASA_(i-1:i+1)": True,
            "SASA_(i:i+2)": True,
            "ddG": True,
            "dTm": True,
            "ddG_NXS": True,
            "dTm_NXS": True,
            "ddG_NXT": True,
            "dTm_NXT": True,
            "MutScore": True,
            "MutScore_NXS": True,
            "MutScore_NXT": True,
        }
        for c, h in higher_is_better.items():
            if c in df.columns and not h:
                df[c] = -df[c]
        for col in ["Site", "SS", "Clash"]:
            if col in df.columns:
                del df[col]
        df.set_index("Mutation", inplace=True)

        scaler = MinMaxScaler()
        df_normalized = pd.DataFrame(
            scaler.fit_transform(df.values),
            columns=df.columns,
            index=df.index
        )

        plt.figure(figsize=(8.5, 7))
        three_color_cmap = mcolors.LinearSegmentedColormap.from_list(
            "three_color_cmap",
            ["#91C8F6", "#FDF6ED", "#F6A6A1"],
            N=256
        )
        sns.heatmap(
            df_normalized,
            cmap=three_color_cmap,
            cbar=True,
            annot_kws={"size": 4},
            vmin=0, vmax=1
        )

        plt.title("Single Point Mutations", fontsize=12)
        plt.rcParams["font.size"] = 12
        plt.rcParams["axes.linewidth"] = 0.5
        plt.rcParams["xtick.major.width"] = 0.5
        plt.rcParams["ytick.major.width"] = 0.5
        plt.rcParams["xtick.major.size"] = 3
        plt.rcParams["ytick.major.size"] = 3
        plt.rcParams["xtick.direction"] = "out"
        plt.rcParams["ytick.direction"] = "out"
        plt.rcParams["axes.labelpad"] = 5
        plt.rcParams["xtick.labeltop"] = False
        plt.rcParams["xtick.labelbottom"] = True
        plt.rcParams["ytick.labelleft"] = True
        plt.rcParams["ytick.labelright"] = False
        plt.rcParams["xtick.major.pad"] = 3.5
        plt.rcParams["ytick.major.pad"] = 3.5
        plt.rcParams["axes.grid"] = False
        plt.savefig(out_file, dpi=300, bbox_inches="tight")

    def _compress_ranges(self, positions: list) -> str:
        """Compress sorted integer positions into range strings like '1-5, 8, 10-15'."""
        if not positions:
            return ""
        ranges = []
        start = positions[0]
        end = positions[0]
        for p in positions[1:]:
            if p == end + 1:
                end = p
            else:
                ranges.append(f"{start}-{end}" if start != end else str(start))
                start = end = p
        ranges.append(f"{start}-{end}" if start != end else str(start))
        return ", ".join(ranges)

    def _generate_heatmap_base64(self, df_file: str) -> str:
        """Generate heatmap as base64-encoded PNG string from CSV results."""
        import base64
        from io import BytesIO

        df = pd.read_csv(df_file)
        higher_is_better = {
            "ConservationScore": False, "CouplingScore": False,
            "SASA_i": True, "SASA_i+1": True, "SASA_i+2": True,
            "SASA_(i-1:i+1)": True, "SASA_(i:i+2)": True,
            "ddG": True, "dTm": True,
            "ddG_NXS": True, "dTm_NXS": True, "ddG_NXT": True, "dTm_NXT": True,
            "MutScore": True, "MutScore_NXS": True, "MutScore_NXT": True,
        }
        for c, h in higher_is_better.items():
            if c in df.columns and not h:
                df[c] = -df[c]
        for col in ["Site", "SS", "Clash"]:
            if col in df.columns:
                del df[col]
        df.set_index("Mutation", inplace=True)

        scaler = MinMaxScaler()
        df_normalized = pd.DataFrame(
            scaler.fit_transform(df.values),
            columns=df.columns,
            index=df.index,
        )

        fig, ax = plt.subplots(figsize=(8.5, max(3, len(df) * 0.35 + 1.5)))
        three_color_cmap = mcolors.LinearSegmentedColormap.from_list(
            "three_color_cmap", ["#91C8F6", "#FDF6ED", "#F6A6A1"], N=256,
        )
        sns.heatmap(
            df_normalized, cmap=three_color_cmap, cbar=True,
            annot_kws={"size": 4}, vmin=0, vmax=1, ax=ax,
        )
        ax.set_title("Single Point Mutations", fontsize=12)

        buf = BytesIO()
        fig.savefig(buf, format="png", dpi=200, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        return base64.b64encode(buf.read()).decode("utf-8")

    def generate_prefilter_report(
        self,
    ):
        """
        Generate a self-contained HTML report from prefilter results.
        """
        
        df = pd.read_csv(self.df_file)

        # --- heatmap ---
        heatmap_b64 = self._generate_heatmap_base64(self.df_file)

        # --- sequence with highlights ---
        line_width = 50
        seq_lines = []
        for start in range(0, len(self.query_sequence), line_width):
            end = min(start + line_width, len(self.query_sequence))
            pos_label = f'<span class="pos-num">{start + 1:>5}</span> '
            chars = []
            for i in range(start, end):
                pos = i + 1
                aa = self.query_sequence[i]
                if pos in self.editable_regions:
                    chars.append(
                        f'<span class="aa editable" title="Pos {pos} &#10004; editable">{aa}</span>'
                    )
                else:
                    chars.append(
                        f'<span class="aa non-edit" title="Pos {pos} &#10008; non-editable">{aa}</span>'
                    )
                if (i - start + 1) % 10 == 0 and i + 1 != end:
                    chars.append('<span class="gap"> </span>')
            end_label = f' <span class="pos-num">{end}</span>'
            seq_lines.append(pos_label + "".join(chars) + end_label)
        sequence_html = "\n".join(seq_lines)

        # --- editable regions summary ---
        editable_sorted = sorted(self.editable_regions)
        non_editable = sorted(set(range(1, len(self.query_sequence) + 1)) - self.editable_regions)
        editable_str = self._compress_ranges(editable_sorted)
        non_editable_str = self._compress_ranges(non_editable)

        # --- results table ---
        table_df = df.copy()
        table_df["Clash"] = table_df["Clash"].apply(
            lambda x: "None" if pd.isna(x) or x == "None" or x == "[]" else str(x)
        )
        table_html = table_df.to_html(
            index=False, classes="data-table", border=0,
            float_format=lambda x: f"{x:.3f}",
        )

        # --- assemble HTML ---
        html = f"""<!DOCTYPE html>
    <html lang="en">
    <head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SugarSwitch Prefilter Report</title>
    <style>
    :root {{
        --bg: #f5f7fa;
        --card: #ffffff;
        --accent: #4a76a8;
        --accent-light: #e8f0fe;
        --edit-bg: #d4edda;
        --edit-border: #28a745;
        --nonedit-bg: #f8f9fa;
        --text: #2c3e50;
        --muted: #6c757d;
        --border: #dee2e6;
        --shadow: 0 2px 8px rgba(0,0,0,0.08);
    }}
    * {{ margin: 0; padding: 0; box-sizing: border-box; }}
    body {{
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
        background: var(--bg);
        color: var(--text);
        line-height: 1.6;
        padding: 20px;
    }}
    .container {{ max-width: 1200px; margin: 0 auto; }}
    header {{
        background: linear-gradient(135deg, #4a76a8 0%, #6a9fd8 100%);
        color: #fff;
        padding: 28px 36px;
        border-radius: 12px;
        margin-bottom: 24px;
        box-shadow: var(--shadow);
    }}
    header h1 {{ font-size: 24px; font-weight: 700; }}
    header p {{ opacity: 0.9; margin-top: 4px; font-size: 14px; }}
    .card {{
        background: var(--card);
        border-radius: 10px;
        padding: 24px 28px;
        margin-bottom: 20px;
        box-shadow: var(--shadow);
    }}
    .card h2 {{
        font-size: 17px;
        font-weight: 600;
        color: var(--accent);
        margin-bottom: 14px;
        padding-bottom: 8px;
        border-bottom: 2px solid var(--accent-light);
    }}
    .info-grid {{
        display: grid;
        grid-template-columns: 140px 1fr;
        gap: 8px 16px;
        font-size: 14px;
    }}
    .info-grid .label {{ font-weight: 600; color: var(--muted); }}
    .info-grid .value {{ word-break: break-all; }}

    /* sequence viewer */
    .seq-viewer {{
        font-family: "SF Mono", "Fira Code", "Cascadia Code", Consolas, monospace;
        font-size: 13px;
        line-height: 2;
        white-space: pre;
        overflow-x: auto;
        padding: 16px;
        background: #fafbfc;
        border: 1px solid var(--border);
        border-radius: 8px;
    }}
    .pos-num {{ color: var(--muted); font-size: 11px; user-select: none; }}
    .aa {{
        display: inline-block;
        width: 1.1em;
        text-align: center;
        border-radius: 3px;
        cursor: default;
        transition: transform 0.1s;
    }}
    .aa:hover {{ transform: scale(1.3); z-index: 1; position: relative; }}
    .editable {{
        background: var(--edit-bg);
        color: #155724;
        font-weight: 600;
        border-bottom: 2px solid var(--edit-border);
    }}
    .non-edit {{
        background: var(--nonedit-bg);
        color: #888;
    }}
    .gap {{ width: 0.5em; display: inline-block; }}

    .legend {{
        display: flex;
        gap: 20px;
        margin-top: 12px;
        font-size: 13px;
    }}
    .legend-item {{ display: flex; align-items: center; gap: 6px; }}
    .legend-swatch {{
        display: inline-block;
        width: 16px; height: 16px;
        border-radius: 3px;
    }}
    .swatch-edit {{ background: var(--edit-bg); border: 1px solid var(--edit-border); }}
    .swatch-nonedit {{ background: var(--nonedit-bg); border: 1px solid var(--border); }}

    /* regions */
    .region-box {{
        font-family: "SF Mono", Consolas, monospace;
        font-size: 13px;
        padding: 10px 14px;
        border-radius: 6px;
        margin-bottom: 10px;
        word-break: break-word;
    }}
    .region-edit {{ background: #d4edda; border-left: 4px solid #28a745; }}
    .region-nonedit {{ background: #fff3cd; border-left: 4px solid #ffc107; }}
    .region-label {{ font-weight: 600; margin-bottom: 2px; font-size: 13px; }}

    .stats {{
        display: flex;
        gap: 16px;
        margin-bottom: 14px;
        flex-wrap: wrap;
    }}
    .stat {{
        background: var(--accent-light);
        padding: 8px 16px;
        border-radius: 8px;
        font-size: 14px;
    }}
    .stat strong {{ color: var(--accent); }}

    /* table */
    .table-wrap {{ overflow-x: auto; }}
    .data-table {{
        width: 100%;
        border-collapse: collapse;
        font-size: 13px;
        white-space: nowrap;
    }}
    .data-table th {{
        background: var(--accent);
        color: #fff;
        padding: 10px 12px;
        text-align: left;
        font-weight: 600;
        position: sticky;
        top: 0;
    }}
    .data-table td {{
        padding: 8px 12px;
        border-bottom: 1px solid var(--border);
    }}
    .data-table tr:hover {{ background: var(--accent-light); }}
    .data-table tr:nth-child(even) {{ background: #f8f9fa; }}
    .data-table tr:nth-child(even):hover {{ background: var(--accent-light); }}

    /* heatmap */
    .heatmap-img {{
        max-width: 100%;
        height: auto;
        border-radius: 6px;
        border: 1px solid var(--border);
    }}
    footer {{
        text-align: center;
        color: var(--muted);
        font-size: 12px;
        margin-top: 24px;
        padding: 16px;
    }}
    </style>
    </head>
    <body>
    <div class="container">

    <header>
        <h1>SugarSwitch &mdash; Prefilter Report</h1>
        <p>N-linked glycosylation site screening results</p>
    </header>

    <!-- Input Info -->
    <div class="card">
        <h2>Input Information</h2>
        <div class="info-grid">
            <div class="label">FASTA File</div>
            <div class="value">{self.input_fasta_file}</div>
            <div class="label">Sequence Length</div>
            <div class="value">{len(self.query_sequence)} residues</div>
            <div class="label">Editable Sites</div>
            <div class="value">{len(editable_sorted)} / {len(self.query_sequence)}</div>
            <div class="label">Candidates</div>
            <div class="value">{len(df)} mutations evaluated</div>
        </div>
    </div>

    <!-- Query Sequence -->
    <div class="card">
        <h2>Query Sequence</h2>
        <div class="seq-viewer">
{sequence_html}
        </div>
        <div class="legend">
            <div class="legend-item"><span class="legend-swatch swatch-edit"></span> Editable</div>
            <div class="legend-item"><span class="legend-swatch swatch-nonedit"></span> Non-editable (functional / conserved / co-evolving / low-SASA)</div>
        </div>
    </div>

    <!-- Editable Regions -->
    <div class="card">
        <h2>Editable Regions</h2>
        <div class="stats">
            <div class="stat"><strong>{len(editable_sorted)}</strong> editable positions</div>
            <div class="stat"><strong>{len(non_editable)}</strong> filtered out</div>
        </div>
        <div class="region-box region-edit">
            <div class="region-label">Editable positions:</div>
            {editable_str}
        </div>
        <div class="region-box region-nonedit">
            <div class="region-label">Non-editable positions:</div>
            {non_editable_str}
        </div>
    </div>

    <!-- Results Table -->
    <div class="card">
        <h2>Ranked Mutation Results</h2>
        <div class="table-wrap">
            {table_html}
        </div>
    </div>

    <!-- Heatmap -->
    <div class="card">
        <h2>Heatmap</h2>
        <img class="heatmap-img" src="data:image/png;base64,{heatmap_b64}" alt="Single Point Mutations Heatmap">
    </div>

    <footer>
        Generated by SugarSwitch &bull; Prefilter Pipeline
    </footer>

    </div>
    </body>
    </html>"""

        with open(self.output_html, "w") as f:
            f.write(html)


class designer_report:
    def __init__(
        self,
        input_fasta_file: str,
        wt_seq: str,
        asn_sites: dict,
        designed_seqs: List[str],
        pos_probs_list: List[dict],
        output_html: str,
    ):
        """
        Generate an HTML report for halludesign_esm results.

        Args:
            input_fasta_file: Path to input FASTA file.
            wt_seq: Wild-type sequence (with X replaced by NP).
            asn_sites: Dict mapping chain_id -> list of Asn positions, e.g. {"A": [10, 50]}.
            designed_seqs: List of designed sequences.
            pos_probs_list: List of dicts, each mapping position -> probability for a designed sequence.
            output_html: Output HTML file path.
        """
        self.input_fasta_file = input_fasta_file
        self.wt_seq = wt_seq
        self.asn_sites = asn_sites
        self.designed_seqs = designed_seqs
        self.pos_probs_list = pos_probs_list
        self.output_html = output_html

    def _format_seq_html(self, seq: str, pos_probs: dict = None, label: str = "Sequence", highlight_positions: set = None):
        """Render a sequence as HTML with optional position highlighting and probability tooltips."""
        line_width = 50
        seq_lines = []
        for start in range(0, len(seq), line_width):
            end = min(start + line_width, len(seq))
            pos_label = f'<span class="pos-num">{start + 1:>5}</span> '
            chars = []
            for i in range(start, end):
                pos = i + 1
                aa = seq[i]
                if pos_probs and pos in pos_probs:
                    prob = pos_probs[pos]
                    prob_pct = f"{prob * 100:.1f}%"
                    css_class = "aa pos-hit" if prob >= 0.5 else "aa pos-miss"
                    chars.append(
                        f'<span class="{css_class}" data-pos="{pos}" data-prob="{prob:.4f}" data-pct="{prob_pct}">{aa}</span>'
                    )
                elif highlight_positions and pos in highlight_positions:
                    chars.append(
                        f'<span class="aa asn-site" data-pos="{pos}" data-info="Possible glycosylation site">{aa}</span>'
                    )
                else:
                    chars.append(
                        f'<span class="aa aa-default" data-pos="{pos}">{aa}</span>'
                    )
                if (i - start + 1) % 10 == 0 and i + 1 != end:
                    chars.append('<span class="gap"> </span>')
            end_label = f' <span class="pos-num">{end}</span>'
            seq_lines.append(pos_label + "".join(chars) + end_label)
        return "\n".join(seq_lines)

    def generate_designer_report(self):
        """Generate a self-contained HTML report for halludesign_esm results."""

        # Collect all Asn positions across chains
        all_asn_positions = set()
        for chain_id, positions in self.asn_sites.items():
            all_asn_positions.update(positions)
        asn_sites_str = "; ".join(
            f"Chain {chain}: {', '.join(str(p) for p in sorted(positions))}"
            for chain, positions in self.asn_sites.items()
        )

        # WT sequence HTML
        wt_seq_html = self._format_seq_html(self.wt_seq, highlight_positions=all_asn_positions, label="WT Sequence")

        # Designed sequences HTML
        designed_cards = []
        for idx, (d_seq, d_probs) in enumerate(zip(self.designed_seqs, self.pos_probs_list)):
            seq_html = self._format_seq_html(d_seq, pos_probs=d_probs, label=f"Design {idx + 1}")

            # Mutation summary
            mutations = []
            for i, (wt_aa, d_aa) in enumerate(zip(self.wt_seq, d_seq)):
                if wt_aa != d_aa:
                    mutations.append(f"{wt_aa}{i+1}{d_aa}")
            mut_str = ", ".join(mutations) if mutations else "None"

            # Prob summary table
            prob_rows = ""
            for pos in sorted(d_probs.keys()):
                prob = d_probs[pos]
                motif = d_seq[pos-1:pos+2] if pos-1+3 <= len(d_seq) else d_seq[pos-1:]
                pred = "Positive" if prob >= 0.5 else "Negative"
                prob_rows += f"<tr><td>{pos}</td><td>{motif}</td><td>{prob:.4f}</td><td>{pred}</td></tr>\n"

            designed_cards.append(f"""
    <div class="card">
        <h2>Design {idx + 1}</h2>
        <div class="info-grid" style="margin-bottom:14px;">
            <div class="label">Mutations</div>
            <div class="value">{mut_str}</div>
        </div>
        <div class="seq-viewer">
{seq_html}
        </div>
        <div class="legend">
            <div class="legend-item"><span class="legend-swatch swatch-pos-hit"></span> Candidate pos (prob &ge; 0.5)</div>
            <div class="legend-item"><span class="legend-swatch swatch-pos-miss"></span> Candidate pos (prob &lt; 0.5)</div>
            <div class="legend-item"><span class="legend-swatch swatch-default"></span> Other residues</div>
        </div>
        <h3 style="margin-top:16px;font-size:14px;color:#4a76a8;">Glycosylation Predictions</h3>
        <div class="table-wrap" style="margin-top:8px;">
            <table class="data-table">
                <thead><tr><th>Position</th><th>Motif</th><th>Probability</th><th>Prediction</th></tr></thead>
                <tbody>{prob_rows}</tbody>
            </table>
        </div>
    </div>""")

        designed_html = "\n".join(designed_cards)

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>SugarSwitch Designer Report</title>
<style>
:root {{
    --bg: #f5f7fa;
    --card: #ffffff;
    --accent: #4a76a8;
    --accent-light: #e8f0fe;
    --text: #2c3e50;
    --muted: #6c757d;
    --border: #dee2e6;
    --shadow: 0 2px 8px rgba(0,0,0,0.08);
}}
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
    background: var(--bg);
    color: var(--text);
    line-height: 1.6;
    padding: 20px;
}}
.container {{ max-width: 1200px; margin: 0 auto; }}
header {{
    background: linear-gradient(135deg, #4a76a8 0%, #6a9fd8 100%);
    color: #fff;
    padding: 28px 36px;
    border-radius: 12px;
    margin-bottom: 24px;
    box-shadow: var(--shadow);
}}
header h1 {{ font-size: 24px; font-weight: 700; }}
header p {{ opacity: 0.9; margin-top: 4px; font-size: 14px; }}
.card {{
    background: var(--card);
    border-radius: 10px;
    padding: 24px 28px;
    margin-bottom: 20px;
    box-shadow: var(--shadow);
}}
.card h2 {{
    font-size: 17px;
    font-weight: 600;
    color: var(--accent);
    margin-bottom: 14px;
    padding-bottom: 8px;
    border-bottom: 2px solid var(--accent-light);
}}
.info-grid {{
    display: grid;
    grid-template-columns: 160px 1fr;
    gap: 8px 16px;
    font-size: 14px;
}}
.info-grid .label {{ font-weight: 600; color: var(--muted); }}
.info-grid .value {{ word-break: break-all; }}

.seq-viewer {{
    font-family: "SF Mono", "Fira Code", "Cascadia Code", Consolas, monospace;
    font-size: 13px;
    line-height: 2;
    white-space: pre;
    overflow-x: auto;
    padding: 16px;
    background: #fafbfc;
    border: 1px solid var(--border);
    border-radius: 8px;
}}
.pos-num {{ color: var(--muted); font-size: 11px; user-select: none; }}
.aa {{
    display: inline-block;
    width: 1.1em;
    text-align: center;
    border-radius: 3px;
    cursor: default;
    transition: transform 0.1s;
}}
.aa:hover {{ transform: scale(1.3); z-index: 1; position: relative; }}
.aa-default {{ background: #f8f9fa; color: #888; }}
.asn-site {{
    background: #fff3cd;
    color: #856404;
    font-weight: 600;
    border-bottom: 2px solid #ffc107;
}}
.pos-hit {{
    background: #d4edda;
    color: #155724;
    font-weight: 600;
    border-bottom: 2px solid #28a745;
}}
.pos-miss {{
    background: #f8d7da;
    color: #721c24;
    font-weight: 600;
    border-bottom: 2px solid #dc3545;
}}
.gap {{ width: 0.5em; display: inline-block; }}

.legend {{
    display: flex;
    gap: 20px;
    margin-top: 12px;
    font-size: 13px;
}}
.legend-item {{ display: flex; align-items: center; gap: 6px; }}
.legend-swatch {{
    display: inline-block;
    width: 16px; height: 16px;
    border-radius: 3px;
}}
.swatch-asn {{ background: #fff3cd; border: 1px solid #ffc107; }}
.swatch-pos-hit {{ background: #d4edda; border: 1px solid #28a745; }}
.swatch-pos-miss {{ background: #f8d7da; border: 1px solid #dc3545; }}
.swatch-default {{ background: #f8f9fa; border: 1px solid var(--border); }}

.table-wrap {{ overflow-x: auto; }}
.data-table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 13px;
    white-space: nowrap;
}}
.data-table th {{
    background: var(--accent);
    color: #fff;
    padding: 10px 12px;
    text-align: left;
    font-weight: 600;
    position: sticky;
    top: 0;
}}
.data-table td {{
    padding: 8px 12px;
    border-bottom: 1px solid var(--border);
}}
.data-table tr:hover {{ background: var(--accent-light); }}
.data-table tr:nth-child(even) {{ background: #f8f9fa; }}
.data-table tr:nth-child(even):hover {{ background: var(--accent-light); }}
footer {{
    text-align: center;
    color: var(--muted);
    font-size: 12px;
    margin-top: 24px;
    padding: 16px;
}}
#seq-tooltip {{
    position: fixed;
    pointer-events: none;
    background: #2c3e50;
    color: #fff;
    padding: 6px 12px;
    border-radius: 6px;
    font-size: 13px;
    font-family: "SF Mono", Consolas, monospace;
    white-space: nowrap;
    z-index: 9999;
    box-shadow: 0 4px 12px rgba(0,0,0,0.25);
    opacity: 0;
    transition: opacity 0.15s;
}}
#seq-tooltip.visible {{ opacity: 1; }}
</style>
</head>
<body>
<div id="seq-tooltip"></div>
<div class="container">

<header>
    <h1>SugarSwitch &mdash; Designer Report</h1>
    <p>ESM-guided glycosylation sequence design results</p>
</header>

<!-- Input Info -->
<div class="card">
    <h2>Input Information</h2>
    <div class="info-grid">
        <div class="label">FASTA File</div>
        <div class="value">{self.input_fasta_file}</div>
        <div class="label">WT Sequence Length</div>
        <div class="value">{len(self.wt_seq)} residues</div>
        <div class="label">Asn Sites</div>
        <div class="value">{asn_sites_str}</div>
        <div class="label">Number of Designs</div>
        <div class="value">{len(self.designed_seqs)}</div>
    </div>
</div>

<!-- WT Sequence -->
<div class="card">
    <h2>Wild-Type Sequence</h2>
    <div class="seq-viewer">
{wt_seq_html}
    </div>
    <div class="legend">
        <div class="legend-item"><span class="legend-swatch swatch-asn"></span> Asn glycosylation site</div>
        <div class="legend-item"><span class="legend-swatch swatch-default"></span> Other residues</div>
    </div>
</div>

<!-- Designed Sequences -->
{designed_html}

<footer>
    Generated by SugarSwitch &bull; Hallucination-designer Pipeline
</footer>

</div>
<script>
(function() {{
    const tip = document.getElementById('seq-tooltip');
    document.addEventListener('mouseover', function(e) {{
        const el = e.target.closest('.aa');
        if (!el) return;
        const pos = el.dataset.pos;
        if (!pos) return;
        let text = 'Pos ' + pos;
        if (el.dataset.prob) {{
            text += ' | prob: ' + el.dataset.prob + ' (' + el.dataset.pct + ')';
        }} else if (el.dataset.info) {{
            text += ' | ' + el.dataset.info;
        }}
        tip.textContent = text;
        tip.classList.add('visible');
    }});
    document.addEventListener('mousemove', function(e) {{
        if (tip.classList.contains('visible')) {{
            tip.style.left = (e.clientX + 12) + 'px';
            tip.style.top = (e.clientY - 32) + 'px';
        }}
    }});
    document.addEventListener('mouseout', function(e) {{
        if (e.target.closest('.aa')) {{
            tip.classList.remove('visible');
        }}
    }});
}})();
</script>
</body>
</html>"""

        with open(self.output_html, "w") as f:
            f.write(html)
