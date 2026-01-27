'''
-------------------- 骨架设计 --------------------
找到位点--> DONE
定位一级序列和三级结构位置，计算一级序列和三级结构距离--> DONE
根据位点数量决定threshold，得到contact map--> DONE
在contact map中寻找所有的n连通--> DONE
结合上一步的得分表，所有n连通中总分最高的n连通作为最终的采样中心sites--> DONE
以采样位点为中心，前后各2个氨基酸共5个位点重设计，装入5-7个氨基酸（pattern随机）--> DONE

-------------------- 序列设计 --------------------
---------- 方法一 ----------
Boltz backbone hallucination幻觉生成若干骨架结构，比对WT结构筛选骨架--> DONE
ligandmpnn设计序列，设置--redesigned_residues，设置NXS/T中的X为P-negative偏置，剩余区域无偏置--> DONE
结构预测，和WT结构对比，设置filter
恢复糖链结构
---------- 方法二 ----------
ESM-LoRA-GLY guided hallucination sequence design
'''
from Bio import SeqIO, pairwise2
from Bio.PDB import Superimposer, PPBuilder
from pathlib import Path
import numpy as np
import pandas as pd
import sys, os
import subprocess
import random
import json
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.util import *
from config import basic_configs

GLY_MOTIFS = ["NXTXX", "NXSXX", "XNXTX", "XNXSX", "XXNXT", "XXNXS",
              "NXTXXX", "NXSXXX", "XNXTXX", "XNXSXX", "XXNXTX", "XXNXSX", "XXXNXT", "XXXNXS",
              "NXTXXXX", "NXSXXXX", "XNXTXXX", "XNXSXXX", "XXNXTXX", "XXNXSXX", "XXXNXTX", "XXXNXSX", "XXXXNXT", "XXXXNXS"]

def find_cliques(contact_map, editable_sites, n):
    N = contact_map.shape[0]
    results = set()

    def backtrack(start, current):

        if len(current) == n:
            sites = sorted(editable_sites[i] for i in current)
            results.add(tuple(sites))
            return

        for next_idx in range(start, N):

            if all(contact_map[next_idx, i] and contact_map[i, next_idx]
                   for i in current):
                backtrack(next_idx + 1, current + [next_idx])

    backtrack(0, [])
    return results

def sample_sites(
    structure_file: str,
    scoring_df: str, 
    chain_id: str,
    num_sites_per_comb: int = 3,
):
    editable_sites = np.array(sorted(pd.read_csv(scoring_df)["Site"].tolist()))
    ca_coords = StructureLoader.get_ca_coords(structure_file=structure_file, chain_id=chain_id)[editable_sites - 1]
    contact_1d = np.abs(editable_sites[:, None] - editable_sites[None, :])
    contact_3d = np.linalg.norm(ca_coords[:, None, :] - ca_coords[None, :, :], axis=-1)

    dist_1d = np.max(contact_1d) / num_sites_per_comb
    dist_3d = np.max(contact_3d) / num_sites_per_comb

    contact_map = (contact_1d >= dist_1d) & (contact_3d >= dist_3d)
    sites_combs = find_cliques(contact_map, editable_sites, num_sites_per_comb)

    df = pd.read_csv(scoring_df)
    sites_scores = {sites: round(np.sum(np.array([df.loc[df["Site"] == s]["SASA"].values[0] for s in sites])), 4) for sites in sites_combs}
    sites_scores = sorted(sites_scores.items(), key=lambda x: x[1], reverse=True)
    sampled_sites = sites_scores[0][0]
    
    return sampled_sites

def backbone_design(
    input_fasta_file: str,
    output_dir: str,
    modify_chain_id: str,
    sampled_sites: list[int],
    backbone_nums: int = 10,
):
    modify_seq_id = (ord(modify_chain_id) - ord("A") + 1)
    seqs = {}
    count_l = count_r = 0

    for rec in SeqIO.parse(input_fasta_file, "fasta"):
        seq_num = int(rec.description.split("copies:")[1])
        count_r += seq_num
        
        if count_l < modify_seq_id and count_r >= modify_seq_id:
            sampled_sites = sorted(sampled_sites, reverse=True)
            seq = str(rec.seq)
            for s in sampled_sites:
                seq = [seq[:s-1-2], seq[s-1-2:s-1+3], seq[s-1+3:]]
                seq[1] = random.sample(GLY_MOTIFS, 1)[0]
                seq = "".join(seq)
        else:
            seq = str(rec.seq)
        
        seqs[rec.id] = (seq, seq_num)
        count_l = count_r
    
    backbone_fasta_file = os.path.join(output_dir, "design", "backbone.fasta")
    with open(backbone_fasta_file, "w") as f:
        for k, (s, n) in seqs.items():
            f.write(f">{k} copies:{n}\n")
            f.write(f"{s}\n")
    
    # MSA
    msa = MsaFileGenerator(input_fasta_file=backbone_fasta_file)
    with open(backbone_fasta_file, "r") as f:
        query = f.read()
    msa.run_mmseqs2(x=query, prefix=f"{output_dir}/design/msa")
    
    # Boltz backbone hallucination
    chain_ptr = 0
    out = []
    for _, (seq, num) in seqs.items():
        ids = FlowList(CHAIN_IDS[chain_ptr : chain_ptr + num])
        out.append({
            "protein": {
                "id": ids,
                "sequence": seq,
                "msa": f"{output_dir}/design/msa/uniref_{len(out)+1}.a3m",
            }
        })
        chain_ptr += num

    with open(f"{output_dir}/design/backbone.yaml", "w") as f:
        yaml.safe_dump({"version": 1, "sequences": out}, f, sort_keys=False, indent=2)
    
    cmd = [
        "boltz", "predict", f"{output_dir}/design/backbone.yaml",
        "--cache", f"../boltz_ckpt",
        "--output_format", "pdb",
        "--out_dir", os.path.join(output_dir, "design"),
        "--diffusion_samples", str(backbone_nums),
    ]
    subprocess.run(cmd, check=True)
    os.makedirs(os.path.join(output_dir, "design", "backbones"), exist_ok=True)
    for i in range(backbone_nums):
        subprocess.run(f"mv {output_dir}/design/boltz_results_backbone/predictions/backbone/backbone_model_{i}.pdb {output_dir}/design/backbones/backbone_{i}.pdb", shell=True)
    subprocess.run(f"rm -r {output_dir}/design/boltz_results_backbone", shell=True)
    subprocess.run(f"rm -f {output_dir}/design/backbone.yaml", shell=True)

def unk_to_gly(
    backbone_pdb_file: str,
):
    with open(backbone_pdb_file, "r") as f:
        text = f.read()
    text = text.replace("UNK", "GLY")
    with open(backbone_pdb_file, "w") as f:
        f.write(text)

def align_backbone(
    wt_pdb_file: str,
    mut_pdb_file: str,
    chain_id: str,
):
    wt_structure, _ = StructureLoader.load_structure(structure_file=wt_pdb_file)
    mut_structure, _ = StructureLoader.load_structure(structure_file=mut_pdb_file)
    wt_model = wt_structure[0][chain_id]
    mut_model = mut_structure[0][chain_id]
    ppb = PPBuilder()
    wt_seq = str(ppb.build_peptides(wt_model)[0].get_sequence())
    mut_seq = str(ppb.build_peptides(mut_model)[0].get_sequence())

    aln = pairwise2.align.globalxx(wt_seq, mut_seq)[0]
    aln1, aln2 = aln.seqA, aln.seqB

    atoms1 = []
    atoms2 = []
    reslist1 = [res for res in wt_model if "CA" in res]
    reslist2 = [res for res in mut_model if "CA" in res]

    i1 = i2 = 0

    for a, b in zip(aln1, aln2):
        if a != "-" and b != "-":
            atoms1.append(reslist1[i1]["CA"])
            atoms2.append(reslist2[i2]["CA"])
        if a != "-":
            i1 += 1
        if b != "-":
            i2 += 1
            
    sup = Superimposer()
    sup.set_atoms(atoms1, atoms2)
    rmsd = sup.rms

    return rmsd

def backbone_filter(
    wt_pdb_file: str,
    output_dir: str,
    chain_id: str
):
    min_rmsd = 10000.
    min_rmsd_mut = None
    for mut in os.listdir(os.path.join(output_dir, "design", "backbones")):
        rmsd = align_backbone(wt_pdb_file=wt_pdb_file, mut_pdb_file=os.path.join(output_dir, "design", "backbones", mut), chain_id=chain_id)
        if rmsd < min_rmsd:
            min_rmsd_mut = mut
            min_rmsd = rmsd
    
    return min_rmsd_mut, min_rmsd
    
def ligandmpnn_redesign(
    input_fasta_file: str,
    backbone_fasta_file: str,
    backbone_pdb_file: str,
    modify_chain_id: str,
    output_dir: str,
    fixed_hotspots: list[int] = None,
):
    unk_to_gly(backbone_pdb_file)
    backbone_seq = FastaLoader.get_sequence(sequence_file=backbone_fasta_file, chain_id=modify_chain_id)
    wt_seq = FastaLoader.get_sequence(sequence_file=input_fasta_file, chain_id=modify_chain_id)

    modify_chain_id = ["A", "B", "C", "D"]
    fixed_aa_idx = [list(re.finditer(wt_seq[i-1:i+3], backbone_seq))[0].start()+1 for i in fixed_hotspots] if fixed_hotspots else list()
    fixed_aa_idx.extend([m.start() + 1 for m in re.finditer("NX", backbone_seq)])
    fixed_aa_idx.extend([m.start() + 3 for m in re.finditer("NX", backbone_seq)])
    fixed_aa_idx = " ".join([f"{c}{i}" for i in fixed_aa_idx for c in modify_chain_id])

    seq_x_idx = " ".join([f"{c}{i + 1}" for i, aa in enumerate(backbone_seq) if aa == "X" for c in modify_chain_id])
    omit_pro_idx = {f"{c}{i}": "P" for i in [m.start() + 2 for m in re.finditer("NX", backbone_seq)] for c in modify_chain_id}
    temp_file = os.path.join(output_dir, "design", "omit_AA_per_residue.json")
    with open(temp_file, "w") as f:
        json.dump(omit_pro_idx, f, indent=2)

    cmd = [
        "python", "../LigandMPNN/run.py",
        "--seed", "42",
        "--model_type", "soluble_mpnn",
        "--checkpoint_soluble_mpnn", "../LigandMPNN/model_params/solublempnn_v_48_002.pt",
        "--pdb_path", backbone_pdb_file,
        "--out_folder", os.path.join(output_dir, "design"),
        # "--redesigned_residues", seq_x_idx,
        "--fixed_residues", fixed_aa_idx,
        "--omit_AA_per_residue", temp_file,
        "--batch_size", "10",
        "--number_of_batches", "1",
        "--verbose", "0",
    ]
    subprocess.run(cmd, check=True)
    subprocess.run(f"mv {output_dir}/design/seqs/{Path(backbone_pdb_file).stem}.fa {output_dir}/design/redesigned.fasta", shell=True)
    subprocess.run(f"rm -r {output_dir}/design/seqs", shell=True)
    subprocess.run(f"rm -r {output_dir}/design/packed", shell=True)
    subprocess.run(f"rm -f {temp_file}", shell=True)

if __name__ == "__main__":
    
    prot = "FTL"
    sampled_sites = sample_sites(
        structure_file=f"/sdata2/WORK/Xuyi/{prot}_test/{prot}.pdb",
        scoring_df=f"/sdata2/WORK/Xuyi/{prot}_test/{prot}_single_points.csv",
        chain_id="A",
        num_sites_per_comb=5,
    )
    print(sampled_sites)
    os.makedirs(os.path.join(f"/sdata2/WORK/Xuyi/{prot}_test", "design"), exist_ok=True)
    backbone_design(
        input_fasta_file=f"/sdata2/WORK/Xuyi/{prot}_test/{prot}.fasta", 
        output_dir=f"/sdata2/WORK/Xuyi/{prot}_test",
        modify_chain_id="A", 
        sampled_sites=sampled_sites,
    )
    mut, rmsd = backbone_filter(
        wt_pdb_file=f"/sdata2/WORK/Xuyi/{prot}_test/{prot}.pdb",
        output_dir=f"/sdata2/WORK/Xuyi/{prot}_test",
        chain_id="A",
    )
    ligandmpnn_redesign(
        input_fasta_file=f"/sdata2/WORK/Xuyi/{prot}_test/{prot}.fasta",
        backbone_fasta_file=f"/sdata2/WORK/Xuyi/{prot}_test/design/backbone.fasta",
        backbone_pdb_file=f"/sdata2/WORK/Xuyi/{prot}_test/design/backbones/{mut}",
        modify_chain_id="A",
        output_dir=f"/sdata2/WORK/Xuyi/{prot}_test",
        fixed_hotspots=[155,156,158,159,161,164,165],
    )
