from pathlib import Path
import click
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from src.prefilters import update_infer, run_prefilters
from src.glycodesign_bb_hallucination import halludesign_bb
from src.glycodesign_esm_hallucination import halludesign_esm

@click.group()
def sugarswitch():
    return

@click.command()
@click.option("--input", type=str, help="fasta file for inference", required=True)
@click.option("--out_dir", default="./output", type=str, help="infer result dir", required=True)
def prefilter(input, out_dir):

    os.makedirs(out_dir, exist_ok=True)

    update_infer(
        input_fasta_file=input,
        output_dir=out_dir,
    )

    run_prefilters(
        input_fasta_file=input,
        input_structure_file=f"{out_dir}/{Path(input).name.split('.')[0]}.pdb",
        output_dir=out_dir,
    )

@click.command()
@click.option("--input", type=str, help="fasta file for inference", required=True)
@click.option("--out_dir", default="./output", type=str, help="infer result dir", required=True)
def design(input, out_dir):
    
    # Method 1: SiteSelection -> MotifInsertion -> BoltzBackboneHallucination -> LigandmpnnSeqDesign
    halludesign_bb(
        wt_structure_file=f"{out_dir}/{Path(input).name.split('.')[0]}.pdb",
        output_dir=out_dir,
    )
    
    # Method 2: SiteSelection -> MotifInsertion -> EsmSeqHallucination
    halludesign_esm(
        input_fasta_file=input,
        wt_structure_file=f"{out_dir}/{Path(input).name.split('.')[0]}.pdb",
        output_dir=out_dir,
    )

sugarswitch.add_command(prefilter)
sugarswitch.add_command(design)

if __name__ == "__main__":
    sugarswitch()