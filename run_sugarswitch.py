from pathlib import Path
import click
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from src.prefilters import run_prefilters
from src.designers import halludesign_esm

@click.group()
def sugarswitch():
    return

@click.command()
@click.option("--input", type=str, help="fasta file for inference", required=True)
@click.option("--out_dir", default="./output", type=str, help="infer result dir", required=True)
def prefilter(input, out_dir):

    os.makedirs(out_dir, exist_ok=True)
    run_prefilters(
        input_fasta_file=input,
        output_dir=out_dir,
    )

@click.command()
@click.option("--input", type=str, help="fasta file for inference", required=True)
@click.option("--out_dir", default="./output", type=str, help="infer result dir", required=False)
@click.option("--num_designs", default=1, type=int, help="number of designs to generate", required=True)
@click.option("--num_gly_sites", default=5, type=int, help="number of glycosylation sites to design", required=True)
def designer(input, out_dir, num_designs, num_gly_sites):
    
    halludesign_esm(
        input_fasta_file=input,
        wt_structure_file=f"{out_dir}/{Path(input).name.split('.')[0]}.pdb",
        output_dir=out_dir,
        num_designs=num_designs,
        num_gly_sites=num_gly_sites,
    )

@click.command()
@click.option("--input", type=str, help="", required=True)
@click.option("--out_dir", default="./output", type=str, help="infer result dir", required=True)
def ssbuilder(input, out_dir):
    
    pass

sugarswitch.add_command(prefilter)
sugarswitch.add_command(designer)
sugarswitch.add_command(ssbuilder)

if __name__ == "__main__":
    sugarswitch()