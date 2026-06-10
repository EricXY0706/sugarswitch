from pathlib import Path
import click
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

@click.group()
def sugarswitch():
    return

@click.command()
@click.option("--input", type=str, help="fasta file for inference", required=True)
@click.option("--input_structure_file", type=str, help="pdb file for inference", required=False)
@click.option("--out_dir", default="./output", type=str, help="infer result dir", required=True)
@click.option("--name", type=str, help="job name", required=False)
@click.option("--chain_id", type=str, help="Chain ID to be glycosyalted", required=True)
@click.option("--functional_hotspots", default="[]", type=str, help="List of positions to avoid modification in string format, e.g. '[1,2,3,'4-10']' ", required=False)
@click.option("--gpu_id", default=0, type=int, help="GPU ID to run the pipeline", required=False)
def prefilter(input, input_structure_file, out_dir, name, chain_id, functional_hotspots, gpu_id):
    
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    from src.prefilters import run_prefilters
    run_prefilters(
        input_fasta_file=input,
        input_structure_file=input_structure_file,
        output_dir=out_dir, 
        name=name,
        protein_chain_id=chain_id,
        functional_hotspots=functional_hotspots,
    )

@click.command()
@click.option("--input", type=str, help="fasta file for inference", required=True)
@click.option("--out_dir", default="./output", type=str, help="infer result dir", required=False)
@click.option("--name", type=str, help="job name", required=False)
@click.option("--chain_id", type=str, help="Chain ID to be glycosyalted", required=True)
@click.option("--num_designs", default=1, type=int, help="number of designs to generate", required=True)
@click.option("--num_gly_sites", default=5, type=int, help="number of glycosylation sites to design", required=True)
@click.option("--num_steps", default=100, type=int, help="number of hallucination steps", required=False)
@click.option("--add_pll_loss", default=True, type=bool, help="whether add pseudo log likelihood loss", required=False)
@click.option("--gpu_id", default=0, type=int, help="GPU ID to run the pipeline", required=False)
def designer(input, out_dir, name, chain_id, num_designs, num_gly_sites, num_steps, add_pll_loss, gpu_id):
    
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    from src.designers import halludesign_esm
    halludesign_esm(
        input_fasta_file=input,
        output_dir=out_dir,
        name=name,
        chain_id=chain_id,
        num_designs=num_designs,
        num_gly_sites=num_gly_sites,
        n_steps=num_steps,
        add_pll_loss=add_pll_loss
    )

@click.command()
@click.option("--input", type=str, help="fasta file for inference", required=True)
@click.option("--out_dir", default="./output", type=str, help="infer result dir", required=True)
def ssbuilder(input, out_dir):
    
    pass

@click.command()
@click.option("--input", type=str, help="fasta file for inference", required=True)
@click.option("--input_structure_file", type=str, help="pdb file for inference", required=False)
@click.option("--out_dir", default="./output", type=str, help="infer result dir", required=False)
@click.option("--name", type=str, help="job name", required=False)
@click.option("--chain_id", type=str, help="Chain ID to be glycosyalted", required=True)
@click.option("--functional_hotspots", default="[]", type=str, help="List of positions to avoid modification in string format, e.g. '[1,2,3,'4-10']' ", required=False)
@click.option("--num_designs", default=1, type=int, help="number of designs to generate", required=False)
@click.option("--num_gly_sites", default=5, type=int, help="number of glycosylation sites to design", required=False)
@click.option("--num_steps", default=100, type=int, help="number of hallucination steps", required=False)
@click.option("--add_pll_loss", default=True, type=bool, help="whether add pseudo log likelihood loss", required=False)
@click.option("--gpu_id", default=0, type=int, help="GPU ID to run the pipeline", required=False)
def pipeline(input, input_structure_file, out_dir, name, chain_id, functional_hotspots, num_designs, num_gly_sites, num_steps, add_pll_loss, gpu_id):
    
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    from src.prefilters import run_prefilters
    from src.designers import halludesign_esm
    
    run_prefilters(
        input_fasta_file=input,
        input_structure_file=input_structure_file,
        output_dir=out_dir, 
        name=name,
        protein_chain_id=chain_id,
        functional_hotspots=functional_hotspots, 
    )
    halludesign_esm(
        input_fasta_file=input,
        output_dir=out_dir,
        name=name,
        chain_id=chain_id,
        num_designs=num_designs,
        num_gly_sites=num_gly_sites,
        n_steps=num_steps,
        add_pll_loss=add_pll_loss,
    )
    
sugarswitch.add_command(prefilter)
sugarswitch.add_command(designer)
sugarswitch.add_command(ssbuilder)
sugarswitch.add_command(pipeline)

if __name__ == "__main__":
    sugarswitch()