from prefilters import update_infer, run_prefilters
from pathlib import Path
import click
import os

@click.group()
def sugarswitch():
    return

@click.command()
@click.option("--input", type=str, help="fasta file for inference", required=True)
@click.option("--out_dir", default="./output", type=str, help="infer result dir", required=True)
def prefilter(input, out_dir):

    os.makedirs(out_dir, exist_ok=True)

    # update_infer(
    #     input_fasta_file=input,
    #     output_dir=out_dir,
    # )

    run_prefilters(
        input_fasta_file=input,
        input_structure_file=f"{out_dir}/{Path(input).name.split('.')[0]}.pdb",
        input_alignment_file=f"{out_dir}/msa/uniref.a3m",
        output_dir=out_dir
    )

sugarswitch.add_command(prefilter)


if __name__ == "__main__":
    sugarswitch()