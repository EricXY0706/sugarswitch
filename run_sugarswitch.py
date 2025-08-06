from prefilters import update_infer, run_prefilters
import click
import os

@click.group()
def sugarswitch():
    return

@click.command()
@click.option("--input", type=str, help="fasta file for inference")
@click.option("--input_structure", type=str, help="pdb or cif file for inference")
@click.option("--out_dir", default="./output", type=str, help="infer result dir")
def prefilter(input, input_structure, out_dir):

    os.makedirs(out_dir, exist_ok=True)

    # update_infer(
    #     input_fasta_file=input,
    #     output_dir=out_dir,
    # )

    run_prefilters(
        input_fasta_file=input,
        input_structure_file=input_structure,
        input_alignment_file=f"{out_dir}/msa/uniref.a3m",
        output_dir=out_dir
    )

sugarswitch.add_command(prefilter)


if __name__ == "__main__":
    sugarswitch()