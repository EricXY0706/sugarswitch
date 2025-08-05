from prefilters import run_prefilters
import click

@click.group()
def sugarswitch():
    return

@click.command()
@click.option("--input_fasta", type=str, help="fasta file for inference")
@click.option("--input_structure", type=str, help="pdb or cif file for inference")
@click.option("--input_alignment", type=str, help="a3m file for inference")
@click.option("--out_dir", default="./output", type=str, help="infer result dir")
def prefilter(input_fasta, input_structure, input_alignment, out_dir):

    run_prefilters(
        input_fasta_file=input_fasta,
        input_structure_file=input_structure,
        input_alignment_file=input_alignment,
        output_dir=out_dir
    )

sugarswitch.add_command(prefilter)


if __name__ == "__main__":
    sugarswitch()