from pathlib import Path
import click
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

@click.group()
def sugarswitch():
    return

@click.command()
@click.option("--input_fasta_file", type=str, help="fasta file for inference", required=True)
@click.option("--input_structure_file", type=str, help="pdb file for inference", required=False)
@click.option("--out_dir", default="./output", type=str, help="infer result dir", required=True)
@click.option("--name", type=str, help="job name", required=False)
@click.option("--chain_id", type=str, help="Chain ID to be glycosyalted", required=True)
@click.option("--functional_hotspots", default="[]", type=str, help="List of positions to avoid modification in string format, e.g. '[1,2,3,'4-10']' ", required=False)
@click.option("--enable_glycan_grafting", default=True, type=bool, help="Whether enables glycan grafting to access cavity volumn to hold glycans", required=False)
def prefilter(input_fasta_file, input_structure_file, out_dir, name, chain_id, functional_hotspots, enable_glycan_grafting):
    
    from src.prefilters import run_prefilters
    run_prefilters(
        input_fasta_file=input_fasta_file,
        input_structure_file=input_structure_file,
        output_dir=out_dir, 
        name=name,
        protein_chain_id=chain_id,
        functional_hotspots=functional_hotspots,
        enable_glycan_grafting=enable_glycan_grafting,
    )

@click.command()
@click.option("--input_fasta_file", type=str, help="fasta file for inference", required=True)
@click.option("--out_dir", default="./output", type=str, help="infer result dir", required=True)
def ssbuilder(input_fasta_file, out_dir):
    
    pass

@click.command()
@click.option("--input_fasta_file", type=str, help="fasta file for inference", required=True)
@click.option("--input_structure_file", type=str, help="pdb file for inference", required=False)
@click.option("--out_dir", default="./output", type=str, help="infer result dir", required=False)
@click.option("--name", type=str, help="job name", required=False)
@click.option("--chain_id", type=str, help="Chain ID to be glycosyalted", required=True)
@click.option("--functional_hotspots", default="[]", type=str, help="List of positions to avoid modification in string format, e.g. '[1,2,3,'4-10']' ", required=False)
@click.option("--enable_glycan_grafting", default=True, type=bool, help="Whether enables glycan grafting to access cavity volumn to hold glycans", required=False)
@click.option("--num_patterns", default=1, type=int, help="number of combinations of sites", required=False)
@click.option("--num_candidates_per_pattern", default=10, type=int, help="number of candidates per combination of sites", required=False)
@click.option("--num_designs_per_pattern", default=1, type=int, help="number of final designs out of candidates per combination of sites", required=False)
@click.option("--num_gly_sites", default=3, type=int, help="number of glycosylation sites to design", required=False)
@click.option("--add_pll_loss", default=True, type=bool, help="whether add pseudo log likelihood loss", required=False)
@click.option("--predict_structure", default=True, type=bool, help="whether predict designed glycoprotein structure", required=False)
def designer(input_fasta_file, input_structure_file, out_dir, name, chain_id, functional_hotspots, num_patterns, num_candidates_per_pattern, num_designs_per_pattern, num_gly_sites, add_pll_loss, enable_glycan_grafting, predict_structure):
    
    from src.prefilters import run_prefilters
    from src.designers import halludesign_esm
    
    conservation_df, coupling_stength, interaction_dict, rsasa_index_dict = run_prefilters(
        input_fasta_file=input_fasta_file,
        input_structure_file=input_structure_file,
        output_dir=out_dir, 
        name=name,
        protein_chain_id=chain_id,
        functional_hotspots=functional_hotspots,
        enable_glycan_grafting=enable_glycan_grafting,
    )
    halludesign_esm(
        input_fasta_file=input_fasta_file,
        output_dir=out_dir,
        conservation_df=conservation_df,
        coupling_stength=coupling_stength,
        interaction_dict=interaction_dict,
        rsasa_index_dict=rsasa_index_dict,
        name=name,
        chain_id=chain_id,
        num_patterns=num_patterns,
        num_candidates_per_pattern=num_candidates_per_pattern,
        num_designs_per_pattern=num_designs_per_pattern,
        num_gly_sites=num_gly_sites,
        add_pll_loss=add_pll_loss,
        predict_structure=predict_structure,
    )
    
sugarswitch.add_command(prefilter)
sugarswitch.add_command(ssbuilder)
sugarswitch.add_command(designer)

if __name__ == "__main__":
    sugarswitch()