#!/bin/bash

# Run_CT2_Virus_Enriched.sh

# Used with phage pipeline snakemake

# Read in the CT2 template path and environment path.
templatePath=$(cat ./workflow/scripts/ct2Temp.txt)
environmentPath=$(cat ./workflow/scripts/ct2Env.txt)

eval "$(conda shell.bash hook)"
conda activate ${environmentPath}/cenote-taker2_env/

python ${environmentPath}/Cenote-Taker2/run_cenote-taker2.py \
	--contigs ./results/contig_dictionary/assembly_1kb.fasta \
	--run_title ct2_annotated_contig_dictionary \
	--template_file ${templatePath} \
	--prune_prophage False \
	--mem 32 \
	--cpu 16 \
	--virus_domain_db standard \
	--circ_minimum_hallmark_genes 0 \
	--lin_minimum_hallmark_genes 1
