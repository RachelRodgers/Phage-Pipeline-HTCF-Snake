#----- tBLASTx Contig Dictionary against the Metagenomic Gut Virus Catalog (MGV) -----#

# query_MGV.smk

#----- Snakemake Set Up -----#
configfile: "./config/phage_pipeline_config.yaml"

# Tools
BLAST = config["Tools"]["BLAST"]

# Databases
MGV = config["DBs"]["MGV"]

rule build_MGV_blast_db:
	"""
	Generate a nucleotide blastDB from the Metagenomic Gut Virus Catalog file.
	"""
	input:
		{MGV}
	params:
		prefix = "MGV",
		out_dir = os.path.join("results", "MGV_nucleotide_db")
	output:
		os.path.join("results", "MGV_nucleotide_db", "MGV.ndb")
	shell:
		"""
		{BLAST}/makeblastdb \
			-in {input} \
			-out {params.prefix} \
			-dbtype nucl
		
		mv MGV.* {params.out_dir}
		"""

rule query_MGV:
	"""
	Run tBLASTx search of filtered contig dictionary against MGV to get annotations.
	"""
	input:
		dbFile = os.path.join("results", "MGV_nucleotide_db", "MGV.ndb"),
		query = os.path.join("results", "contig_dictionary", "assembly_1kb.fasta")
	params:
		blast_dir = os.path.join("results", "MGV_nucleotide_db"),
		prefix = "MGV"
	output:
		os.path.join("results", "MGV_tBLASTx_results", "MGV_contig_dictionary_tBLASTx_results.tab")
	threads: 8
	shell:
		"""
		{BLAST}/tblastx \
			-db {params.blast_dir}/{params.prefix} \
			-query {input.query} \
			-out {output} \
			-outfmt "6 qseqid sseqid pident length evalue bitscore" \
			-evalue 0.001 \
			-num_threads {threads}
		"""
