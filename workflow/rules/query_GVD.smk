#----- tBLASTx Contig Dictionary against the Human Gut Virome DB (GVD) -----#

# query_GVD.smk

#----- Snakemake Set Up -----#
configfile: "./config/phage_pipeline_config.yaml"

# Tools
BLAST = config["Tools"]["BLAST"]

# Databases
GVD = config["DBs"]["GVD"]

rule build_GVD_blast_db:
	"""
	Generate a nucleotide blastDB from the GVD file.
	"""
	input:
		{GVD}
	params:
		prefix = "GVD",
		out_dir = os.path.join("results", "GVD_nucleotide_db")
	output:
		os.path.join("results", "GVD_nucleotide_db", "GVD.ndb")
	shell:
		"""
		{BLAST}/makeblastdb \
			-in {input} \
			-out {params.prefix} \
			-dbtype nucl

		mv GVD.* {params.out_dir}
		"""

rule query_GVD:
	"""
	Run tBLASTx search of filtered contig dictionary against GVD to get annotations.
	"""
	input:
		dbFile = os.path.join("results", "GVD_nucleotide_db", "GVD.ndb"),
		query = os.path.join("results", "contig_dictionary", "assembly_1kb.fasta")
	params:
		blast_dir = os.path.join("results", "GVD_nucleotide_db"),
		prefix = "GVD"
	output:
		os.path.join("results", "GVD_tBLASTx_results", "GVD_contig_dictionary_tBLASTx_results.tab")
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
