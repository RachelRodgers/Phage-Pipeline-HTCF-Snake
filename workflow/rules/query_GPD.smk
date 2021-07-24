#----- tBLASTx Contig Dictionary against the Gut Phage DB (GPD) -----#

# query_GPD.smk

#----- Snakemake Set Up -----#
configfile: "./config/phage_pipeline_config.yaml"

# Tools
BLAST = config["Tools"]["BLAST"]

# Databases
GPD = config["DBs"]["GPD"]

rule build_GPD_blast_db:
	"""
	Generate a nucleotide blastDB from the Gut Phage Database file.
	"""
	input:
		{GPD}
	params:
		prefix = "GPD",
		out_dir = os.path.join("results", "GPD_nucleotide_db")
	output:
		os.path.join("results", "GPD_nucleotide_db", "GPD.ndb")
	shell:
		"""
		{BLAST}/makeblastdb \
			-in {input} \
			-out {params.prefix} \
			-dbtype nucl
		
		mv GPD.* {params.out_dir}
		"""

rule query_GPD:
	"""
	Run tBLASTx search of filtered contig dictionary against GPD to get annotations.
	"""
	input:
		dbFile = os.path.join("results", "GPD_nucleotide_db", "GPD.ndb"),
		query = os.path.join("results", "contig_dictionary", "assembly_1kb.fasta")
	params:
		blast_dir = os.path.join("results", "GPD_nucleotide_db"),
		prefix = "GPD"
	output:
		os.path.join("results", "GPD_tBLASTx_results", "GPD_contig_dictionary_tBLASTx_results.tab")
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
