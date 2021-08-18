#----- Untranslated (blastn) MMSeqs2 Search of Probable Phage Contigs against NT -----#

# query_NT.smk

#----- Snakemake Set Up -----#
configfile: "./config/phage_pipeline_config.yaml"

# Tools
MMSEQS = config["Tools"]["MMSeqs2"]

# Databases
NT = config["DBs"]["NT"]

rule convert_contigs_to_sequenceDB:
	"""
	Convert the phage_contigs_for_NT.fasta file into an mmseqsDB
	"""
	input:
		os.path.join("results", "phage_contigs_for_NT", "phage_contigs_for_NT.fasta")
	params: 
		prefix = "contig_queryDB",
		directory = os.path.join("results", "mmseqs_NT_results", "contig_queryDB")
	output:
		os.path.join("results", "mmseqs_NT_results", "contig_queryDB", "contig_queryDB.index")
	resources:
		mem_mb = 64000
	shell:
		"""
		ml {MMSEQS}
		mmseqs createdb {input} {params.prefix}
		mkdir -p {params.directory}
		mv {params.prefix}* {params.directory}
		"""

rule query_NT:
	"""
	Search the phage contigs against NT with mmseqs2 (nucleotide search).
	(TODO) Once all the NT_resultsDB* files are written out, may want to move into subdirectory of
	mmseqs_NT_results to keep things more organized.
	"""
	input:
		idx = os.path.join("results", "mmseqs_NT_results", "contig_queryDB", "contig_queryDB.index"),
		queryDB = os.path.join("results", "mmseqs_NT_results", "contig_queryDB", "contig_queryDB"),
		targetDB = {NT}
	params:
		alnDB = os.path.join("results", "mmseqs_NT_results", "NT_resultsDB"),
		tmp = directory(os.path.join("results", "mmseqs_NT_results", "tmp_NT_search"))
	output:
		idx = os.path.join("results", "mmseqs_NT_results", "NT_resultsDB.index")
	threads: 16
	resources:
		mem_mb = 200000
	shell:
		"""
		ml {MMSEQS}
		mmseqs search \
			{input.queryDB} \
			{input.targetDB} \
			{params} \
			{params.tmp} \
			-a \
			-e 1.000E-10 \
			--start-sens 1 \
			--sens-steps 3 \
			-s 7 \
			--search-type 3 \
			--threads {threads}
		"""

rule convert_NT_results:
	"""
	Convert the result database created in rule query_NT into a BLAST tab formatted file
	"""
	input:
		idx = os.path.join("results", "mmseqs_NT_results", "NT_resultsDB.index"),
		queryDB = os.path.join("results", "mmseqs_NT_results", "contig_queryDB", "contig_queryDB"),
		targetDB = {NT},
		resultDB = os.path.join("results", "mmseqs_NT_results", "NT_resultsDB")
	output:
		os.path.join("results", "mmseqs_NT_results", "NT_resultsDB.m8")
	threads: 16
	shell:
		"""
		ml {MMSEQS}
		mmseqs convertalis \
			{input.queryDB} \
			{input.targetDB} \
			{input.resultDB} \
			{output} \
			--threads {threads}
		"""
