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

rule search_NT:
	"""
	Search the phage contigs against NT with mmseqs2 (nucleotide search).
	TODO: Option for force-reuse upon failure
	"""
	input:
		queryIdx = os.path.join("results", "mmseqs_NT_results", "contig_queryDB", "contig_queryDB.index"),
		targetDB = {NT}
	params:
		queryDB = os.path.join("results", "mmseqs_NT_results", "contig_queryDB", "contig_queryDB"),
		resultDB = os.path.join("results", "mmseqs_NT_results", "NT_resultsDB"),
		outputDB = os.path.join("results", "mmseqs_NT_results", "NT_search_results"),
		tmp = directory(os.path.join("results", "mmseqs_NT_results", "tmp_NT_search"))
	output:
		os.path.join("results", "mmseqs_NT_results", "NT_search_results", "NT_resultsDB.index")
	threads: 16
	resources:
		mem_mb = 245000
	shell:
		"""
		ml {MMSEQS}
		mmseqs search \
			{params.queryDB} \
			{input.targetDB} \
			{params.resultDB} \
			{params.tmp} \
			-a \
			-e 1.000E-10 \
			--start-sens 1 \
			--sens-steps 3 \
			-s 7 \
			--search-type 3 \
			--threads {threads}
		mkdir -p {params.outputDB}
		mv {params.resultDB}* {params.outputDB}
		"""

rule extract_best_NT_hit:
	"""
	Extract the best hit from the NT search
	"""
	input:
		resultIdx = os.path.join("results", "mmseqs_NT_results", "NT_search_results", "NT_resultsDB.index")
	params:
		resultDB = os.path.join("results", "mmseqs_NT_results", "NT_search_results", "NT_resultsDB"),
		bestResultDB = os.path.join("results", "mmseqs_NT_results", "NT_search_results", "NT_search_bestResultDB")
	output:
		os.path.join("results", "mmseqs_NT_results", "NT_search_results", "NT_search_bestResultDB.index")
	threads: 16
	resources:
		mem_mb = 64000
	shell:
		"""
		ml {MMSEQS}
		mmseqs filterdb \
			{params.resultDB} \
			{params.bestResultDB} \
			--extract-lines 1 \
			--threads {threads}
		"""

rule convert_best_NT_results_to_m8:
	"""
	Convert the output of NT_search_bestResultDB to m8 file
	"""
	input:
		bestResultIdx = os.path.join("results", "mmseqs_NT_results", "NT_search_results", "NT_search_bestResultDB.index"),
		targetDB = {NT}
	params:
		queryDB = os.path.join("results", "mmseqs_NT_results", "contig_queryDB", "contig_queryDB"),
		bestResultDB = os.path.join("results", "mmseqs_NT_results", "NT_search_results", "NT_search_bestResultDB")
	output:
		os.path.join("results", "mmseqs_NT_results", "NT_search_results", "NT_results_bestHit.m8")
	threads:
		16
	resources:
		mem_mb = 64000
	shell:
		"""
		ml {MMSEQS}
		mmseqs convertalis \
			{params.queryDB} \
			{input.targetDB} \
			{params.bestResultDB} \
			{output} \
			--threads {threads} 
		"""
