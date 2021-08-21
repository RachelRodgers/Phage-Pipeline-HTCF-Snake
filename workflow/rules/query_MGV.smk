#----- MMSeqs2 Translated Seach (tBLASTx) Contig Dictionary against the Metagenomic Gut Virus (MGV) Catalog -----#

# query_MGV.smk

#----- Snakemake Set Up -----#
configfile: "./config/phage_pipeline_config.yaml"

# Tools
MMSEQS = config["Tools"]["MMSeqs2"]

rule mmseqs_search_MGV:
	"""
	Query contig dictionary sequences against MGV with translated (tblastx) type search.
	"""
	input:
		queryIdx = os.path.join("results", "mmseqs_phageDB_results", "contig_dictionary_queryDB", "contig_dictionary_queryDB.index"),
		targetIdx = os.path.join("results", "mmseqs_phageDB_results", "MGV_targetDB", "MGV_targetDB.index"),
		queryDB = os.path.join("results", "mmseqs_phageDB_results", "contig_dictionary_queryDB", "contig_dictionary_queryDB"),
		targetDB = os.path.join("results", "mmseqs_phageDB_results", "MGV_targetDB", "MGV_targetDB")
	params:
		resultDB = os.path.join("results", "mmseqs_phageDB_results", "MGV_search_resultDB"),
		outputDB = os.path.join("results", "mmseqs_phageDB_results", "MGV_search_results"),
		tmp = directory(os.path.join("results", "mmseqs_phageDB_results", "tmp_MGV_search"))
	output:
		os.path.join("results", "mmseqs_phageDB_results", "MGV_search_results", "MGV_search_resultDB.index")
	threads: 16
	resources:
		mem_mb = 64000
	shell:
		"""
		ml {MMSEQS}
		mmseqs search \
			{input.queryDB} \
			{input.targetDB} \
			{params.resultDB} \
			{params.tmp} \
			-a \
			--start-sens 1 \
			--sens-steps 3 \
			-s 7 \
			--search-type 2 \
			--threads {threads}
		mkdir -p {params.outputDB}
		mv {params.resultDB}* {params.outputDB}
		"""

rule extract_best_MGV_hit:
	"""
	Extract the best hit from the MGV search
	"""
	input:
		resultIdx = os.path.join("results", "mmseqs_phageDB_results", "MGV_search_results", "MGV_search_resultDB.index")
	params:
		resultDB = os.path.join("results", "mmseqs_phageDB_results", "MGV_search_results", "MGV_search_resultDB"),
		bestResultDB = os.path.join("results", "mmseqs_phageDB_results", "MGV_search_results", "MGV_search_bestResultDB")
	threads: 16
	resources:
		mem_mb = 64000
	output:
		os.path.join("results", "mmseqs_phageDB_results", "MGV_search_results", "MGV_search_bestResultDB.index")
	shell:
		"""
		ml {MMSEQS}
		mmseqs filterdb \
			{params.resultDB} \
			{params.bestResultDB} \
			--extract-lines 1 \
			--threads {threads}
		"""

rule convert_best_MGV_results_to_m8:
	"""
	Conver the output of MGV_search_bestResultDB to m8 file
	"""
	input:
		bestResultIdx = os.path.join("results", "mmseqs_phageDB_results", "MGV_search_results", "MGV_search_bestResultDB.index"),
		queryDB = os.path.join("results", "mmseqs_phageDB_results", "contig_dictionary_queryDB", "contig_dictionary_queryDB"),
		targetDB = os.path.join("results", "mmseqs_phageDB_results", "MGV_targetDB", "MGV_targetDB")
	params:
		bestResultDB = os.path.join("results", "mmseqs_phageDB_results", "MGV_search_results", "MGV_search_bestResultDB")
	output:
		os.path.join("results", "mmseqs_phageDB_results", "MGV_search_results", "MGV_results_bestHit.m8")
	threads: 16
	resources:
		mem_mb = 64000
	shell:
		"""
		ml {MMSEQS}
		mmseqs convertalis \
			{input.queryDB} \
			{input.targetDB} \
			{params.bestResultDB} \
			{output} \
			--threads {threads}
		"""