#----- MMSeqs2 Translated Seach (tBLASTx) Contig Dictionary against the Human Gut Virome DB (GVD) -----#

# query_GVD.smk

#----- Snakemake Set Up -----#
configfile: "./config/phage_pipeline_config.yaml"

# Tools
MMSEQS = config["Tools"]["MMSeqs2"]

rule mmseqs_search_GVD:
	"""
	Query contig dictionary sequences against GVD with translated (tblastx) type search.
	"""
	input:
		queryIdx = os.path.join("results", "mmseqs_phageDB_results", "contig_dictionary_queryDB", "contig_dictionary_queryDB.index"),
		targetIdx = os.path.join("results", "mmseqs_phageDB_results", "GVD_targetDB", "GVD_targetDB.index")
	params:
		queryDB = os.path.join("results", "mmseqs_phageDB_results", "contig_dictionary_queryDB", "contig_dictionary_queryDB"),
		targetDB = os.path.join("results", "mmseqs_phageDB_results", "GVD_targetDB", "GVD_targetDB"),
		resultDB = os.path.join("results", "mmseqs_phageDB_results", "GVD_search_resultDB"),
		outputDB = os.path.join("results", "mmseqs_phageDB_results", "GVD_search_results"),
		tmp = directory(os.path.join("results", "mmseqs_phageDB_results", "tmp_GVD_search"))
	output:
		os.path.join("results", "mmseqs_phageDB_results", "GVD_search_results", "GVD_search_resultDB.index")
	threads: 16
	resources:
		mem_mb = 64000
	shell:
		"""
		ml {MMSEQS}
		mmseqs search \
			{params.queryDB} \
			{params.targetDB} \
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

rule extract_best_GVD_hit:
	"""
	Extract the best hit from the GVD search
	"""
	input:
		resultIdx = os.path.join("results", "mmseqs_phageDB_results", "GVD_search_results", "GVD_search_resultDB.index")
	params:
		resultDB = os.path.join("results", "mmseqs_phageDB_results", "GVD_search_results", "GVD_search_resultDB"),
		bestResultDB = os.path.join("results", "mmseqs_phageDB_results", "GVD_search_results", "GVD_search_bestResultDB")
	threads: 16
	resources:
		mem_mb = 64000
	output:
		os.path.join("results", "mmseqs_phageDB_results", "GVD_search_results", "GVD_search_bestResultDB.index")
	shell:
		"""
		ml {MMSEQS}
		mmseqs filterdb \
			{params.resultDB} \
			{params.bestResultDB} \
			--extract-lines 1 \
			--threads {threads}
		"""

rule convert_best_GVD_results_to_m8:
	"""
	Conver the output of GVD_search_bestResultDB to m8 file
	"""
	input:
		bestResultIdx = os.path.join("results", "mmseqs_phageDB_results", "GVD_search_results", "GVD_search_bestResultDB.index")
	params:
		queryDB = os.path.join("results", "mmseqs_phageDB_results", "contig_dictionary_queryDB", "contig_dictionary_queryDB"),
		targetDB = os.path.join("results", "mmseqs_phageDB_results", "GVD_targetDB", "GVD_targetDB"),
		bestResultDB = os.path.join("results", "mmseqs_phageDB_results", "GVD_search_results", "GVD_search_bestResultDB")
	output:
		os.path.join("results", "mmseqs_phageDB_results", "GVD_search_results", "GVD_results_bestHit.m8")
	threads: 16
	resources:
		mem_mb = 64000
	shell:
		"""
		ml {MMSEQS}
		mmseqs convertalis \
			{params.queryDB} \
			{params.targetDB} \
			{params.bestResultDB} \
			{output} \
			--threads {threads}
		"""
