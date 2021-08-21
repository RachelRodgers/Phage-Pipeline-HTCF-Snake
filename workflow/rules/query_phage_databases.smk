#----- Query Contig Dictionary against Phage Databases using MMSeqs -----#

# query_phage_databases.smk

#----- Snakemake Set Up -----#
configfile: "./config/phage_pipeline_config.yaml"

# Tools
MMSEQS = config["Tools"]["MMSeqs2"]

# Databases
GPD = config["DBs"]["GPD"]
GVD = config["DBs"]["GVD"]
MGV = config["DBs"]["MGV"]

rule convert_contig_dict_to_queryDB:
	"""
	Convert the contig dictionary sequences into a queryDB for mmseqs search
	"""
	input:
		os.path.join("results", "contig_dictionary", "assembly_1kb.fasta")
	params:
		prefix = "contig_dictionary_queryDB",
		directory = os.path.join("results", "mmseqs_phageDB_results", "contig_dictionary_queryDB")
	output:
		os.path.join("results", "mmseqs_phageDB_results", "contig_dictionary_queryDB", "contig_dictionary_queryDB.index")
	shell:
		"""
		ml {MMSEQS}
		mmseqs createdb {input} {params.prefix}
		mkdir -p {params.directory}
		mv {params.prefix}* {params.directory}
		"""

rule convert_GPD_to_queryDB:
	"""
	Convert the GPD database file into a targetDB for mmseqs search
	"""
	input:
		{GPD}
	params:
		prefix = "GPD_targetDB",
		directory = os.path.join("results", "mmseqs_phageDB_results", "GPD_targetDB")
	output:
		os.path.join("results", "mmseqs_phageDB_results", "GPD_targetDB", "GPD_targetDB.index")
	shell:
		"""
		ml {MMSEQS}
		mmseqs createdb {input} {params.prefix}
		mkdir -p {params.directory}
		mv {params.prefix}* {params.directory}
		"""

rule convert_GVD_to_queryDB:
	"""
	Convert the GVD database file into a targetDB for mmseqs search
	"""
	input:
		{GVD}
	params:
		prefix = "GVD_targetDB",
		directory = os.path.join("results", "mmseqs_phageDB_results", "GVD_targetDB")
	output:
		os.path.join("results", "mmseqs_phageDB_results", "GVD_targetDB", "GVD_targetDB.index")
	shell:
		"""
		ml {MMSEQS}
		mmseqs createdb {input} {params.prefix}
		mkdir -p {params.directory}
		mv {params.prefix}* {params.directory}
		"""

rule convert_MGV_to_queryDB:
	"""
	Convert the MGV database file into a targetDB for mmseqs search
	"""
	input:
		{MGV}
	params:
		prefix = "MGV_targetDB",
		directory = os.path.join("results", "mmseqs_phageDB_results", "MGV_targetDB")
	output:
		os.path.join("results", "mmseqs_phageDB_results", "MGV_targetDB", "MGV_targetDB.index")
	shell:
		"""
		ml {MMSEQS}
		mmseqs createdb {input} {params.prefix}
		mkdir -p {params.directory}
		mv {params.prefix}* {params.directory}
		"""

rule mmseqs_search_GPD:
	"""
	Query contig dictionary sequences against GPD with translated (tblastx) type search.
	"""
	input:
		queryIdx = os.path.join("results", "mmseqs_phageDB_results", "contig_dictionary_queryDB", "contig_dictionary_queryDB.index"),
		GPDIdx = os.path.join("results", "mmseqs_phageDB_results", "GPD_targetDB", "GPD_targetDB.index"),
		queryDB = os.path.join("results", "mmseqs_phageDB_results", "contig_dictionary_queryDB", "contig_dictionary_queryDB"),
		targetDB = os.path.join("results", "mmseqs_phageDB_results", "GPD_targetDB", "GPD_targetDB")
	params:
		resultDB = os.path.join("results", "mmseqs_phageDB_results", "GPD_search_resultDB"),
		outputDB = os.path.join("results", "mmseqs_phageDB_results", "GPD_search_results"),
		tmp = directory(os.path.join("results", "mmseqs_phageDB_results", "tmp_GPD_search"))
	output:
		os.path.join("results", "mmseqs_phageDB_results", "GPD_search_results", "GPD_search_resultDB.index")
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

rule extract_best_GPD_hit:
	"""
	Extract the best hit from the GPD search
	"""
	input:
		resultIdx = os.path.join("results", "mmseqs_phageDB_results", "GPD_search_results", "GPD_search_resultDB.index")
	params:
		resultDB = os.path.join("results", "mmseqs_phageDB_results", "GPD_search_results", "GPD_search_resultDB"),
		bestResultDB = os.path.join("results", "mmseqs_phageDB_results", "GPD_search_results", "GPD_search_bestResultDB")
	threads: 16
	resources:
		mem_mb = 64000
	output:
		os.path.join("results", "mmseqs_phageDB_results", "GPD_search_results", "GPD_search_bestResultDB.index")
	shell:
		"""
		ml {MMSEQS}
		mmseqs filterdb \
			{params.resultDB} \
			{params.bestResultDB} \
			--extract-lines 1 \
			--threads {threads}
		"""

rule convert_best_GPD_results_to_m8:
	"""
	Conver the output of GPD_search_bestResultDB to m8 file
	"""
	input:
		bestResultIdx = os.path.join("results", "mmseqs_phageDB_results", "GPD_search_results", "GPD_search_bestResultDB.index"),
		queryDB = os.path.join("results", "mmseqs_phageDB_results", "contig_dictionary_queryDB", "contig_dictionary_queryDB"),
		targetDB = os.path.join("results", "mmseqs_phageDB_results", "GPD_targetDB", "GPD_targetDB")
	params:
		bestResultDB = os.path.join("results", "mmseqs_phageDB_results", "GPD_search_results", "GPD_search_bestResultDB")
	output:
		os.path.join("results", "mmseqs_phageDB_results", "GPD_search_results", "GPD_results_bestHit.m8")
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
