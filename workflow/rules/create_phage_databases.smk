#----- Generate MMSeqs DBs from Contig Dictionary & Phage DBs -----#

# create_phage_databases.smk

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
		os.path.join("results", "contig_dictionary_filtered", "assembly_1kb.fasta")
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
