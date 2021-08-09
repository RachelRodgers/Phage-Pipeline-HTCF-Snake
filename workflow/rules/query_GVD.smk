#----- tBLASTx Contig Dictionary against the Human Gut Virome DB (GVD) -----#

# query_GVD.smk

#----- Snakemake Set Up -----#
configfile: "./config/phage_pipeline_config.yaml"

# Tools
MMSEQS = config["Tools"]["MMSeqs2"]

# Databases
GVD = config["DBs"]["GVD"]

rule mmseqs_GVD:
	"""
	Run mmseqs2 easy-search of contig dictionary against GVD
	"""
	input:
		queryFile = os.path.join("results", "contig_dictionary", "assembly_1kb.fasta"),
		targetFile = {GVD}
	output:
		tmpDir = directory(temp(os.path.join("results", "mmseqs_GVD_results_tmp"))),
		alignmentFile = os.path.join("results", "mmseqs_GVD_results", "mmseqs_GVD_results.m8")
	threads: 16
	resources:
		mem_mb = 64000
	shell:
		"""
		module load {MMSEQS}
		mmseqs easy-search \
			{input.queryFile} \
			{input.targetFile} \
			{output.alignmentFile} \
			{output.tmpDir} \
			--search-type 2 \
			--threads {threads}
		"""
