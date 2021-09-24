#----- Untranslated (blastn) MMSeqs2 Search of Probable Phage Contigs against NT -----#

# query_NT.smk

#----- Snakemake Set Up -----#
configfile: "./config/phage_pipeline_config.yaml"

# Tools
MMSEQS = config["Tools"]["MMSeqs2"]

# Databases
NT = config["DBs"]["NT"]

rule mmseqs_NT:
	"""
	Run mmseqs2 easy-search of probable phage contigs against NT
	"""
	input:
		queryFile = os.path.join("results", "phage_contigs_for_NT", "phage_contigs_for_NT.fasta"),
		targetFile = {NT}
	output:
		tmpDir = directory(temp(os.path.join("results", "mmseqs_NT_results_tmp"))),
		alignmentFile = os.path.join("results", "mmseqs_NT_results", "mmseqs_NT_results.m8")
	threads: 16
	resources:
		mem_mb = 96000
	shell:
		"""
		module load {MMSEQS}
		mmseqs easy-search \
			{input.queryFile} \
			{input.targetFile} \
			{output.alignmentFile} \
			{output.tmpDir} \
			-e 1.000E-10 \
			--search-type 3 \
			--threads {threads}
		"""
