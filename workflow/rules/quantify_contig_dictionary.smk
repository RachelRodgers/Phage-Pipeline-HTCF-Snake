#----- Build Reference From Contig Dictionary for Quantification By Mapping -----#

# quantify_contig_dictionary.smk

#----- Snakemake Set Up -----#
configfile: "./config/phage_pipeline_config.yaml"

# Paths
HECATOMB_DIR = config["Paths"]["Hecatomb"]

# Tools
BBTOOLS = config["Tools"]["BBTools"]
R = config["Tools"]["R"]

rule build_reference_from_contig_dictionary:
	"""
	Build reference from contig dictionary (concatenated individual assemblies) 
	for BBMap for the quantification by mapping (RPKM step)
	"""
	input:
		os.path.join("results", "contig_dictionary_filtered", "assembly_1kb.fasta")
	output:
		os.path.join("results", "ref", "genome", "1", "summary.txt")
	shell:
		"""
		bash {BBTOOLS}bbmap.sh \
			ref={input} \
			build=1 \
			path=./results
		"""
	
rule quantification_by_mapping:
	"""
	Map QC'd reads to contig dictionary with BBMap
	"""
	input:
		r1 = os.path.join(HECATOMB_DIR, "results", "QC", "step_7", "{sample}_R1.s7.out.fastq"),
                r2 = os.path.join(HECATOMB_DIR, "results", "QC", "step_7", "{sample}_R2.s7.out.fastq"),
		summary = os.path.join("results", "ref", "genome", "1", "summary.txt")
	output:
		sam = os.path.join("results", "quantification", "{sample}.aln.sam.gz"),
		rpkm = os.path.join("results", "quantification", "{sample}.rpkm"),
		covstats = os.path.join("results", "quantification", "{sample}.covstats")
	threads: 8
	shell:
		"""
		bash {BBTOOLS}bbmap.sh \
			path=./results \
			build=1 \
			in={input.r1} \
			in2={input.r2} \
			out={output.sam} \
			rpkm={output.rpkm} \
			covstats={output.covstats} \
			kfilter=22 \
			subfilter=15 \
			maxindel=80 \
			ambiguous=random \
			physcov=t \
			t={threads}
		"""

# rule generate_contig_abundances:
# 	"""
# 	Calculate contig abundances using TPM
# 	"""
# 	input:
# 		rpkm = expand(os.path.join("results", "quantification", "{sample}.rpkm"), sample = SAMPLES),
# 		covstats = expand(os.path.join("results", "quantification", "{sample}.covstats"), sample = SAMPLES)
# 	output:
# 		os.path.join("results", "contig_abundance_table", "contig_abundance_table.txt")
# 	shell:
# 		"""
# 		ml {R}
# 		Rscript ./workflow/scripts/Generate_Contig_Abundance_Table.R
# 		"""
