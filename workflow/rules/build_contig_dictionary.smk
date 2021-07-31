#----- Build Contig Dictionary -----#

# build_contig_dictionary.smk

#----- Snakemake Set Up -----#
configfile: "./config/phage_pipeline_config.yaml"

# Paths
HECATOMB_DIR = config["Paths"]["Hecatomb"]

# Tools
BBTOOLS = config["Tools"]["BBTools"]
METASPADES = config["Tools"]["metaSPAdes"]
MEGAHIT = config["Tools"]["megahit"]
SEQKIT = config["Tools"]["SeqKit"]
FLYE = config["Tools"]["Flye"]

rule get_r2_singletons:
	"""
	Split R2 singletons
	"""
	input:
		os.path.join(HECATOMB_DIR, "results", "QC", "step_7", "{sample}_singletons.s7.out.fastq")
	output:
		os.path.join(HECATOMB_DIR, "results", "QC", "step_7", "{sample}_singletons_R2.out.fastq")
	shell:
		"""
		if grep -q '2:N:' {input}
		then
			grep -A 3 '2:N:' {input} | sed '/^--$/d' > {output}
		else
			touch {output}
		fi
		"""

rule bbnorm:
	"""
	Digital normalization with BBNorm
	"""
	input:
		r1 = os.path.join(HECATOMB_DIR, "results", "QC", "step_7", "{sample}_R1.s7.out.fastq"),
		r2 = os.path.join(HECATOMB_DIR, "results", "QC", "step_7", "{sample}_R2.s7.out.fastq"),
		r1Singletons = os.path.join(HECATOMB_DIR, "results", "QC", "step_7", "{sample}_singletons_R1.out.fastq"),
		r2Singletons = os.path.join(HECATOMB_DIR, "results", "QC", "step_7", "{sample}_singletons_R2.out.fastq")
	output:
		r1 = os.path.join("results", "bbnorm", "{sample}_R1.norm.out.fastq"),
		r2 = os.path.join("results", "bbnorm", "{sample}_R2.norm.out.fastq")
	threads: 4
	shell:
		"""
		ml {BBTOOLS}
		bbnorm.sh \
			in={input.r1} \
			in2={input.r2} \
			extra={input.r1Singletons},{input.r2Singletons} \
			out={output.r1} \
			out2={output.r2} \
			target=20 \
			mindepth=2 \
			t={threads}
		"""

rule megahit:
	"""
	Contig assembly per sample with MEGAHIT
	"""
	input:
		r1 = os.path.join("results", "bbnorm", "{sample}_R1.norm.out.fastq"),
		r2 = os.path.join("results", "bbnorm", "{sample}_R2.norm.out.fastq"),
		r1Singletons = os.path.join(HECATOMB_DIR, "results", "QC", "step_7", "{sample}_singletons_R1.out.fastq"),
		r2Singletons = os.path.join(HECATOMB_DIR, "results", "QC", "step_7", "{sample}_singletons_R2.out.fastq")
	output:
		temp_dir = temp(directory(os.path.join("results", "megahit", "{sample}_TEMP"))),
		contigs = os.path.join("results", "megahit", "{sample}", "final.contigs.fa")
	shell:
		"""
		ml {MEGAHIT}
		megahit \
			-1 {input.r1} \
			-2 {input.r2} \
			-r {input.r1Singletons},{input.r2Singletons} \
			-o {output.temp_dir} \
			--verbose

		mv {output.temp_dir}/final.contigs.fa {output.contigs}
		"""

rule rename_contigs:
	"""
	Edit contig headers so each sample has unique contig names with BBTools rename.sh.
	Prevents flye from crashing.
	"""
	input:
		os.path.join("results", "megahit", "{sample}", "final.contigs.fa")
	output:
		os.path.join("results", "renamed_contigs", "{sample}_contigs.fasta")
	shell:
		"""
		ml {BBTOOLS}
		rename.sh \
			in={input} \
			out={output} \
			prefix={wildcards.sample} \
			addprefix=t
		"""

rule filter_contigs:
	"""
	Filter megahit-assembled contigs to a minimum of 1kb with Seqkit
	"""
	input:
		os.path.join("results", "renamed_contigs", "{sample}_contigs.fasta")
	output:
		os.path.join("results", "filtered_contigs", "{sample}_contigs_1kb.fasta")
	shell:
		"""
		ml {SEQKIT}
		seqkit \
			seq {input} \
			-m 1000 \
			-o {output}
		"""

rule concatenate_renamed_contigs:
        """
        Concatenate all individual (renamed) contig assemblies into one file for flye.
        Hopefully this prevents issues with samples producing 0 contigs.
        """
        input:
                expand(os.path.join("results", "filtered_contigs", "{sample}_contigs_1kb.fasta"), sample = SAMPLES)
        output:
                os.path.join("results", "concatenated_contigs", "concatenated_contigs_1kb.fasta")
	shell:
		"cat {input} > {output}"

rule assemble_contig_dictionary:
	"""
	Generate contig dictionary with Flye
	"""
	input:
		os.path.join("results", "concatenated_contigs", "concatenated_contigs_1kb.fasta")
	params:
		directory(os.path.join("results", "contig_dictionary"))
	output:
		os.path.join("results", "contig_dictionary", "assembly.fasta")
	threads: 8
	shell:
		"""
		ml {FLYE}
		flye \
			--subassemblies {input} \
			-o {params} \
			-t {threads} \
			--meta
		"""

rule filter_contig_dictionary:
	"""
	Filter the contigs in the contig dictionary to a minimum of 1kb with Seqkit.
	I'm not sure how it's possible for flye to generate contigs < 1kb considering
	I only gave it files of 1kb contigs but that's what seems to be happening anyway...
	"""
	input:
		os.path.join("results", "contig_dictionary", "assembly.fasta")
	output:
		os.path.join("results", "contig_dictionary", "assembly_1kb.fasta")
	shell:
                """
                ml {SEQKIT}
                seqkit \
                        seq {input} \
                        -m 1000 \
                        -o {output}
                """
