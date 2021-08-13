#----- Select Probable Phage Contigs from the MMSeqs Results -----#

# select_phage_contigs_for_NT.smk

rule select_phage_contigs_for_NT:
	input:
		GVD = os.path.join("results", "mmseqs_GVD_results", "mmseqs_GVD_results.m8"),
		GPD = os.path.join("results", "mmseqs_GPD_results", "mmseqs_GPD_results.m8"),
		MGV = os.path.join("results", "mmseqs_MGV_results", "mmseqs_MGV_results.m8")
	output:
		os.path.join("results", "phage_contigs_for_NT", "phage_contigs_for_NT.txt")
	shell:
		"""
		module load {R}
		Rscript ./workflow/scripts/Select_Phage_Contigs_for_NT.R
		"""
