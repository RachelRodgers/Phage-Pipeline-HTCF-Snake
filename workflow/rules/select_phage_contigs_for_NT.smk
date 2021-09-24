#----- Select Probable Phage Contigs from the MMSeqs Results -----#

# select_phage_contigs_for_NT.smk

rule select_phage_contigs_for_NT:
	input:
		GPD = os.path.join("results", "mmseqs_phageDB_results", "GPD_search_results", "GPD_results_bestHit.m8"),
		GVD = os.path.join("results", "mmseqs_phageDB_results", "GVD_search_results", "GVD_results_bestHit.m8"),
		MGV = os.path.join("results", "mmseqs_phageDB_results", "MGV_search_results", "MGV_results_bestHit.m8")
	output:
		os.path.join("results", "phage_contigs_for_NT", "phage_contigs_for_NT.txt")
	shell:
		"""
		module load {R}
		Rscript ./workflow/scripts/Select_Phage_Contigs_for_NT.R
		"""
