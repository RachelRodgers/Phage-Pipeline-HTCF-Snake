# Select_Phage_Contigs_for_NT.R

# Store the contig names which had a hit in any of the three phage DBsso they 
#   may be selected from the contig dictionary fasta file (assembly_1k.fasta) 
#   then subsequently searched against the mmseqs-formatted NT database.  
#   This NT search is to remove any contigs with a high identity to human.

options(warn = -1) # Suppress warning messages for clarity

source("./workflow/scripts/snakemake_helpers/snakemake_helpers.R")

requiredPackages <- c("tidyverse")

for (package in requiredPackages) {
  TryInstall(package)
  
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    stop(cat("FATAL: Problem loading R package:", package, 
             "(Generate_Contig_Abundance_Table.R)\n\n"),
         call. = FALSE)
  }
}

filePaths <- "./results/mmseqs_phageDB_results/"

# Currently we have searched against three phage DBs, but ultimately this may
#   be limited to one and negate all the following looping functions.

dataSets <- c("GPD" = "GPD_search_results/GPD_results_bestHit.m8", 
              "GVD" = "GVD_search_results/GVD_results_bestHit.m8", 
              "MGV" = "MGV_search_results/MGV_results_bestHit.m8")

resFiles <- paste0(filePaths, dataSets)
names(resFiles) <- names(dataSets)

m8Columns <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
               "qstart", "qend", "sstart", "send", "evalue", "bitscore")

resData <- map(resFiles, 
               ~ read.delim(.x, header = FALSE, stringsAsFactors = FALSE,
                            colClasses = c("character", rep("NULL", 11))))

resDataUnique <- map(resData, ~ unique(.x))

# bind all elements in resData and select unique
uniqueContigs <- unique(unlist(resDataUnique))

# write out these unique contigs to file
dir.create(path = "./results/phage_contigs_for_NT/", showWarnings = FALSE)

write.table(x = uniqueContigs, 
            file = "./results/phage_contigs_for_NT/phage_contigs_for_NT.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

#----- Save session information -----#

print("Generate_Contig_Abundance_Table.R: Saving session info (retain for debugging).\n")

workingDirectory <- getwd()
savePath <- paste(workingDirectory, "/results/R_session_info/", sep = "")
dir.create(path = savePath, showWarnings = FALSE)
saveFile <- file(paste(savePath, 
                       "Select_Phage_Contigs_for_NT_session_info.txt", 
                       sep = ""))
writeLines(capture.output(sessionInfo()), saveFile)
close(saveFile)
