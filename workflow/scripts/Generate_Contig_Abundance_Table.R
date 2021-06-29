# Generate_Contig_Abundance_Table.R

# Calculate TPM for each sample and add coverage statistics for each contig.
#   Combine TPM/stats for each contig for every sample into a "long-format"
#   output table that can be used for analyses.

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

dataFilePath <- "./results/quantification/"

#----- Read-in CovStats Files -----#

covstatsFiles <- list.files(path = dataFilePath, pattern = "*.covstats")
covstatsFilePaths <- paste0(dataFilePath, covstatsFiles)

covstatsSampleNames <- map_chr(.x = covstatsFiles,
                               .f = ~ str_remove(.x, ".covstats$"))
names(covstatsFilePaths) <- covstatsSampleNames

covstatsFileData <- map(.x = covstatsFilePaths,
                        .f = ~ read.delim(file = .x, header = TRUE,
                                          stringsAsFactors = FALSE,
                                          check.names = FALSE))

#----- Process RPKM Files -----#

rpkmFiles <- list.files(path = dataFilePath, pattern = "*.rpkm")
rpkmFilePaths <- paste0(dataFilePath, rpkmFiles)

sampleNames <- map_chr(.x = rpkmFiles, .f = ~ str_remove(.x, ".rpkm$"))
names(rpkmFilePaths) <- sampleNames

# read in each rpkm file, ignoring the first four comment lines (each file has 
#   a total of 5 comment lines including column names)
rpkmFileData <- map(.x = rpkmFilePaths,
                    .f = ~ read.delim(file = .x, skip = 4, header = TRUE,
                                      stringsAsFactors = FALSE,
                                      check.names = FALSE))

# Calculate Transcripts per Million
# For each sample:
#   1. Divide reads for each contig by contig length in kilobases to get RPK 
#     (reads per thousand)
#   2. Calculate the size factor for the sample by summing the RPK for each
#     contig and dividing by one million.
#   3. Calculate TPM for each contig by dividing each contig's RPK by the size
#     factor.

tpmDFList <- vector(mode = "list", length = length(rpkmFileData))

for (i in 1:length(rpkmFileData)) {
  
  currentData <- rpkmFileData[[i]]
  currentSample <- names(rpkmFileData)[i]
  
  currentCovstatsData <- covstatsFileData[[currentSample]] %>% 
    rename(contig_id = `#ID`) %>% 
    select(contig_id, Avg_fold, Ref_GC, Covered_percent, Covered_bases, 
           Median_fold)
  
  rpk <- currentData %>% 
    select(`#Name`, Length, Reads, RPKM, FPKM) %>% 
    mutate(RPK = Reads/(Length/1000)) %>% 
    rename(contig_id = `#Name`)
  
  sizeFactor <- sum(rpk$RPK/1000000)
  
  tpm <- rpk %>% 
    mutate(TPM = RPK/sizeFactor,
           sample_id = currentSample) %>% 
    select(sample_id, everything())
  
  tpmFinal <- merge(tpm, currentCovstatsData, by = "contig_id") %>% 
    select(sample_id, contig_id, everything())
  
  tpmDFList[[i]] <- tpmFinal
  names(tpmDFList)[i] <- currentSample
  
}

# Merge all the individual TPM data frames into one large table
longTPM <- Reduce(f = rbind, x = tpmDFList)

dir.create(path = "./results/contig_abundance_table/", showWarnings = FALSE)

write.table(x = longTPM, 
            file = "./results/contig_abundance_table/contig_abundance_table.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

#----- Save session information -----#

print("Generate_Contig_Abundance_Table.R: Saving session info (retain for debugging).\n")

workingDirectory <- getwd()
savePath <- paste(workingDirectory, "/results/R_session_info/", sep = "")
dir.create(path = savePath, showWarnings = FALSE)
saveFile <- file(paste(savePath, 
                       "Generate_Contig_Abundance_Table_session_info.txt", 
                       sep = ""))
writeLines(capture.output(sessionInfo()), saveFile)
close(saveFile)
