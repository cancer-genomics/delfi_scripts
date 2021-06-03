#--------------------------------------------------------------------------
# The purpose of this script is to generate mitochondrial representation
# statistics for 19 additional samples that were requested on Nov 12, 2018
#--------------------------------------------------------------------------

library(readr)
library(ggplot2)

#-------------------------------------------------------------------
# Reading in the list of samples that stephen is currently using
# with  DELFI (Nov 2, 2018)
#-------------------------------------------------------------------------------------------------------------------
sample_table <- readLines("/dcl01/scharpf1/data/dbruhm/stephen/pa-pipeline/summary/samples-to-summarize-nov12.txt")
sample_table <- gsub("_Amended", "", sample_table)
#-------------------------------------------------------------------------------------------------------------------


#----------------------------------------------------
# Compiling the mitochondria stats from the Bowtie2
# realignment into a table
#------------------------------------------------------------------------------------------------
files <- list.files("/dcl01/scharpf1/data/dbruhm/stephen/mito_realignment/outDir",
                     pattern = "reads_mapped_to_mito.rds", full.names = TRUE, recursive = TRUE)

rm(mito_data)
samps <- NULL
for (i in seq_along(files)) {
  mito.list <- readRDS(files[i])
  ifelse(exists("mito_data"), mito_data <- rbind(mito_data, as.data.frame(mito.list)), mito_data <- as.data.frame(mito.list))
}
rownames(mito_data) <- NULL
mito_data$Sample <- sapply(strsplit(as.character(mito_data$Sample), split = "_", fixed = TRUE), `[`, 1)
mito_data <- mito_data[match(sample_table, mito_data$Sample),]
#------------------------------------------------------------------------------------------------


#-------------------------------------------------------
# Adding to mito_data the total number of fragments in
# each sample from Stephen's GAlignmentPairs objects
#---------------------------------------------------------------------------------------------
fragDir <- "/dcl01/scharpf1/data/dbruhm/stephen/mito_realignment/total_fragments"
fragFiles <- list.files(fragDir)
mito_data$total_fragments <- NA
for (i in 1:nrow(mito_data)) {
  hit <- grep(paste0(mito_data$Sample[i], "_"), fragFiles)
  stopifnot(length(hit) == 1)
  mito_data$total_fragments[i] <- readRDS(file.path(fragDir, fragFiles[hit]))
}
mito_data$Percent_mito_frag = mito_data$Number_of_Fragments / mito_data$total_fragments * 100
#---------------------------------------------------------------------------------------------


#-------------------------------------
# Saving the mitochondria statistics
#--------------------------------------------------------------------------------------------------
colnames(mito_data) <- c("Sample", "MT Fragments", "MT Reads", "MT Fragment Length (Median)",
                         "MT Fragment Length (Mode)", "MT Fragment Length (Lower Quartile)",
                         "MT Fragment Length (Upper Quartile)", "Total Fragments", "Percent MT Fragments")
mito_data$Sample[which(mito_data$Sample == "PGDX7348P1")] <- "PGDX7348P1_Amended"
write.table(mito_data, file = "/dcl01/scharpf1/data/dbruhm/stephen/mito_realignment/stats/mt_stats_nov12.txt",
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
#--------------------------------------------------------------------------------------------------
