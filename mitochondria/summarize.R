#--------------------------------------------------------------------------
# This script generates the mitochondrial representation statistic called 
# 'Percent MT Fragments; as well as other statistics relateed to
# reads mapping to the mitochondrial genome / mitochondrial fragment lengths.
#--------------------------------------------------------------------------

library(readr)
library(ggplot2)


#----------------------------------------------------
# Compiling the mitochondria stats from the Bowtie2
# realignment into a table
#------------------------------------------------------------------------------------------------
files <- list.files(<PATH THAT WAS ASSIGNED BY USER TO THE outDir VARIABLE IN mito_statistics.sh>,
                     pattern = "reads_mapped_to_mito.rds", full.names = TRUE, recursive = TRUE)

rm(mito_data)
samps <- NULL
for (i in seq_along(files)) {
  mito.list <- readRDS(files[i])
  ifelse(exists("mito_data"), mito_data <- rbind(mito_data, as.data.frame(mito.list)), mito_data <- as.data.frame(mito.list))
}
rownames(mito_data) <- NULL
mito_data$Sample <- sapply(strsplit(as.character(mito_data$Sample), split = "_", fixed = TRUE), `[`, 1)
#------------------------------------------------------------------------------------------------


#-------------------------------------------------------
# Adding to mito_data the total number of fragments in
# each sample
#---------------------------------------------------------------------------------------------
fragDir <- <THE FILE PATH ASSIGNED TO THE outDir VARIABLE IN calc_frag_from_galp.R>
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
write.table(mito_data, file = <A FILE PATH TO WRITE THE OUTPUT TABLE TO>,
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
#--------------------------------------------------------------------------------------------------
