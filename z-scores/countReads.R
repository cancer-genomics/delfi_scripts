library(PlasmaTools)

#----------------------------
# Get command line arguments
#----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
bamdir <- args[2]
binGR <- args[3]
outDir <- args[4]
#----------------------------------------------------


#---------------------------
# Import reads into GRanges
#------------------------------------------------------------------------------
data("filters.hg19", package = "PlasmaTools")
bins <- readRDS(binGR)
bamfile <- file.path(bamdir, id)
reads_dir <- file.path(outDir, "gal")
if (!dir.exists(reads_dir)){
    dir.create(reads_dir)
}
sample <- gsub(".bam", "", id)
reads_file <- file.path(reads_dir, paste0(sample, ".rds"))
if(!file.exists(reads_file)){
    original_reads <- getReads(bamfile)
    saveRDS(original_reads, reads_file)
}
#------------------------------------------------------------------------------


#--------------
# Filter reads
#-----------------------------------------------------------
filtered_reads <- filterReads(original_reads, filters.hg19)
#-----------------------------------------------------------


#------------------------
# Count and GC-normalize
#-------------------------------------------------------------
normalized_counts <- countAndNormalize(filtered_reads, bins)
counts_dir <- file.path(outDir, "binCounts")

if (!dir.exists(counts_dir)){
    dir.create(counts_dir)
}

counts_file <- file.path(counts_dir, paste0(sample, ".rds"))
saveRDS(normalized_counts, counts_file)
#-------------------------------------------------------------
