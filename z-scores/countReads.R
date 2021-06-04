library(Rsamtools)
library(GenomicAlignments)
library(GenomicRanges)
library(matrixStats)
library(GenomeInfoDb)
library(data.table)
library(dplyr)
library(biovizBase)

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
bins <- readRDS(binGR)
bamfile <- file.path(bamdir, id)
reads_dir <- file.path(outDir, "gal")
if (!dir.exists(reads_dir)){
    dir.create(reads_dir)
}
sample <- gsub(".bam", "", id)
reads_file <- file.path(reads_dir, paste0(sample, ".rds"))

getReads <- function(bamfile, bamparams = NULL) {
    if (is.null(bamparams)) {
        bamparams <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE,
                                                 isSecondaryAlignment = FALSE,
                                                 isUnmappedQuery = FALSE),
                              mapqFilter = 30)
    }
    gal <- readGAlignments(bamfile, param = bamparams)
    reads <- granges(gal)
    return(reads)
}

if(!file.exists(reads_file)){
    original_reads <- getReads(bamfile)
    saveRDS(original_reads, reads_file)
}
#------------------------------------------------------------------------------


#--------------
# Filter reads
#-----------------------------------------------------------
load(<PATH TO filters.hg19.rda>)

filterReads <- function(reads, filtered.regions) {
    reads <- reads[!overlapsAny(reads, filtered.regions)]
    return(reads)
}

filtered_reads <- filterReads(original_reads, filters.hg19)
#-----------------------------------------------------------


#------------------------
# Count and GC-normalize
#-------------------------------------------------------------
gcCorrect <- function(granges) {
    lower <- 0
    upper <- quantile(granges$counts.mult, .99)
    trend.points <- granges$counts.mult > lower & granges$counts.mult < upper
    trend.counts <- granges[trend.points]$counts.mult
    trend.gc <- granges[trend.points]$gc
    num.sample.points <- min(length(trend.counts), 10E3L)
    samp <- sample(1:length(trend.counts), num.sample.points)
    # pad counts and GC
    include <- c(which(bins$gc==min(bins$gc)), which(bins$gc==max(bins$gc)))
    trend.counts <- c(granges[include]$counts.mult, trend.counts[samp])
    trend.gc <- c(granges[include]$gc, trend.gc[samp])
    initial.trend <- loess(trend.counts ~ trend.gc)
    i <- seq(min(bins$gc), max(bins$gc), by = 0.001)
    final.trend <- loess(predict(initial.trend, i) ~ i)
    counts.pred <- predict(final.trend, granges$gc)
    return(counts.pred)
}

countAndNormalize <- function(reads, bins) {
    bins$counts <- countOverlaps(bins, reads)
    ### normalize for filtered bases
    ### weight is to get count per expected width of basepairs
    weight <- (width(bins)[1]) / (width(bins) - bins$filtered.bases)
    bins$counts.mult <- bins$counts * weight
    bins$counts.mult <- log2(bins$counts.mult + 1)
    bins$loess.pred <- gcCorrect(bins)
    bins$adjusted <- bins$counts.mult - bins$loess.pred
    return(bins)
}

normalized_counts <- countAndNormalize(filtered_reads, bins)
counts_dir <- file.path(outDir, "binCounts")

if (!dir.exists(counts_dir)){
    dir.create(counts_dir)
}

counts_file <- file.path(counts_dir, paste0(sample, ".rds"))
saveRDS(normalized_counts, counts_file)
#-------------------------------------------------------------
