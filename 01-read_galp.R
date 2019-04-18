library(GenomicAlignments)
library(GenomicRanges)
library(getopt)

### Used for getting information from shell script
args <- commandArgs(trailingOnly = TRUE)
hh <- paste(unlist(args), collapse = " ")
listoptions <- unlist(strsplit(hh, "--"))[-1]
options.args <- sapply(listoptions, function(x) {
    unlist(strsplit(x, " "))[-1]
})
options.names <- sapply(listoptions, function(x) {
    option <- unlist(strsplit(x, " "))[1]
})
names(options.args) <- unlist(options.names)
id <- options.args[1]
bamdir <- options.args[2]
galpdir <- options.args[3]
###

### Read GAlignmentPairs
bamfile <- file.path(bamdir, id)
indexed.bam <- gsub("$", ".bai", bamfile)
if (!file.exists(indexed.bam)) {
    indexBam(bamfile)
}

param <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE,
                                         isSecondaryAlignment = FALSE,
                                         isUnmappedQuery = FALSE),
                      mapqFilter = 30)
sample <- gsub(".bam", "", id)

galp.file <- file.path(galpdir, paste0(sample, ".rds"))
galp <- readGAlignmentPairs(bamfile, param = param)
saveRDS(galp, galp.file)
