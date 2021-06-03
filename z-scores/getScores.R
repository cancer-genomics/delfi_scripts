library(PlasmaTools)
library(ggplot2)
library(fs)

#--------
# Input
#--------------------------------------------------------------------------------------------------
dataDir <- "countReads.R"
##
## User must specify a plain text file providing the names of the non-cancer files created by countReads.R
##
healthyList <- "healthy-list.txt"
binSize="50kb"
#--------------------------------------------------------------------------------------------------

binpath <- file.path(dataDir, "binCounts")
binfiles <- list.files(binpath, full.names = TRUE)
bincounts <- lapply(binfiles, function(x) readRDS(x))

healthy.files <- readLines(healthyList)
healthy.index <- which(basename(binfiles) %in% healthy.files)

#--------------------
# Calculate Z scores
#--------------------------------------------------------------------
if (!dir.exists(file.path(dataDir, "scores"))) {
  dir.create(file.path(dataDir, "scores"), showWarnings = FALSE)
}

zscore.file <- file.path(dataDir, "scores", "zscores.rds")
zscores.loo <- getZscores(bincounts, healthy.index, loo = TRUE)
sample_names <- basename(binfiles)
names(zscores.loo) <- gsub(".rds$", "", sample_names)
saveRDS(zscores.loo, zscore.file)
