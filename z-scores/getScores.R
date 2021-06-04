library(Rsamtools)
library(GenomicAlignments)
library(GenomicRanges)
library(matrixStats)
library(GenomeInfoDb)
library(data.table)
library(dplyr)
library(biovizBase)
library(fs)

#--------
# Input
#--------------------------------------------------------------------------------------------------
dataDir <- "countReads.R"
##
## User must specify a plain text file providing the names of the healthy control files created by countReads.R
##
healthyList <- "healthy-list.txt"
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

getArmMeans <- function(bincounts.list) {
    arm_means.list <- lapply(1:length(bincounts.list), function(i){
        b <- bincounts.list[[i]]
        df <- data.table(arm = as.factor(b$arm), adjusted = b$adjusted)
        df$adjusted <- df$adjusted - median(df$adjusted, na.rm = TRUE)
        df <- df %>%
            mutate(total = sum(adjusted)) %>%
            group_by(arm) %>%
            summarize(armmean = mean(adjusted))
        df$armmean <- df$armmean - mean(df$armmean, na.rm = TRUE)
        return(df$armmean)
    })
    return(arm_means.list)
}

getZscores <- function(bincounts.list, normal.indices, loo = FALSE) {
    gr.list <- getArmMeans(bincounts.list)
    zscores.list <- lapply(1:length(gr.list), function(i){
        if (loo == FALSE | !(i %in% normal.indices)) {
            keep <- normal.indices
        } else {
            keep <- setdiff(normal.indices, i)
        }
        normals <- gr.list[keep]
        normals <- do.call("cbind", normals)
        normals.mean <- rowMeans(as.matrix(normals))
        normals.sd <- rowSds(as.matrix(normals))
        zscores <- lapply(gr.list, function(gr){
            (gr - normals.mean) / normals.sd
        })
        return(zscores[[i]])
    })
    return(zscores.list)
}

zscores.loo <- getZscores(bincounts, healthy.index, loo = TRUE)
sample_names <- basename(binfiles)
names(zscores.loo) <- gsub(".rds$", "", sample_names)
saveRDS(zscores.loo, zscore.file)
#--------------------------------------------------------------------
