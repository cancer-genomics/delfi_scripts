library(tidyverse)
tx <-as.numeric(Sys.getenv("SGE_TASK_ID"))

gc.correct <- function(coverage, bias) {
    i <- seq(min(bias, na.rm=TRUE), max(bias, na.rm=TRUE), by = 0.001)
    coverage.trend <- loess(coverage ~ bias)
    coverage.model <- loess(predict(coverage.trend, i) ~ i)
    coverage.pred <- predict(coverage.model, bias)
    coverage.corrected <- coverage - coverage.pred + median(coverage)
}

fragpath <- "../fragments"
fragfiles <- list.files(fragpath, pattern=".rds",full.name=TRUE)
fragfile <- fragfiles[tx]

id <- strsplit(basename(fragfile), "\\.")[[1]][1]

outdir <- "." ####
filename <- file.path(outdir, paste0(id, "_bin_100kb.rds"))
if(file.exists(filename)) q('no')

library(GenomicRanges)
library(rtracklayer)
library(Homo.sapiens)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Rsamtools)
class(Homo.sapiens)
library(devtools)
library(biovizBase)
load("./filters.hg19.rda")

library(RCurl)
ABurl <- getURL('https://raw.githubusercontent.com/Jfortin1/HiC_AB_Compartments/master/data/hic_compartments_100kb_ebv_2014.txt', ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)

AB <- read.table(textConnection(ABurl), header=TRUE)
AB <- makeGRangesFromDataFrame(AB, keep.extra.columns=TRUE)

chromosomes <- GRanges(paste0("chr", 1:22),
                       IRanges(0, seqlengths(Hsapiens)[1:22]))

tcmeres <- gaps.hg19[grepl("centromere|telomere", gaps.hg19$type)]

arms <- GenomicRanges::setdiff(chromosomes, tcmeres)
arms <- arms[-c(25,27,29,41,43)]

armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")

arms$arm <- armlevels
AB <- AB[-queryHits(findOverlaps(AB, gaps.hg19))]
AB <- AB[queryHits(findOverlaps(AB, arms))]
AB$arm <- armlevels[subjectHits(findOverlaps(AB, arms))]

seqinfo(AB) <- seqinfo(Hsapiens)[seqlevels(seqinfo(AB))]
AB <- trim(AB)
AB$gc <- GCcontent(Hsapiens, AB)

## These bins had no coverage
AB <- AB[-c(8780, 13665)]
fragments <- readRDS(fragfile)
# 
### Filters
fragments <- fragments[-queryHits(findOverlaps(fragments, filters.hg19))]
w.all <- width(fragments)

fragments <- fragments[which(w.all >= 100 & w.all <= 220)]
w <- width(fragments)

frag.list <- split(fragments, w)

counts <- sapply(frag.list, function(x) countOverlaps(AB, x))
if(min(w) > 100) {
    m0 <- matrix(0, ncol=min(w) - 100, nrow=nrow(counts),
                 dimnames=list(rownames(counts), 100:(min(w)-1)))
    counts <- cbind(m0, counts)
}

olaps <- findOverlaps(fragments, AB)
bin.list <- split(fragments[queryHits(olaps)], subjectHits(olaps))
bingc <- rep(NA, length(bin.list))
bingc[unique(subjectHits(olaps))] <- sapply(bin.list, function(x) mean(x$gc))

### Get modes
Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}
modes <- Mode(w)
medians <- median(w)
q25 <- quantile(w, 0.25)
q75 <- quantile(w, 0.75)

short <- rowSums(counts[,1:51])
long <- rowSums(counts[,52:121])
ratio <- short/long
short.corrected=gc.correct(short, bingc)
long.corrected=gc.correct(long, bingc)
nfrags.corrected=gc.correct(short+long, bingc)
ratio.corrected=gc.correct(ratio, bingc)

AB$short <- short
AB$long <- long
AB$short.corrected <- short.corrected
AB$long.corrected <- long.corrected
AB$nfrags.corrected <- nfrags.corrected
AB$ratio.corrected <- ratio.corrected

AB$mode <- modes
AB$mean <- round(mean(w), 2)
AB$median <- medians
AB$quantile.25 <- q25
AB$quantile.75 <- q75
AB$frag.gc <- bingc

for(i in 1:ncol(counts)) elementMetadata(AB)[,colnames(counts)[i]] <- counts[,i]

saveRDS(AB, filename)
q('no')
