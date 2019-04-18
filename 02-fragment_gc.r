tx <-as.numeric(Sys.getenv("SGE_TASK_ID"))
galpdir <- "/dcl01/scharpf1/data/galignmentpairs/Low_Coverage_WGS"

files <- list.files(galpdir, full.names=TRUE)
file <- files[tx]

names <- strsplit(basename(file), "_")[[1]][1]

outdir <- "/dcl01/scharpf1/data/galignmentpairs/Low_Coverage_WGS/fragments"
# if(file.exists(file.path(outdir, paste0(names, "_frags.rds")))) q('no')

library(GenomicAlignments)
library(GenomicRanges)
library(Rsamtools)
library(devtools)
library(Homo.sapiens)
library(biovizBase)
library(BSgenome.Hsapiens.UCSC.hg19)
class(Homo.sapiens)


galp <- readRDS(file)
frags <- granges(keepSeqlevels(galp, paste0("chr", 1:22), pruning.mode="coarse"),
             on.discordant.seqnames="drop")

## filter outliers
w.all <- width(frags)
q.all <- quantile(w.all, c(0.001, 0.999))
frags <- frags[which(w.all > q.all[1] & w.all < q.all[2])]

gcs <- GCcontent(Hsapiens, unstrand(frags))
frags$gc <- gcs

saveRDS(frags, file.path(outdir, paste0(names, "_frags.rds")) )
q('no')
