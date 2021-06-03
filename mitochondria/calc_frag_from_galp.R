library(GenomicAlignments)
library(ggplot2)
library(magrittr)

args <- commandArgs(trailingOnly = TRUE)
bamFile <- args[1]

sampName <- gsub(".bam", "", basename(bamFile))

f.path <- file.path("<PATH TO GALIGNMENT PAIRS DATA>", paste0(sampName, ".rds"))

gr <- readRDS(f.path) %>% granges(on.discordant.seqnames = "drop")
out.vec <- length(gr)

outDir <- "./"

saveRDS(out.vec, file.path(outDir, paste0(sampName, ".rds")))
