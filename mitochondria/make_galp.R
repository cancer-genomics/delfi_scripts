args <- commandArgs(trailingOnly = TRUE)
bamPath <- args[1]
outDir <- args[2]
sampleName <- gsub(".bam$", "", basename(bamPath))

library(GenomicAlignments)
library(ggplot2)

galp <- readGAlignmentPairs(bamPath, param=ScanBamParam(what="mapq"))

# Keeping reads where read 1 mapq and read2 mapq are >= 30
first.pass <- which(mcols(first(galp))$mapq >= 30)
second.pass <- which(mcols(second(galp))$mapq >= 30)
score.pass <- intersect(first.pass, second.pass)

if (length(score.pass) == 0) {
  stop(paste0("There are no MT reads with Phred>=30 in ", sampleName))
}

galp <- galp[score.pass]

# Keeping reads where both matches align to chrM
first.m <- which(seqnames(first(galp)) == "chrM")
second.m <- which(seqnames(second(galp)) == "chrM")
m.hits <- intersect(first.m, second.m)

if (length(m.hits) == 0) {
  stop(paste0("There are no mitochondrial reads in ", sampleName))
}

galp <- galp[m.hits]

frag <- granges(galp)
frag.len <- width(frag)

# Calculating the total number of MT fragments that pass QC
nFrag <- length(frag)

# Excluding fragments that span the position 16571 in the MT genome for plotting
# Their width is overestimated
hits <- which(frag.len > 15000)
if (length(hits) > 0) {
  frag.len <- frag.len[-hits]
}


p <- ggplot() +
       geom_density(aes(frag.len)) +
       theme_bw() +
       ggtitle(paste0(sampleName, "\n", "Fragments: ", nFrag)) +
       theme(plot.title = element_text(face = "italic")) +
       scale_x_continuous(breaks = seq(0,500,50)) +
       xlim(c(0,500))

setwd(file.path(outDir, "plots"))
pdf(file.path(outDir, "plots", paste0(sampleName, "_mito_density.pdf")), height = 5, width = 8)
p
dev.off()


# Calculating the number of MT reads total
mt.tot <- length(first.m) + length(second.m)

# Calculating the median, IQR, and mode fragment lengths
mt.med <- median(frag.len)
mt.lower <- quantile(frag.len)["25%"]
mt.upper <- quantile(frag.len)["75%"]
dens <- density(frag.len)
mt.mode <- dens$x[which.max(dens$y)]

out.list <- list(Sample = sampleName,
                 Number_of_Fragments = nFrag,
                 Number_of_Reads = mt.tot,
                 Fragment_Length_Median = mt.med,
                 Fragment_Length_Mode = mt.mode,
                 Fragment_Length_Lower_Quartile = mt.lower,
                 Fragment_Length_Upper_Quartile = mt.upper)

saveRDS(out.list, file = file.path(outDir, "stats", paste0(sampleName, "_reads_mapped_to_mito.rds")))
