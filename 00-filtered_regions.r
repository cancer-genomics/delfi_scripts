library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(tidyverse)
library(httr)

### hg19 gaps & blacklisted regions
genome <- "hg19"
mySession <- browserSession()
genome(mySession) <- genome
gaps <- getTable(ucscTableQuery(mySession, track="gap"))
gaps.hg19 <- GRanges(gaps$chrom, IRanges(gaps$chromStart,
                                     gaps$chromEnd),
                 type=gaps$type)
gaps.hg19 <- keepSeqlevels(gaps.hg19, paste0("chr", c(1:22, "X", "Y")),
                           pruning.mode="coarse")
hsapiens <- BSgenome.Hsapiens.UCSC.hg19::Hsapiens
seqinfo(gaps.hg19) <- seqinfo(hsapiens)[seqlevels(gaps.hg19),]
save(gaps.hg19, file="gaps.hg19.rda")
# devtools::use_data(gaps.hg19, overwrite = TRUE)

blacklisted.file <- httr::content(GET("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz"))
blacklisted.tib <- read_tsv(gzcon(rawConnection(blacklisted.file)),
                            col_names=c("seqnames", "start",
                                        "end", "name", "score"))
blacklisted.tib <- blacklisted.tib %>% mutate(start=start+1)
filters.hg19 <- makeGRangesFromDataFrame(blacklisted.tib,
                                           keep.extra.columns=TRUE)
filters.hg19 <- keepSeqlevels(filters.hg19, paste0("chr", c(1:22, "X", "Y")),
                           pruning.mode="coarse")
seqinfo(filters.hg19) <- seqinfo(Hsapiens)[seqlevels(filters.hg19),]
save(filters.hg19, file="filters.hg19.rda")
# devtools::use_data(filters.hg19, overwrite = TRUE)
