library(tidyverse)
library(multidplyr)
library(GenomicRanges)
library(readxl)

df.fr <- readRDS("../inst/extdata/bins_100kbcompartments.rds")
master <- read_csv("sample_reference.csv")

df.fr2 <- inner_join(df.fr, master, by=c("sample"="WGS ID"))

hic.eigen <- (df.fr2 %>% filter(sample=="PGDX10346P1"))$hic.eigen

armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
df.fr2$arm <- factor(df.fr2$arm, levels=armlevels)

### combine
df.fr2 <- df.fr2 %>% group_by(sample, arm) %>%
    mutate(combine = ifelse(grepl("p", arm), ceiling((1:length(arm))/50),
                           ceiling(rev((1:length(arm))/50) )))

df.fr3 <- df.fr2 %>% group_by(sample, chromosome, arm, combine) %>%
    summarize(short2=sum(short),
              long2=sum(long),
              short.corrected2=sum(short.corrected),
              long.corrected2=sum(long.corrected),
              hic.eigen=mean(hic.eigen),
              gc=mean(gc),
              ratio2=mean(ratio),
              ratio.corrected2=mean(ratio.corrected),
              nfrags2=sum(nfrags),
              nfrags.corrected2=sum(nfrags.corrected),
              domain = median(as.integer(domain)),
              short.var=var(short.corrected),
              long.var=var(long.corrected),
              nfrags.var=var(nfrags.corrected),
              mode_size=unique(mode_size),
              mean_size=unique(mean_size),
              median_size=unique(median_size),
              q25_size=unique(q25_size),
              q75_size=unique(q75_size),
              start=start[1],
              end=rev(end)[1],
              binsize = n())
### assign bins
df.fr3 <- inner_join(df.fr3, master, by=c("sample"="WGS ID"))
df.fr3 <- df.fr3 %>% mutate(type = gsub(" Cancer|carcinoma", "", `Patient Type`, ignore.case=TRUE))

df.fr3 <- df.fr3 %>% filter(binsize==50)
df.fr3 <- df.fr3 %>% group_by(sample) %>% mutate(bin = 1:length(sample))

hic.eigen <- (df.fr3 %>% filter(sample=="PGDX10346P1"))$hic.eigen
saveRDS(df.fr3, "../inst/extdata/bins/bins_5mbcompartments.rds")
