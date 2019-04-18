library(tidyverse)
library(readxl)

master <- read_csv("sample_reference.csv")
df.fr3 <- readRDS("../inst/extdata/bins_5mbcompartments.rds")

healthy.median <- df.fr3 %>%
    group_by(bin) %>% 
    summarize(median.cov=median(nfrags2, na.rm=TRUE),
              median.short=median(short2, na.rm=TRUE),
              median.long=median(long2, na.rm=TRUE),
              median.ratio=median(ratio2, na.rm=TRUE),
              median.corrected.cov=median(nfrags.corrected2, na.rm=TRUE),
              median.corrected.short=median(short.corrected2, na.rm=TRUE),
              median.corrected.long=median(long.corrected2, na.rm=TRUE),
              median.corrected.ratio=median(ratio.corrected2, na.rm=TRUE),
              median.corrected.ratio2=median(short.corrected2/long.corrected2, na.rm=TRUE))

summary.df <- df.fr3 %>% ungroup() %>% group_by(sample, type) %>%
    summarize(cov.cor=cor(nfrags2, healthy.median$median.cov, method="pearson", use="complete.obs"),
              short.cor=cor(short2, healthy.median$median.short, method="pearson", use="complete.obs"),
              long.cor=cor(long2, healthy.median$median.long, method="pearson", use="complete.obs"),
              ratio.cor=cor(ratio2, healthy.median$median.ratio, method="pearson", use="complete.obs"),
              cov.corrected.cor=cor(nfrags.corrected2, healthy.median$median.corrected.cov, method="pearson", use="complete.obs"),
              short.corrected.cor=cor(short.corrected2, healthy.median$median.corrected.short, method="pearson", use="complete.obs"),
              long.corrected.cor=cor(long.corrected2, healthy.median$median.corrected.long, method="pearson", use="complete.obs"),
              ratio.corrected.cor=cor(ratio.corrected2, healthy.median$median.corrected.ratio, method="pearson", use="complete.obs"),
              ratio2.corrected.cor=cor(short.corrected2/long.corrected2, healthy.median$median.corrected.ratio2, method="pearson", use="complete.obs"),
              nfrags = sum(nfrags2),
              mode_size=unique(mode_size),
              mean_size=unique(mean_size),
              median_size=unique(median_size),
              q25_size=unique(q25_size),
              q75_size=unique(q75_size),
              hqbases_analyzed = 100*sum(nfrags)*2,
              coverage = hqbases_analyzed/(504*5e6)
              )

summary.df <- inner_join(summary.df, master, by=c("sample"="WGS ID"))
summary.df$`type` = relevel(as.factor(summary.df$`type`), "Healthy")

saveRDS(summary.df, "../inst/extdata/summary_tibble.rds")
