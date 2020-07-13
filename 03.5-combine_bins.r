library(tidyverse)
library(GenomicRanges)
bindir <- "../bins_100kb"

metadata <- read_csv("sample_reference.csv")
ids <- metadata %>% select(`WGS ID`) %>% unlist()
files <- file.path(bindir, paste0(ids, "_bin_100kb.rds"))

bins.list <- lapply(files, readRDS)
tib.list <- lapply(bins.list, as_tibble)
names(tib.list) <- ids
tib.list <- map2(tib.list, names(tib.list), ~ mutate(.x, id = .y)) %>%
    bind_rows() %>% select(id, everything())

tib.list <- tib.list %>% select(-matches("X"))
saveRDS(tib.list, "bins_100kbcompartments.rds")
