library(tidyverse)
library(caret)
library(pROC)

df.fr3 <- readRDS("../inst/extdata/bins_5mbcompartments.rds")
summary.df <- readRDS("../inst/extdata/summary_tibble.rds")

features.cov <- df.fr3  %>% ungroup() %>%
    select(nfrags.corrected2, sample, bin) %>%
    spread(sample, nfrags.corrected2) %>%
    select(-bin) %>% 
    na.omit() %>%
    scale() %>%
    t() %>%
    as.data.frame()

features.short <- df.fr3  %>% ungroup() %>%
    select(short.corrected2, sample, bin) %>%
    spread(sample, short.corrected2) %>%
    select(-bin) %>% 
    na.omit() %>%
    scale() %>%
    t() %>%
    as.data.frame()

features.sl <- cbind(features.cov, features.short)
colnames(features.sl) <- c(paste0("total", 1:498), paste0("short", 1:498))
features.sl$type <- ifelse(summary.df$type == "Healthy", "Healthy", "Cancer")

features <- cbind(features.sl,
             as.matrix(summary.df %>% ungroup() %>%
                       select(contains("Z Score"))))
features$mito <- -log10(summary.df$"% of Mapped Reads Mapping to Mitochondria")

ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 10,
                     verboseIter = FALSE,
                     savePredictions=TRUE,
                     classProbs=TRUE,
#                      preProcOptions=list(thres = 0.90),
                     summaryFunction = twoClassSummary)
set.seed(1234)
model_gbm <- caret::train(type ~ .,
                               data = features,
                               method = 'gbm',
                               tuneGrid=data.frame(n.trees=150,
                                                   interaction.depth=3,
                                                   shrinkage=0.1,
                                                   n.minobsinnode=10),
                               preProcess = c("corr", "nzv"),
                         trControl = ctrl)

####### Only short/total coverage
set.seed(1234)
model_sl <- caret::train(type ~ .,
                               data = features.sl,
                               method = 'gbm',
                               tuneGrid=data.frame(n.trees=150, interaction.depth=3,
                                                   shrinkage=0.1,
                                                   n.minobsinnode=10),
                               preProcess = c("corr", "nzv"),
                         trControl = ctrl)

###### Only z-scores
features.z <- summary.df %>% ungroup() %>% select(contains("Z Score"))
features.z$type <- ifelse(summary.df$type == "Healthy", "Healthy", "Cancer")

set.seed(1234)
model_z <- caret::train(type ~ .,
                               data = features.z,
                               method = 'gbm',
                               tuneGrid=data.frame(n.trees=150,
                                                   interaction.depth=3,
                                                   shrinkage=0.1,
                                                   n.minobsinnode=10),
                               preProcess = c("corr", "nzv"),
                         trControl = ctrl)

#### Save
models.list <- list("all"=model_gbm, "SL"=model_sl, "z"=model_z)
saveRDS(models.list, "../inst/extdata/models_list.rds")

pred.tbl <- model_gbm$pred %>% filter(n.trees==150, interaction.depth==3) %>%
group_by(rowIndex) %>% dplyr::summarize(obs=obs[1], Cancer=mean(Cancer))
pred.tbl$sample <- rownames(features)
pred.tbl <- inner_join(pred.tbl, summary.df)

## 95% specificity
cutoff <- (pred.tbl %>% filter(type=="Healthy") %>%
           arrange(desc(Cancer)))$Cancer[11]
cutoff98 <- (pred.tbl %>% filter(type=="Healthy") %>%
             arrange(desc(Cancer)))$Cancer[5]
## 90% specificity cutoff to be used in tissue prediction.
cutoff90 <- (pred.tbl %>% filter(type=="Healthy") %>%
             arrange(desc(Cancer)))$Cancer[21]

pred.tbl <- pred.tbl %>%
    mutate(detected95 = ifelse(Cancer > cutoff, "Detected", "Not detected"),
       detected98 = ifelse(Cancer > cutoff98, "Detected", "Not detected"),
       detected90 = ifelse(Cancer > cutoff90, "Detected", "Not detected"),
       stage = gsub("A|B|C", "", `Stage at Diagnosis`))

write.csv(inner_join(summary.df %>% select(-contains("Z Score")), pred.tbl %>%
                     select(rowIndex, sample, stage, Cancer, detected95, detected98),
                     by=c("sample"="sample")),"../inst/extdata/predictions_gbm.csv",
          row.names=FALSE)
