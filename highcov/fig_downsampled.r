library(tidyverse)
library(GenomicRanges)
library(cowplot)

df.fr <- readRDS("downsampled_5mb.rds")

df.fr %>% group_by(type) %>% filter(bin==1) %>% summarize(n=n())
df.fr$type <- relevel(factor(df.fr$type), "Healthy")
tissue <- c(Healthy = "Healthy (n=30)",
            Lung = "Lung (n=8)")

mytheme <- theme_classic(base_size=12) + theme(
                               axis.text.x = element_blank(),
                               axis.ticks.x=element_blank(),
                               strip.text.x = element_text(size = 17),
#                                strip.text.y = element_text(size=23),
                               strip.text.y = element_blank(),
                               axis.title.x = element_text(face="bold", size=17),
                               axis.title.y = element_blank(),
                               axis.text.y = element_text(size = 20),
                               plot.title = element_text(size=15),
                               legend.title = element_text(size=10),
                               legend.text = element_text(size=10),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
							   strip.background=element_rect(fill="white", color="white"))

arm <- df.fr %>% group_by(arm) %>%
  summarize(n=n()) %>%
  mutate(arm = as.character(arm))
small.arms <- setNames(c("", "12q", "", "16q",
                         "", "17q", "", "18q",
                         "", "19", "", "20",
                         "21", "22"),
                       c("12p", "12q", "16p", "16q",
                         "17p", "17q", "18p", "18q",
                         "19p", "19q", "20p", "20q",
                         "21q", "22q"))
arm.labels <- setNames(arm$arm, arm$arm)
arm.labels[names(small.arms)] <- small.arms

g <- ggplot(df.fr, aes(x=bin, y=ratio.centered, group=sample)) + geom_line(color="gray50", size=0.75)
# g <- g + stat_summary(aes(group = type), geom = "line", fun.y = median, size = 1.5, color="steelblue")
g <- g + labs(x="", y="Fragmentation profile\n", color="")
g <- g + facet_grid(type~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(type=tissue,
                                                                                            arm=arm.labels))
g <- g + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
g <- g + mytheme

g1 <- ggplot(df.fr, aes(x=bin, y=ratio.centered25, group=sample)) + geom_line(color="gray50", size=0.75)
# g1 <- g1 + stat_summary(aes(group = type), geom = "line", fun.y = median, size = 1.5, color="steelblue")
g1 <- g1 + labs(x="", y="Fragmentation profile\n", color="")
g1 <- g1 + facet_grid(type~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(type=tissue,
                                                                                              arm=arm.labels))
g1 <- g1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
g1 <- g1 + mytheme

g2 <- ggplot(df.fr, aes(x=bin, y=ratio.centered10, group=sample)) + geom_line(color="gray50", size=0.75)
# g2 <- g2 + stat_summary(aes(group = type), geom = "line", fun.y = median, size = 1.5, color="steelblue")
g2 <- g2 + labs(x="", y="Fragmentation profile\n", color="")
g2 <- g2 + facet_grid(type~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(type=tissue,
                                                                                              arm=arm.labels))
g2 <- g2 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
g2 <- g2 + mytheme

g3 <- ggplot(df.fr, aes(x=bin, y=ratio.centered05, group=sample)) + geom_line(color="gray50", size=0.75)
# g3 <- g3 + stat_summary(aes(group = type), geom = "line", fun.y = median, size = 1.5, color="steelblue")
g3 <- g3 + labs(x="", y="Fragmentation profile\n", color="")
g3 <- g3 + facet_grid(type~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(type=tissue,
                                                                                              arm=arm.labels))
g3 <- g3 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
g3 <- g3 + mytheme

g4 <- ggplot(df.fr, aes(x=bin, y=ratio.centered02, group=sample)) + geom_line(color="gray50", size=0.75)
# g4 <- g4 + stat_summary(aes(group = type), geom = "line", fun.y = median, size = 1.5, color="steelblue")
g4 <- g4 + labs(x="", y="Fragmentation profile\n", color="")
g4 <- g4 + facet_grid(type~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(type=tissue,
                                                                                              arm=arm.labels))
g4 <- g4 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
g4 <- g4 + mytheme

g5 <- ggplot(df.fr, aes(x=bin, y=ratio.centered01, group=sample)) + geom_line(color="gray50", size=0.75)
# g5 <- g5 + stat_summary(aes(group = type), geom = "line", fun.y = median, size = 1.5, color="steelblue")
g5 <- g5 + labs(x="", y="Fragmentation profile\n", color="")
g5 <- g5 + facet_grid(type~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(type=tissue,
                                                                                              arm=arm.labels))
g5 <- g5 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
g5 <- g5 + mytheme

plx <- plot_grid(g, g1,g2,g3,g4,g5, nrow=6, rel_heights=1, align="v")
save_plot("downsampled_profiles.pdf", plx, ncol=1, nrow=6, base_height=8, base_width=40, limitsize=FALSE)

cor.df <- df.fr %>% group_by(sample, type) %>% summarize(
                                               cor100 = cor(ratio, ratio),
                                               cor25 = cor(ratio, ratio25),
                                               cor10 = cor(ratio, ratio10),
                                               cor05 = cor(ratio, ratio05),
                                               cor02 = cor(ratio, ratio02),
                                               cor01 = cor(ratio, ratio01),
                                               sum100 = sum(nfrags),
                                               sum25 = sum(nfrags25),
                                               sum10 = sum(nfrags10),
                                               sum05 = sum(nfrags05),
                                               sum05 = sum(nfrags05),
                                               sum01 = sum(nfrags01))


cor.wide <- cor.df %>% gather(`Downsample`, Correlation, cor100:cor01)
cor.wide$Downsample = factor(cor.wide$Downsample, levels=c("cor100", "cor25",
                                                           "cor10", "cor05",
                                                           "cor02", "cor01"),
                             labels=c("100%", "25%", "10%", "5%", "2%", "1%"))

p <- ggplot(cor.wide, aes(Downsample, Correlation, group=sample)) + geom_point()
p <- p + geom_line() + facet_grid(.~type) + ylim(c(0, 1))
p <- p + xlab("Percent downsampled") + ylab("Correlation to high coverage")
ggsave("downsampled_correlations.pdf", p, width=12)
