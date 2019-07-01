library(tidyverse)
library(GenomicRanges)
library(cowplot)

AB.chr <- readRDS("AB_chr1.rds")
df.chr1 <- readRDS("chr1_tibble.rds")

x <- 1:length(AB.chr)
m1 <- smooth.spline(x, AB.chr$nsomediff.plasma - median(AB.chr$nsomediff.plasma))
p1 <- predict(m1)

m2 <- smooth.spline(x, AB.chr$nsomediff.lymph - median(AB.chr$nsomediff.lymph))
p2 <- predict(m2)

tissue <- c(Healthy = "Healthy (n=30)",
            Lung = "Lung (n=8)",
            Lymphocyte = "Lymphocyte (n=2)")

mytheme <- theme_classic(base_size=12) + theme( axis.text.x = element_blank(),
                               axis.ticks.x=element_blank(),
                               strip.text.y = element_blank(),
                               strip.text.x = element_text(face="bold", size=20),
                               axis.title.x = element_text(face="bold", size=50),
                               axis.title.y = element_text(face="bold", size=50),
                               axis.text.y = element_text(size = 20),
                               plot.title = element_text(size=15),
							   legend.position = "none",
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank())

g <- ggplot(df.chr1, aes(x=bin, y=ratio.centered, group=sample)) +
    geom_line( size=1,  color="gray40")
g <- g + labs(x="", y="Fragmentation profile", color="")
g <- g + facet_grid(type~., switch="x",space="free_x", scales="free",
                    labeller=labeller(type=tissue))
g <- g + coord_cartesian(xlim = NULL, ylim=c(-0.07, 0.12), expand = TRUE)
g <- g +  scale_x_continuous(expand = c(0, 0))
g <- g +  mytheme
# ggsave("~/fig/tmp.pdf", g, width=40, height=10)


AB.chr$plasma.smooth <- p1$y/sd(p1$y)
AB.chr$lymph.smooth <- p2$y/sd(p2$y)
AB.chr$bin <- 1:length(AB.chr)
AB.long <- as_tibble(AB.chr) %>% gather(`Measurement type`, `Measurement`,
                                         c(eigen, plasma.smooth, lymph.smooth)) %>%
           mutate(color=ifelse(Measurement < 0, "gray50", "brickred"))
AB.long$"Measurement type" <- factor(AB.long$"Measurement type", levels=c("plasma.smooth", "lymph.smooth", "eigen"))

p <- ggplot(AB.long, aes(x=bin, y=Measurement, color=color)) + geom_bar(stat="Identity")
p <- p + facet_grid(`Measurement type`~., space="free", scales="free")
p <- p +  scale_x_continuous(expand = c(0, 0)) + xlab("Chromosome 1") + ylab("Nucleosome distances\nof HI-C profile")
p <- p + scale_color_manual(values=c('gray50','firebrick'))  + mytheme

pg <- plot_grid(g, p, nrow=2, rel_heights=c(1, 1.2), align="v")
save_plot("~/fig/tmp.pdf", pg, ncol=1, nrow=2, base_height=8, base_width=45) 

df.fr <- readRDS("downsampled_5mb.rds")

df.fr %>% group_by(type) %>% filter(bin==1) %>% summarize(n=n())
df.fr$type <- relevel(factor(df.fr$type), "Healthy")
tissue <- c(Healthy = "Healthy (n=30)",
            Lung = "Lung (n=8)")

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

gg <- ggplot(df.fr, aes(x=bin, y=ratio.centered, group=sample)) + geom_line(color="gray50", size=1.25)
# g <- g + stat_summary(aes(group = type), geom = "line", fun.y = median, size = 1.5, color="steelblue")
gg <- gg + labs(x="", y="Fragmentation profile", color="")
gg <- gg + facet_grid(type~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(type=tissue,
                                                                                            arm=arm.labels))
gg <- gg + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
gg <- gg  + mytheme

pg <- plot_grid(gg, g, p, nrow=3, rel_heights=c(1, 1, 1.2), align="v",
                label_size=40, hjust=10, labels=c("a", "b", "c"))

pg <- plot_grid(gg,  nrow=1, label_size=40)
save_plot("~/fig/tmp.pdf", pg, ncol=1, nrow=1, base_height=18, base_width=35)
