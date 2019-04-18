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
                               axis.title.x = element_text(face="bold", size=50),
                               axis.title.y = element_blank(),
                               axis.text.y = element_text(size = 20),
                               plot.title = element_text(size=15),
							   legend.position = "none",
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank())

g <- ggplot(df.chr1, aes(x=bin, y=ratio.centered, group=sample)) +
    geom_line( size=1,  color="gray40")
g <- g + labs(x="", y="Fragmentation profile\n", color="")
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
p <- p +  scale_x_continuous(expand = c(0, 0)) + xlab("Chromosome 1")
p <- p + scale_color_manual(values=c('gray50','firebrick'))  + mytheme

pg <- plot_grid(g, p, nrow=2, rel_heights=c(1, 1.2), align="v")
save_plot("~/fig/tmp.pdf", pg, ncol=1, nrow=2, base_height=8, base_width=45)
