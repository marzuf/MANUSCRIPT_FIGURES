# Rscript test_shaman_plot.R

source("shaman_plot_map.R")


sparse_dt <- read.delim("k562_10kb.txt", header=F, stringsAsFactors=F)
head(sparse_dt)
colnames(sparse_dt) <- c("start1", "start2", "score")
sparse_dt <- na.omit(sparse_dt)

chromo <- "chr7"
start <- 149960000
end <- 150600000
plot_range <- data.frame(chromo=chromo, start=start, end=end)

stopifnot(ncol(plot_range) == 3)

png("shaman_test.png")
shaman_plot_map_score_with_annotations(genome="hg19", points_score=sparse_dt, interval_range=plot_range)
dev.off()

# col.neg.Rdata  col.pos.Rdata  col.scores.Rdata