# Rscript shaman_brick_plot.R

source("shaman_plot_map.R")
source("adapt_hicbrick.R")

cat("... load data\n")

sparse_dt <- read.delim("k562_10kb.txt", header=F, stringsAsFactors=F)
head(sparse_dt)
colnames(sparse_dt) <- c("start1", "start2", "score")
sparse_dt <- na.omit(sparse_dt)

cat("... convert data for brick\n")

resol <- 10000
distPlot <- 60
qt_cap <- 0.99

brick_dt <- prep_dt(init_dt=sparse_dt, binResol=resol, plotDist=distPlot, cap_val = qt_cap) 

chromo <- "chr7"
start <- 149960000
end <- 150600000
plot_range <- data.frame(chromo=chromo, start=start, end=end)

my_tads <- data.frame(chr="chr7", start=c(149960001,150120001,150440001), end=c(150120000,150440000,150600000), stringsAsFactors = FALSE)

stopifnot(ncol(plot_range) == 3)

cat("... plot map\n")

map_score <- mz_Brick_vizart_plot_heatmap(File = file.path(tempdir(),
                                                           "shaman_brick_mapOnly.pdf"),
                                          Matrix.df = brick_dt,
                                          tad_dt = my_tads,
                                          # x_coords = "chr3R:1:2900000",
                                          # y_coords = "chr3R:1:2900000",
                                          x_coords = "chr7:149960001:150600000",
                                          y_coords = "chr7:149960001:150600000",
                                          resolution = resol,
                                          colours = "#230C0F",
                                          FUN = Failsafe_log10,
                                          value_cap = qt_cap,
                                          distance = distPlot,
                                          legend_title = "Log10 Hi-C signal",
                                          palette = "Reds",
                                          width = 15,
                                          height = 5,
                                          line_width = 0.8,
                                          cut_corners = TRUE,
                                          rotate = TRUE,
                                          return_object=TRUE)
  
cat("... plot annotation\n")

png("shaman_brick_mapAndTrack.png")
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(5, 1)))
print(map_score, vp=grid::viewport(layout.pos.row=1:3, layout.pos.col=1))
grid::pushViewport(grid::viewport(layout.pos.row=4:5, layout.pos.col=1))
mzAnnotOnly_shaman_plot_map_score_with_annotations(genome="hg19", points_score=sparse_dt, interval_range=plot_range)


dev.off()

png("shaman_brick_trackOnly.png")
mzAnnotOnly_shaman_plot_map_score_with_annotations(genome="hg19", points_score=sparse_dt, interval_range=plot_range)

dev.off()

# col.neg.Rdata  col.pos.Rdata  col.scores.Rdata