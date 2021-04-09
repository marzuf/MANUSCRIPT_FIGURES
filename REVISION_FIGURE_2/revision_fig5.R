cell_line <- "RWPE1"
tad_id <- "chr12_CD194"
tad_chromo <- "chr12"
tad_start <- 54160000
tad_end <- 54440000

# Rscript revision_fig5.R

plotType <- "svg"

myHeight <- 7
myWidth <- 7

plot_around <- 0
resol <- 10000

outFolder <- file.path("REVISION_FIG5_HICPLOT")
dir.create(outFolder)

htc_mat <- read.delim(file.path("extract_hic", paste0(cell_line, "_", tad_chromo, "_obs_KR_", resol/1000, "kb_matrix.txt")), header=F, stringsAsFactors = FALSE)
dim(htc_mat)

plotTit <- paste0(cell_line, " - ", tad_chromo, ":", tad_start, "-", tad_end, " (", tad_id, ")")
subTit <- paste0("Hi-C data - ", resol/1000, "bb (KR - obs)")

source("my_plot_matrix.R")

outFile <- file.path(outFolder, paste0(cell_line, "_", tad_id, "_hicplot.", plotType))

do.call(plotType, list(outFile, height=myHeight, width=myWidth))

my_plot_matrix(mat = htc_mat, 
               tad_coord = c(tad_start,tad_end),
               resolution = resol, 
               bins_around=c(plot_around, plot_around),
               main=plotTit) 
mtext(side=3, text=subTit)

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))