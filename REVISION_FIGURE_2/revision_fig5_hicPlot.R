# RWPE1	PRAD normal vs. tumor	chr12_TAD194	54160001	54440000	HOXC4,HOXC5,HOXC6
# RWPE1	PRAD normal vs. tumor	chr7_TAD424	116080001	116320000	CAV1,CAV2,MET
# RWPE1	PRAD normal vs. tumor	chr17_TAD174	46720001	46880000	HOXB-AS5,HOXB13,PRAC


# Rscript revision_fig5.R
# Rscript revision_fig5_hicPlot.R RWPE1 chr12_CD194 chr12 54160000 54440000 0 20000 png
# Rscript revision_fig5_hicPlot.R RWPE1 chr12_CD194 chr12 54160000 54440000 0 20000 svg
# Rscript revision_fig5_hicPlot.R RWPE1 chr7_CD424 chr7 116080000 116320000 0 20000 png
# Rscript revision_fig5_hicPlot.R RWPE1 chr7_CD424 chr7 116080000 116320000 0 20000 svg
# Rscript revision_fig5_hicPlot.R RWPE1 chr17_CD174 chr17 46720000 46880000 0 20000 png
# Rscript revision_fig5_hicPlot.R RWPE1 chr17_CD174 chr17 46720000 46880000 0 20000 svg



args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0) {
  stopifnot(length(args) == 8)
  cell_line <- args[1]
  tad_id <- args[2]
  tad_chromo <- args[3]
  tad_start <- as.numeric(args[4])
  stopifnot(!is.na(tad_chromo))
  tad_end <-  as.numeric(args[5])
  stopifnot(!is.na(tad_end))
  plot_around <-  as.numeric(args[6])
  stopifnot(!is.na(plot_around))
  resol <- as.numeric(args[7])
  stopifnot(!is.na(resol))
  plotType <- args[8]
  
} else {
  cell_line <- "RWPE1"
  tad_id <- "chr12_CD194"
  tad_chromo <- "chr12"
  tad_start <- 54160001
  tad_end <- 54440000
  plot_around <- 0
  resol <- 20000
  plotType <- "png"
}

if(plotType=="png") {
  myHeight <- 400
  myWidth <- 400
  
}else{
  myHeight <- 6
  myWidth <- 6
  
}

stopifnot(tad_end > tad_start)

stopifnot(tad_end %% resol == 0)
stopifnot( (tad_start-1) %% resol == 0)

outFolder <- file.path("REVISION_FIG5_HICPLOT")
dir.create(outFolder)

htc_mat <- read.delim(file.path("extract_hic", paste0(cell_line, "_", tad_chromo, "_obs_KR_", resol/1000, "kb_matrix.txt")), header=F, stringsAsFactors = FALSE)
dim(htc_mat)

colnames(htc_mat) <- seq(from=0, by=resol, length.out=ncol(htc_mat))
rownames(htc_mat) <- seq(from=0, by=resol, length.out=nrow(htc_mat))

plotTit <- paste0(tad_id, " (", tad_start, "-", tad_end, ")")
subTit <- paste0(cell_line, " Hi-C data - ", resol/1000, "kb (KR - obs)")

source("my_plot_matrix.R")

outFile <- file.path(outFolder, paste0(cell_line, "_", tad_id, "_", resol/1000, "kb_hicplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
# outFile <- file.path(outFolder, paste0(cell_line, "_", tad_id, "_", resol/1000, "kb_hicplot_2.", "pdf"))
# pdf(outFile, height=7, width=7, dpi=10)

my_plot_matrix(mat = htc_mat, 
               tad_coord = c(tad_start,tad_end),
               resolution = resol, 
               bins_around=c(plot_around, plot_around),
               main=plotTit) 
mtext(side=3, text=subTit)

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

