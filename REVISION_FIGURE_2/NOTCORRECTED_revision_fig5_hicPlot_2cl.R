stop("-not updated with correct func\n")
# RWPE1	PRAD normal vs. tumor	chr12_TAD194	54160001	54440000	HOXC4,HOXC5,HOXC6
# RWPE1	PRAD normal vs. tumor	chr7_TAD424	116080001	116320000	CAV1,CAV2,MET
# RWPE1	PRAD normal vs. tumor	chr17_TAD174	46720001	46880000	HOXB-AS5,HOXB13,PRAC


# Rscript revision_fig5.R
# Rscript revision_fig5_hicPlot_2cl.R RWPE1 chr12_CD194 chr12 54160000 54440000 0 20000 png 22Rv1
# Rscript revision_fig5_hicPlot_2cl.R RWPE1 chr12_CD194 chr12 54160000 54440000 0 20000 svg
# Rscript revision_fig5_hicPlot_2cl.R RWPE1 chr7_CD424 chr7 116080001 116320000 0 20000 png
# Rscript revision_fig5_hicPlot_2cl.R RWPE1 chr7_CD424 chr7 116080001 116320000 0 20000 svg
# Rscript revision_fig5_hicPlot_2cl.R RWPE1 chr17_CD174 chr17 46720001 46880000 0 20000 png 22Rv1
# Rscript revision_fig5_hicPlot_2cl.R RWPE1 chr17_CD174 chr17 46720001 46880000 0 20000 svg

require(grid)
require(gridBase)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) > 0) {
  stopifnot(length(args) == 9)
  cell_line1 <- args[1]
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
  cell_line2 <- args[9]
  
} else {
  cell_line1 <- "RWPE1"
  tad_id <- "chr12_CD194"
  tad_chromo <- "chr12"
  tad_start <- 54160000
  tad_end <- 54440000
  plot_around <- 0
  resol <- 20000
  plotType <- "png"
  cell_line2 <- "22Rv1"
}

if(plotType=="png") {
  myHeight <- 400
  myWidth <- 400
  
}else{
  myHeight <- 6
  myWidth <- 6
  
}

myHeight <- myHeight*1.2
myWidth <- myWidth*2

outFolder <- file.path("REVISION_FIG5_HICPLOT_2CL")
dir.create(outFolder)

htc_mat1 <- read.delim(file.path("extract_hic", paste0(cell_line1, "_", tad_chromo, "_obs_KR_", resol/1000, "kb_matrix.txt")), header=F, stringsAsFactors = FALSE)
dim(htc_mat1)
colnames(htc_mat1) <- seq(from=0, by=resol, length.out=ncol(htc_mat1))
rownames(htc_mat1) <- seq(from=0, by=resol, length.out=nrow(htc_mat1))

htc_mat2 <- read.delim(file.path("extract_hic", paste0(cell_line2, "_", tad_chromo, "_obs_KR_", resol/1000, "kb_matrix.txt")), header=F, stringsAsFactors = FALSE)
dim(htc_mat2)
colnames(htc_mat2) <- seq(from=0, by=resol, length.out=ncol(htc_mat2))
rownames(htc_mat2) <- seq(from=0, by=resol, length.out=nrow(htc_mat2))

plotTit1 <- paste0(cell_line1, " - ", resol/1000, "kb")
plotTit2 <- paste0(cell_line2, " - ", resol/1000, "kb")

mainTit <- paste0(tad_id, " (", tad_start, "-", tad_end, ")")

source("my_plot_matrix.R")

outFile <- file.path(outFolder, paste0(cell_line1, "_", cell_line2, "_", tad_id, "_", resol/1000, "kb_hicplot.", plotType))


do.call(plotType, list(outFile, height=myHeight, width=myWidth))

opar <- par(no.readonly=TRUE, mar = c(2.5, 3.1, 1, 2))
# CHANGE
# Start new page with traditional graphics
plot.new()

# pushViewport(viewport(x=0,y=0,width=1,height=1,
#                       just=c(0,0), name='base'))
# pushViewport (viewport(x=0,y=0.5,width=1,height=0.5,just=c(0,0)))
# # CHANGE
# par (new=TRUE, fig=gridFIG())

pushViewport(viewport(x=0,y=0,width=1,height=1,
                      just=c(0,0), name='base'))
# pushViewport (viewport(x=0,y=0.5,width=1,height=0.5,just=c(0,0)))
# create left viewport
pushViewport (viewport(x=0,y=0,width=0.5,height=1,just=c(0,0)))
# CHANGE
par (new=TRUE, fig=gridFIG())

# plot in top half of page
my_plot_matrix(mat = htc_mat2, 
               tad_coord = c(tad_start,tad_end),
               resolution = resol, 
               bins_around=c(plot_around, plot_around),
               # main=plotTit,
               main=plotTit1,
               plotOnly = "lower_tri"
) 

# get toplevel view
seekViewport('base')

# # create lower viewport
# pushViewport(viewport(x=0,y=0,width=1, height=0.5, just=c(0,0)))
# create right viewport
pushViewport(viewport(x=0.5,y=0,width=0.5, height=1, just=c(0,0)))


par(new=TRUE,fig=gridFIG())

my_plot_matrix(mat = htc_mat2, 
               tad_coord = c(tad_start,tad_end),
               resolution = resol, 
               bins_around=c(plot_around, plot_around),
               main=plotTit2,
               plotOnly = "upper_tri"
) 
seekViewport("base")
grid.text(paste0(mainTit), vjust=1,x=0.5,y=1,gp = gpar(fontsize = 20))


foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))
