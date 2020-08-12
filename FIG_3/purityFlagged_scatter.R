
options(scipen=100)


# Rscript purityFlagged_scatter.R


script_name <- "purityFlagged_scatter.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(ggplot2)
require(reshape2)
require(ggpubr)
require(ggrepel)


registerDoMC(4)

source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")


outFolder <- "PURITYFLAGGED_SCATTER"
dir.create(outFolder, recursive = TRUE)

all_cols[all_cols == "red"] <- "brown3"
all_cols[all_cols == "blue"] <- "darkblue"
all_cols[all_cols == "green"] <- "forestgreen"

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

source("../settings.R")
source("../full_dataset_names.R")

myWidthGG <- 5
myHeightGG <- 5

setDir <- "/media/electron"
setDir <- ""

mainFolder <- file.path(runFolder)
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(pipFolder)
stopifnot(dir.exists(pipFolder))

purity_ds <- "aran"
pm <- "CPE"
transfExpr <- "log10"

tadSignifThresh <- 0.01


inFolder <- file.path(runFolder, "CMP_SIGNIF_PURITYFLAGGED_FINAL_RANDOM", purity_ds, pm, transfExpr)
in_dt <- get(load(file.path(inFolder, "out_dt.Rdata")))

all_rd_types <- c("OBSERVED", unique(in_dt$data_type)) 
rd_type = all_rd_types[3]

# obs_dt$ratioFlagged <- obs_dt$nPurityFlagged/obs_dt$nTot  # computed in CMP_SIGNIF_PURITYFLAGGED_FINAL_RANDOM
# all_dt$ratioSignifFlagged <- all_dt$nSignifAndFlagged/all_dt$nSignif # computed in signif_purityflagged_final



for(rd_type in all_rd_types) {
  
  xvar <- ifelse(rd_type =="OBSERVED", "ratioSignifFlagged_obs", "ratioSignifFlagged_rd")
  yvar <- ifelse(rd_type =="OBSERVED", "ratioFlagged_obs", "ratioFlagged_rd")
  stopifnot(xvar %in% colnames(in_dt))
  stopifnot(yvar %in% colnames(in_dt))
  
  if(rd_type != "OBSERVED") {
    plot_dt <-   in_dt[in_dt$data_type == rd_type,]
  } else {
    plot_dt <- in_dt
  }
  
  plot_dt <- plot_dt[,c("dataset", xvar, yvar)]
  plot_dt <- na.omit(unique(plot_dt))
  if(nrow(plot_dt) == 0) next
  stopifnot(!duplicated(plot_dt$dataset))
  plot_dt$dotcols <- all_cols[all_cmps[basename(as.character(plot_dt$dataset))]]
  stopifnot(!is.na(plot_dt$dotcols))
  
  my_x <- plot_dt[,c(xvar)]
  my_y <- plot_dt[,c(yvar)]
  
  outFile <- file.path(outFolder, paste0(yvar, "_", xvar, "_", rd_type, ".", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  par(bty="L")
  plot(
    x = my_x*100,
    y = my_y*100,
    col = plot_dt$dotcols,
    xlab = "% of purity-flagged signif. TADs",
    ylab = "% of purity-flagged TADs",
    pch=16,
    cex=0.7,
    main = paste0("Ratio of purity-flagged TADs"),
    cex.axis=plotCex,
    cex.main = plotCex,
    cex.lab = plotCex
  )
  mtext(side=3, text = paste0("# datsets = ", nrow(plot_dt), "; TAD signif. tresh: adj. p-val <= ", tadSignifThresh), cex  =plotCex, font=3)
  addCorr(x=my_x,y=my_y, bty="n")
  legend("bottomright", legend=cmp_names[names(all_cols)], pch=16,col = all_cols, bty="n")
  # curve(1*x, col="darkgrey", add=T)
  #legend("bottomleft", 
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}








