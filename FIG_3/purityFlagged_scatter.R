
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
  
  save(plot_dt, file=file.path(outFolder, paste0("plot_dt_", rd_type, ".Rdata")), version=2)
  
}



nTot_dt <- in_dt[,c("dataset", "nTot_obs", "nPurityFlagged_obs")]
nTot_dt <- unique(nTot_dt)
nTot_dt$nTot_sub <- nTot_dt$nTot_obs - nTot_dt$nPurityFlagged_obs
nTot_dt <- nTot_dt[order(nTot_dt$nTot_obs, decreasing = TRUE),]

nTot_dt[nTot_dt$dataset == as.character(nTot_dt$dataset)[1],]
# dataset nTot_obs nPurityFlagged_obs nTot_sub
# 157 HMEC_40kb/TCGAbrca_lum_bas     1820                 84     1736
nTot_dt$nTot_obs <- NULL
nTot_dt_m <- melt(nTot_dt, id="dataset")
nTot_dt_m$dataset <- factor(nTot_dt_m$dataset, levels=as.character(nTot_dt$dataset))
stopifnot(!is.na(nTot_dt_m$dataset))

ggbarplot(nTot_dt_m, x="dataset", y="value", fill="variable") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), expand=c(0,10))+
  labs(y="# TADs", x="Datasets", fill="")+
  theme(
    text = element_text(family=fontFamily),
    panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.background = element_rect(fill = "transparent"),
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x =  element_blank(),
    axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
    # axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
    axis.text.x = element_blank(),
    plot.title = element_text(hjust=0.5, size = 16, face="bold"),
    plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
    legend.title = element_text(face="bold")
  )
# min(nTot_dt$nTot_sub)
# [1] 1493
resc_y <- 1400
nTot_dt <- in_dt[,c("dataset", "nTot_obs", "nPurityFlagged_obs")]
nTot_dt <- unique(nTot_dt)
nTot_dt$nTot_sub <- nTot_dt$nTot_obs - nTot_dt$nPurityFlagged_obs - resc_y
nTot_dt <- nTot_dt[order(nTot_dt$nTot_obs, decreasing = TRUE),]
nTot_dt$nTot_obs <- NULL
nTot_dt_m <- melt(nTot_dt, id="dataset")
nTot_dt_m$dataset <- factor(nTot_dt_m$dataset, levels=as.character(nTot_dt$dataset))
stopifnot(!is.na(nTot_dt_m$dataset))

y_v <- scales::pretty_breaks(n = 10)(nTot_dt_m$value)
y_l <- resc_y + scales::pretty_breaks(n = 10)(nTot_dt_m$value)
  
ggbarplot(nTot_dt_m, x="dataset", y="value", fill="variable") +
  labs(y="# TADs", x="Datasets", fill="")+
  # scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_y_continuous(breaks = y_v, labels = y_l, expand=c(0,10))+
  theme(
    text = element_text(family=fontFamily),
    panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.background = element_rect(fill = "transparent"),
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x =  element_blank(),
    axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
    # axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
    plot.title = element_text(hjust=0.5, size = 16, face="bold"),
    axis.text.x = element_blank(),
    plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
    legend.title = element_text(face="bold")
  )

  

nSignif_dt <- in_dt[,c("dataset", "nSignif_obs", "nSignifAndFlagged_obs")]
nSignif_dt <- unique(nSignif_dt)
nSignif_dt$nSignif_sub <- nSignif_dt$nSignif_obs - nSignif_dt$nSignifAndFlagged_obs
nSignif_dt <- nSignif_dt[order(nSignif_dt$nSignif_obs, decreasing = TRUE),]
nSignif_dt$nSignif_obs <- NULL
nSignif_dt_m <- melt(nSignif_dt, id="dataset")
nSignif_dt_m$dataset <- factor(nSignif_dt_m$dataset, levels=as.character(nTot_dt$dataset))
stopifnot(!is.na(nSignif_dt_m$dataset))

ggbarplot(nSignif_dt_m, x="dataset", y="value", fill="variable") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), expand=c(0,1))+
  labs(y="# signif. TADs", x="Datasets", fill="")+
  theme(
    text = element_text(family=fontFamily),
    panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.background = element_rect(fill = "transparent"),
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x =  element_blank(),
    axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
    # axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
    axis.text.x = element_blank(),
    plot.title = element_text(hjust=0.5, size = 16, face="bold"),
    plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
    legend.title = element_text(face="bold")
  )





obs_dt <- in_dt[,c("dataset", "ratioSignif_obs", "ratioFlagged_obs")]
obs_dt <- unique(obs_dt)



plot(x=obs_dt$ratioSignif_obs, y=obs_dt$ratioFlagged_obs, 
     xlab="ratioSignif_obs", ylab="ratioFlagged_obs",
     cex=0.7, pch=16, cex.lab=plotCex, cex.main=plotCex)
addCorr(x=obs_dt$ratioSignif_obs, y=obs_dt$ratioFlagged_obs)

























