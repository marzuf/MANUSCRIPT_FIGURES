options(scipen=100)

SSHFS=F

# Rscript sameFam_sameTADdiffTAD.R

script_name <- "sameFam_sameTADdiffTAD.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(ggplot2)
require(RColorBrewer)
require(reshape2)
# require(gplots)
registerDoMC(4)

source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")


buildTable <- TRUE

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 8)
myWidth <- myHeight
axisCex <- 1.4

myWidthGG <- 12
myHeightGG <- 12

source("../settings.R")
source("../full_dataset_names.R")

outFolder <- "SAMEFAM_SAMETADDIFFTAD"
dir.create(outFolder, recursive = TRUE)

family <- "hgnc"
familyType <- "hgnc_family_short"

hicds="Rao_HCT-116_2017_40kb"
exprds="TCGAcoad_msi_mss"

if(buildTable) {
  all_sameFam_auc_dt <- foreach(hicds=all_hicds, .combine='rbind') %dopar% {
    exprds_dt <- foreach(exprds=all_exprds[[paste0(hicds)]], .combine='rbind') %dopar% {
      auc_file <- file.path(runFolder, "AUC_COEXPRDIST_WITHFAM_SORTNODUP", hicds, paste0(exprds, "_",family), familyType, "auc_values.Rdata")
      if(!file.exists(auc_file)) return(NULL)
      all_auc <- get(load(auc_file))
      auc_sameFamSameTAD <- get(load(file.path(runFolder, "AUC_COEXPRDIST_WITHFAM_SORTNODUP", hicds, paste0(exprds, "_", family), familyType, "auc_sameFamSameTAD_distVect.Rdata")))
      auc_sameFamDiffTAD <- get(load(file.path(runFolder, "AUC_COEXPRDIST_WITHFAM_SORTNODUP", hicds, paste0(exprds, "_", family), familyType, "auc_sameFamDiffTAD_distVect.Rdata")))
      stopifnot(all_auc["auc_sameFamSameTAD_distVect"] ==auc_sameFamSameTAD)
      stopifnot(all_auc["auc_sameFamDiffTAD_distVect"] ==auc_sameFamDiffTAD)
      data.frame(hicds=hicds,exprds=exprds,auc_sameFamSameTAD=auc_sameFamSameTAD,auc_sameFamDiffTAD=auc_sameFamDiffTAD,stringsAsFactors = FALSE)
    }
    exprds_dt
  }
  outFile <- file.path(outFolder, "all_sameFam_auc_dt.Rdata")
  save(all_sameFam_auc_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "all_sameFam_auc_dt.Rdata")
  all_sameFam_auc_dt <- get(load(outFile))
}

all_sameFam_auc_dt$data_type <- gsub(".+_(.+)_40kb", "\\1", all_sameFam_auc_dt$hicds)
all_sameFam_auc_dt$data_type[all_sameFam_auc_dt$hicds %in% all_obs_hicds] <- "OBSERVED"
# unique(all_sameFam_auc_dt$data_type)
all_sameFam_auc_dt$sameTAD_diffTAD_ratio <- all_sameFam_auc_dt$auc_sameFamSameTAD/all_sameFam_auc_dt$auc_sameFamDiffTAD

all_sameFam_auc_dt$hicds_lab <- gsub("_RANDOM.+_40kb", "_40kb", all_sameFam_auc_dt$hicds)
all_sameFam_auc_dt$hicds_lab <- gsub("_PERMUT.+_40kb", "_40kb", all_sameFam_auc_dt$hicds_lab)
all_sameFam_auc_dt$dataset <- file.path(all_sameFam_auc_dt$hicds_lab, all_sameFam_auc_dt$exprds)

plot_dt <- reshape(all_sameFam_auc_dt[,c("dataset", "sameTAD_diffTAD_ratio", "data_type")], 
                   direction="wide", idvar=c("dataset"), timevar="data_type")

stopifnot(!duplicated(plot_dt$dataset))
plot_dt$dotcols <- all_cols[all_cmps[basename(as.character(plot_dt$dataset))]]
stopifnot(!is.na(plot_dt$dotcols))

plot_var <- "sameTAD_diffTAD_ratio"

rd_to_plot <- c("RANDOMMIDPOS","RANDOMMIDPOSDISC", "RANDOMMIDPOSSTRICT")
rd="RANDOMMIDPOSSTRICT"
for(rd in rd_to_plot) {
  
  xvar <- paste0(plot_var, ".OBSERVED")
  yvar <- paste0(plot_var, ".", rd)
  
  stopifnot(xvar %in% colnames(plot_dt))
  stopifnot(yvar %in% colnames(plot_dt))
  
  tmp_dt <- plot_dt[,c("dataset",xvar, yvar, "dotcols")]
  stopifnot(!is.na(tmp_dt))
  
  my_x <- tmp_dt[,c(xvar)]
  my_y <- tmp_dt[,c(yvar)]
  
  subT <- paste0(sum(my_x>1), " obs. >1; ",  sum(my_y > 1), " rd. > 1")
  
  outFile <- file.path(outFolder, paste0(plot_var, "_", rd, "_vs_", "obs.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  par(bty="L")
  plot(
    x = my_x,
    y = my_y,
    col = tmp_dt$dotcols,
    xlab = "observed",
    ylab = paste0(rd),
    pch=16,
    # cex=0.7,
    main = paste0("Same fam. coexpr ~ dist - same TAD/diffTAD AUC ratio"),
    cex.axis=plotCex,
    cex.main = plotCex,
    cex.lab = plotCex
  )
  abline(v=1, lty=2, col="grey")
  abline(h=1, lty=2, col="grey")
  mtext(side=3, text = paste0("# datasets = ", nrow(tmp_dt), "; ", subT), cex  =plotCex, font=3)
  addCorr(x=my_x,y=my_y, bty="n", legPos = "topleft")
  curve(1*x, col="darkgrey", add=T)
  legend("bottomright", legend=cmp_names[names(all_cols)], pch=16,col = all_cols, bty="n")
  # curve(1*x, col="darkgrey", add=T)
  #legend("bottomleft", 
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  
  
  
  
}






