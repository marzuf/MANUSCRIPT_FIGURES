options(scipen=100)

# v2 version
# -> one single script for all and for cmpType
# -> if share more than 80% of the genes -> merge conserved regions
# -> count as conserved in one dataset at the exprds level
# -> min conserved region

SSHFS=F

# Rscript selected_TADs_and_purity.R


args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 1) {
  purity_ds <- args[1]  
} else{
  purity_ds <- ""
}

script0_name <- "0_prepGeneData"


if(purity_ds == "") {
  purity_plot_name <- "aran"
  # all the ranks are between 1 and 0
} else if(purity_ds == "EPIC") {
  purity_plot_name <- "EPIC"
} else {
  stop("--invalid purity_ds\n")
}

#  purity_ds <- "CPE"
#  purity_plot_name <- "Aran - CPE"

### HARD-CODED - MAIN SETTINGS

corMet <- "pearson"
transfExpr <- "log10"
signifThresh <- 0.01
corrPurityQtThresh <- 0.05
signif_column <- "adjPvalComb"
signifcol <- paste0(signif_column, "_", signifThresh)

# args <- commandArgs(trailingOnly = TRUE)
# if(length(args) == 1) {
#   purity_ds <- args[1]  
#   purity_plot_name <- "EPIC"
# } else{
#   purity_ds <- ""
#   purity_plot_name <- "aran"
# }

script_name <- "selected_TADs_and_purity.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(GenomicRanges)
require(ggplot2)
require(RColorBrewer)
require(reshape2)
# require(gplots)
registerDoMC(4)

source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("plot_lolliTAD_funct.R")
# source("my_heatmap.2.R")

buildTable <- TRUE

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

myHeightHeat <- myHeight * 1.8
myWidthHeat <- myWidth * 1.8

myWidthGG <- 12
myHeightGG <- 12




selected_tads <- c("ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/chr11_TAD390",
                   "ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/chr10_TAD268",
                   "ENCSR489OCU_NCI-H460_40kb/TCGAluad_mutKRAS_mutEGFR/chr10_TAD16",
                   "ENCSR489OCU_NCI-H460_40kb/TCGAluad_mutKRAS_mutEGFR/chr17_TAD162")
                   

outFolder <- file.path("SELECTED_TADS_AND_PURITY", purity_plot_name, transfExpr)
dir.create(outFolder, recursive = TRUE)

purity_file <- file.path("../../v2_Yuanlong_Cancer_HiC_data_TAD_DA/ALLTADS_AND_PURITY", purity_ds, transfExpr, "all_ds_corrPurity_dt.Rdata")
purityData <- get(load(purity_file))
agg_purity <- aggregate(purityCorr~dataset+region, FUN=mean, data=purityData)

agg_purity$regID <- file.path(agg_purity$dataset, agg_purity$region)

result_file <- file.path("../../v2_Yuanlong_Cancer_HiC_data_TAD_DA/CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)

merge_dt <- merge(agg_purity, resultData, by=c("dataset", "region"))
merge_dt$signif <- merge_dt$adjPvalComb <= signifThresh
merge_dt$signif_lab <- ifelse(merge_dt$signif, paste0("adj. p-val <=", signifThresh), paste0("adj. p-val >", signifThresh) )
purityCorrThresh <- as.numeric(quantile(merge_dt$purityCorr[!merge_dt$signif], probs = corrPurityQtThresh ))

subTit <- paste0(corMet, "'s corr.", " - ", purity_plot_name, " data")
plotTit <- paste0("Purity corr. distribution")
myx_lab <- paste0(transfExpr, " expr. and purity correlation (meanTAD)")

agg_purity$regID_lab <- gsub("/chr", "\nchr", agg_purity$regID)

outFile <- file.path(outFolder, paste0("exprPurityCorr_meanTAD_signif_notSignif_selectedTADs_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.5))
# par(mar=c(11,4,4,4))
plot_multiDens(split(merge_dt$purityCorr, merge_dt$signif_lab),
               plotTit = plotTit, my_xlab = myx_lab)
lines(density(merge_dt$purityCorr), col="green")
abline(v=purityCorrThresh, col="darkgrey")
legend("topleft", lty=1, lwd=2, col=c("green", "darkgrey"), bty="n", legend=c("all", paste0(corrPurityQtThresh, "-qt non-signif. TADs\n(=", round(purityCorrThresh, 2), ")")))
abline(v=agg_purity$purityCorr[agg_purity$regID %in% selected_tads], col="blue")
# mtext(side=1, at = agg_purity$purityCorr[agg_purity$regID %in% selected_tads], text=agg_purity$regID_lab[agg_purity$regID %in% selected_tads],  col="blue", las=2, cex=0.6, padj = c(0,1))
text_dt <- agg_purity[agg_purity$regID %in% selected_tads,]
text_dt <- text_dt[order(text_dt$purityCorr),]
y_pos <- seq(from=3, length.out=nrow(text_dt), by=-0.4)
text(x=text_dt$purityCorr,
     y=y_pos, 
      labels=text_dt$regID_lab  ,
     adj=c(1,1),
     col="blue", cex=0.6)
mtext(side=3, text = subTit)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



