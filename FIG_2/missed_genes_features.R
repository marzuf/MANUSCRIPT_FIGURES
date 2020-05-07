
# Rscript missed_genes_features.R

plotType <- "svg"
source("../settings.R")


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

plotType <- "svg"
myHeight <- ifelse(plotType == "png", 480 , 7)
myWidth <- ifelse(plotType == "png", 600, 10)
plotCex <- 1.4

require(flux)
require(foreach)
require(doMC)
registerDoMC(nCpu)

geneSignifThresh <- 0.01
tadSignifThresh <- 0.01

all_ranks_dt <- get(load(file.path(runFolder, "GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata")))


outFolder <- file.path("MISSED_GENES_FEATURES")
dir.create(outFolder, recursive=TRUE)

options(scipen=100)

startTime <- Sys.time()

buildData <- TRUE

hicds <- "ENCSR489OCU_NCI-H460_40kb"
exprds <- "TCGAlusc_norm_lusc"
# all_obs_hicds=all_obs_hicds[1]
if(buildData){
  
  all_dt <- foreach(hicds = all_obs_hicds, .combine='rbind') %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    hicds_dt <- foreach(exprds = all_obs_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      ds_rank_dt <- all_ranks_dt[all_ranks_dt$hicds == hicds & all_ranks_dt$exprds == exprds, ]
      
      tad_rD <- get(load(file.path(pipFolder, hicds, exprds, "8cOnlyRatioDownFastSave_runAllDown", "all_obs_ratioDown.Rdata")))
      tad_rD_dt <- data.frame(
        region = names(tad_rD),
        TAD_rD = as.numeric(tad_rD), 
        stringsAsFactors = FALSE
      )
      tad_meanFC <- get(load(file.path(pipFolder, hicds, exprds, "3_runMeanTADLogFC", "all_meanLogFC_TAD.Rdata")))
      tad_meanFC_dt <- data.frame(
        region = names(tad_meanFC),
        TAD_meanFC = as.numeric(tad_meanFC), 
        stringsAsFactors = FALSE
      )
      
      tad_maxFC <- get(load(file.path(pipFolder, hicds, exprds, "3max_runMaxTADLogFC", "all_maxLogFC_TAD.Rdata")))
      tad_maxFC_dt <- data.frame(
        region = names(tad_maxFC),
        TAD_maxFC = as.numeric(tad_maxFC), 
        stringsAsFactors = FALSE
      )
      tad_varFC <- get(load(file.path(pipFolder, hicds, exprds, "3var_runVarTADLogFC", "all_varLogFC_TAD.Rdata")))
      tad_varFC_dt <- data.frame(
        region = names(tad_varFC),
        TAD_varFC = as.numeric(tad_varFC), 
        stringsAsFactors = FALSE
      )
      
      # !!! NOT ADJUSTED !!!
      tad_meanCorrPval <- get(load(file.path(pipFolder, hicds, exprds, "10sameNbr_runEmpPvalMeanTADCorr", "emp_pval_meanCorr.Rdata")))
      tad_meanCorrPval_dt <- data.frame(
        region = names(tad_meanCorrPval),
        TAD_meanCorrPval = as.numeric(tad_meanCorrPval), 
        stringsAsFactors = FALSE
      )
      tad_meanCorrPval_dt$TAD_meanCorrPval_log10 <- -log10(tad_meanCorrPval_dt$TAD_meanCorrPval)
      
      tad_logFCpval <- get(load(file.path(pipFolder, hicds, exprds, "9_runEmpPvalMeanTADLogFC", "emp_pval_meanLogFC.Rdata")))
      tad_logFCpval_dt <- data.frame(
        region = names(tad_logFCpval),
        TAD_logFCpval = as.numeric(tad_logFCpval), 
        stringsAsFactors = FALSE
      )
      tad_logFCpval_dt$TAD_logFCpval_log10 <- -log10(tad_logFCpval_dt$TAD_logFCpval)
      
      
      # missed_genes_dt <- ds_rank_dt[ds_rank_dt$adj.P.Val > geneSignifThresh &
      #                              ds_rank_dt$tad_adjCombPval <= tadSignifThresh,]
      
      
      # dt <- merge(missed_genes_dt, tad_maxFC_dt, by = c("region"), all.x=TRUE, all.y=FALSE)
      # stopifnot(!is.na(dt))
      # 
      # dt <- merge(dt, tad_varFC_dt, by = c("region"), all.x=T, all.y=F)
      # stopifnot(!is.na(dt))
      
      dt <- merge(ds_rank_dt, tad_maxFC_dt, by = c("region"), all=T)
      stopifnot(!is.na(dt))
      dt <- merge(dt, tad_meanFC_dt, by = c("region"), all=T)
      stopifnot(!is.na(dt))
      dt <- merge(dt, tad_varFC_dt, by = c("region"), all=T)
      stopifnot(!is.na(dt))
      dt <- merge(dt, tad_rD_dt, by = c("region"), all=T)
      stopifnot(!is.na(dt))
      dt <- merge(dt, tad_logFCpval_dt, by = c("region"), all=T)
      stopifnot(!is.na(dt))
      dt <- merge(dt, tad_meanCorrPval_dt, by = c("region"), all=T)
      stopifnot(!is.na(dt))
      
      dt
    }
  }  
  outFile <- file.path(outFolder, paste0("all_dt.Rdata"))
  auc_fract_file <- outFile
  save(all_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  inFile <- file.path(outFolder, paste0("all_dt.Rdata"))
  all_dt <- get(load(inFile))
  # load("COORD_DETECT/all_dt.Rdata")
}

# stop("-ok")

rd_ranges <- seq(0,1,by=0.1)
pval_ranges <- c(seq(0,0.05,by=0.01), 1)
lowFC_thresh <- 0.5
rD_lowThresh <- 1/3
rD_highThresh <- 2/3

all_dt$geneSignif <- all_dt$adj.P.Val <= geneSignifThresh
all_dt$tadSignif <- all_dt$tad_adjCombPval <= tadSignifThresh
all_dt$sameDirFC <- sign(all_dt$logFC) == sign(all_dt$TAD_meanFC)
all_dt$lowFC <- abs(all_dt$logFC) <= lowFC_thresh
all_dt$lowSameDirFC <- all_dt$sameDirFC & all_dt$lowFC
all_dt$lowDiffDirFC <- !all_dt$sameDirFC & all_dt$lowFC
all_dt$missedGeneSignif <- !all_dt$geneSignif & all_dt$tadSignif
all_dt$missedLowFCsameDir <- all_dt$lowSameDirFC & all_dt$missedGeneSignif
all_dt$missedLowFCdiffDir <- all_dt$lowDiffDirFC & all_dt$missedGeneSignif
# plot(all_dt$tad_adjCombPval~all_dt$TAD_rD)

all_dt$dataset <- file.path(hicds, exprds)
nDS <- length(unique(all_dt$dataset))

all_dt$rD_range <- cut(all_dt$TAD_rD, breaks=rd_ranges, include.lowest=TRUE)
stopifnot(!is.na(all_dt$rD_range))

all_dt$pvalCorr_range <- cut(all_dt$TAD_meanCorrPval, breaks=pval_ranges, include.lowest=TRUE)
stopifnot(!is.na(all_dt$pvalCorr_range))

save(all_dt, file=file.path(outFolder, "all_dt2.Rdata"), version=2)

################################
# Boxplot 1: ratioMissedGenes by rD_range (for each dataset ratio of missed genes out of all genes)
################################

ratioMissed_dt <- aggregate(missedGeneSignif ~ hicds + exprds + rD_range, FUN=mean, data=all_dt)

plot_dt <- "ratioMissed_dt"
plot_var <- "missedGeneSignif"

plotDesc <- "Ratio missed genes (over all genes, by dataset)"
subTit <- paste0("gene-level adj. p-val > ",  geneSignifThresh, "TAD-level adj. p-val <= ",  tadSignifThresh)

outFile <- file.path(outFolder, paste0(plot_var, "_ratio.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(as.formula(paste0(plot_var, "~ ", "rD_range")), ylab=paste0(plot_var, " (ratio)"), cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, main = plotDesc,
        xlab="ratioDown range",data = get(plot_dt))
mtext(side=3, text = subTit)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

################################
# Boxplot 2: actual number of missed genes per dataset
################################

nMissed_dt <- aggregate(missedGeneSignif ~ hicds + exprds + rD_range, FUN=sum, data=all_dt)

plot_dt <- "nMissed_dt"
plot_var <- "missedGeneSignif"

plotDesc <- "# missed genes (by dataset)"
subTit <- paste0("gene-level adj. p-val > ",  geneSignifThresh, "; TAD-level adj. p-val <= ",  tadSignifThresh)

outFile <- file.path(outFolder, paste0(plot_var, "_sum.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(as.formula(paste0(plot_var, "~ ", "rD_range")), ylab=paste0(plot_var, " (#)"), cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, main = plotDesc, xlab="ratioDown range",data = get(plot_dt))
mtext(side=3, text = subTit)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

################################
# Boxplot 3: actual number of missed genes with lowFCsameDir per dataset
################################
nMissedLowFCsameDir_dt <- aggregate(missedLowFCsameDir ~ hicds + exprds + rD_range, FUN=sum, data=all_dt)

plot_dt <- "nMissedLowFCsameDir_dt"
plot_var <- "missedLowFCsameDir"

plotDesc <- "# missedLowFCsameDir genes (by dataset)"
subTit <- paste0("gene-level adj. p-val > ",  geneSignifThresh, "; TAD-level adj. p-val <= ",  tadSignifThresh,
                 "; abs(logFC) <= ", lowFC_thresh)

outFile <- file.path(outFolder, paste0(plot_var, "_sum.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(as.formula(paste0(plot_var, "~ ", "rD_range")), ylab=paste0(plot_var, " (#)"), cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, main = plotDesc, xlab="ratioDown range",data = get(plot_dt))
mtext(side=3, text = subTit)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

################################
# Boxplot 4: out of all genes, ratio of lowFCsameDir
################################
ratioMissedLowFCsameDir_dt <- aggregate(missedLowFCsameDir ~ hicds + exprds + rD_range, FUN=mean, data=all_dt)


plotDesc <- "ratio missedLowFCsameDir genes (by dataset)"
subTit <- paste0("gene-level adj. p-val > ",  geneSignifThresh, "; TAD-level adj. p-val <= ",  tadSignifThresh,
                 "; abs(logFC) <= ", lowFC_thresh)


plot_dt <- "ratioMissedLowFCsameDir_dt"
plot_var <- "missedLowFCsameDir"

outFile <- file.path(outFolder, paste0(plot_var, "_ratio.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(as.formula(paste0(plot_var, "~ ", "rD_range")), ylab=paste0(plot_var, " (ratio)"), cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, main = plotDesc, xlab="ratioDown range",data = get(plot_dt))
mtext(side=3, text = subTit)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


################################
# Boxplot 4: out of the missed genes, ratio of lowFCsameDir
################################


missedRatios_dt <- merge(nMissed_dt, nMissedLowFCsameDir_dt, by=c("hicds", "exprds", "rD_range"), all=TRUE)
stopifnot(!is.na(missedRatios_dt))
missedRatios_dt$ratioLowFCsameDir <- missedRatios_dt$missedLowFCsameDir/missedRatios_dt$missedGeneSignif
stopifnot(na.omit(missedRatios_dt)$ratioLowFCsameDir >= 0 & na.omit(missedRatios_dt)$ratioLowFCsameDir <= 1)


plotDesc <- "ratio missedLowFCsameDir/missed (by dataset)"
subTit <- paste0("gene-level adj. p-val > ",  geneSignifThresh, "; TAD-level adj. p-val <= ",  tadSignifThresh,
                 "; abs(logFC) <= ", lowFC_thresh)


plot_dt <- "missedRatios_dt"
plot_var <- "ratioLowFCsameDir"

outFile <- file.path(outFolder, paste0(plot_var, "_missedLowFCsameDir_missed_ratio.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(as.formula(paste0(plot_var, "~ ", "rD_range")), ylab=plot_var, cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, main = plotDesc, xlab="ratioDown range",data = get(plot_dt))
mtext(side=3, text = subTit)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

################################
# Boxplot 4: out of the missed genes, ratio of lowFCdiffDir
################################
nMissedLowFCdiffDir_dt <- aggregate(missedLowFCdiffDir ~ hicds + exprds + rD_range, FUN=sum, data=all_dt)
missedRatiosDiffDir_dt <- merge(nMissed_dt, nMissedLowFCdiffDir_dt, by=c("hicds", "exprds", "rD_range"), all=TRUE)
stopifnot(!is.na(missedRatiosDiffDir_dt))
missedRatiosDiffDir_dt$ratioLowFCdiffDir <- missedRatiosDiffDir_dt$missedLowFCdiffDir/missedRatiosDiffDir_dt$missedGeneSignif
stopifnot(na.omit(missedRatiosDiffDir_dt)$ratioLowFCdiffDir >= 0 & na.omit(missedRatiosDiffDir_dt)$ratioLowFCdiffDir <= 1)

plot_dt <- "missedRatiosDiffDir_dt"
plot_var <- "ratioLowFCdiffDir"

plotDesc <- "ratio missedLowFCdiffDir/missedGeneSignif (by dataset)"
subTit <- paste0("gene-level adj. p-val > ",  geneSignifThresh, "; TAD-level adj. p-val <= ",  tadSignifThresh,
                 "; abs(logFC) <= ", lowFC_thresh)

outFile <- file.path(outFolder, paste0(plot_var, "_missedLowFCdiffDir_missed_ratio.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(as.formula(paste0(plot_var, "~ ", "rD_range")), ylab=plot_var, cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, main = plotDesc, xlab="ratioDown range",data = get(plot_dt))
mtext(side=3, text = subTit)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



################################
# Boxplot 5: total detected genes tad-level
################################
tadSignif_dt <- aggregate(tadSignif~hicds + exprds + rD_range, FUN=sum, data=all_dt)
plot_dt <- "tadSignif_dt"
plot_var <- "tadSignif"

plotDesc <- "# TAD-level signif. genes"
subTit <- paste0("TAD-level adj. p-val <= ",  tadSignifThresh)

outFile <- file.path(outFolder, paste0(plot_var, "_tadSignif_sum.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(as.formula(paste0(plot_var, "~ ", "rD_range")), ylab=paste0(plot_var, " (#)"), cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, main = plotDesc, xlab="ratioDown range",data = get(plot_dt))
mtext(side=3, text = subTit)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



################################
# Boxplot 5: total detected genes gene-level
################################
geneSignif_dt <- aggregate(geneSignif~hicds + exprds + rD_range, FUN=sum, data=all_dt)
plot_dt <- "geneSignif_dt"
plot_var <- "geneSignif"

plotDesc <- "# gene-level signif. genes"
subTit <- paste0("gene-level adj. p-val <= ",  geneSignifThresh)

outFile <- file.path(outFolder, paste0(plot_var, "_geneSignif_sum.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(as.formula(paste0(plot_var, "~ ", "rD_range")), ylab=paste0(plot_var, " (#)"), cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, main = plotDesc, xlab="ratioDown range",data = get(plot_dt))
mtext(side=3, text = subTit)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



################################
# Boxplot 6: out of all TAD signif genes, ratio missed genes
################################

signif_missed_dt <- merge(tadSignif_dt, nMissed_dt, by=c("hicds", "exprds", "rD_range"), all=TRUE)
signif_missed_dt$ratioMissed <- signif_missed_dt$missedGeneSignif/signif_missed_dt$tadSignif
stopifnot(na.omit(signif_missed_dt)$ratioMissed >= 0 & na.omit(signif_missed_dt)$ratioMissed  <= 1)

plot_dt <- "signif_missed_dt"
plot_var <- "ratioMissed"

plotDesc <- "ratio missedGeneSignif/tadSignif (by dataset)"
subTit <- paste0("gene-level adj. p-val > ",  geneSignifThresh, "; TAD-level adj. p-val <= ",  tadSignifThresh)

outFile <- file.path(outFolder, paste0(plot_var, "_missedGeneSignif_tadSignif_ratio.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(as.formula(paste0(plot_var, "~ ", "rD_range")), ylab=plot_var, cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, main = plotDesc, xlab="ratioDown range",data = get(plot_dt))
mtext(side=3, text = subTit)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



################################
# Boxplot 7: out of all TAD signif genes, ratio missedLowFCsameDir
################################

signif_missedLowFCsameDir_dt <- merge(tadSignif_dt, nMissedLowFCsameDir_dt, by=c("hicds", "exprds", "rD_range"), all=TRUE)
signif_missedLowFCsameDir_dt$ratioMissedLowFCsameDir <- signif_missedLowFCsameDir_dt$missedLowFCsameDir/signif_missedLowFCsameDir_dt$tadSignif
stopifnot(na.omit(signif_missedLowFCsameDir_dt)$ratioMissedLowFCsameDir >= 0 & na.omit(signif_missedLowFCsameDir_dt)$ratioMissedLowFCsameDir  <= 1)

plot_dt <- "signif_missedLowFCsameDir_dt"
plot_var <- "ratioMissedLowFCsameDir"

plotDesc <- "ratio missedLowFCsameDir/tadSignif (by dataset)"
subTit <- paste0("gene-level adj. p-val > ",  geneSignifThresh, "; TAD-level adj. p-val <= ",  tadSignifThresh,
                 "; abs(logFC) <= ", lowFC_thresh)
outFile <- file.path(outFolder, paste0(plot_var, "_missedLowFCsameDir_tadSignif_ratio.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(as.formula(paste0(plot_var, "~ ", "rD_range")), ylab=plot_var, cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, main = plotDesc, xlab="ratioDown range",data = get(plot_dt))
mtext(side=3, text = subTit)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


################################
# Boxplot 1: # missed lowFC genes by TAD and rD_range 
################################
nMissedByTAD_dt <- aggregate(missedLowFCsameDir ~ hicds + exprds + rD_range + region, FUN=sum, data=all_dt)


plotDesc <- "# missedLowFCsameDir (by TAD)"
subTit <- paste0("gene-level adj. p-val > ",  geneSignifThresh, "; TAD-level adj. p-val <= ",  tadSignifThresh,
                 "; abs(logFC) <= ", lowFC_thresh)


plot_dt <- "nMissedByTAD_dt"
plot_var <- "missedLowFCsameDir"

outFile <- file.path(outFolder, paste0(plot_var, "_missedLowFCsameDir_byTAD_sum.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
boxplot(as.formula(paste0(plot_var, "~ ", "rD_range")), ylab=paste0(plot_var, " (#)"), cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, main = plotDesc, xlab="ratioDown range",data = get(plot_dt))
mtext(side=3, text = subTit)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



# min1_nMissedByTAD_dt <- nMissedByTAD_dt[nMissedByTAD_dt$missedLowFCdiffDir > 0,]
# min2_nMissedByTAD_dt <- nMissedByTAD_dt[nMissedByTAD_dt$missedLowFCdiffDir > 1,]
# min2_nMissedByTAD_dt$ds_id <- file.path(min2_nMissedByTAD_dt$hicds, min2_nMissedByTAD_dt$exprds, min2_nMissedByTAD_dt$region)
# 
# plot_dt <- "min1_nMissedByTAD_dt"
# plot_var <- "missedLowFCdiffDir"
# 
# plotDesc <- "ratio missedLowFCsameDir/tadSignif"
# subTit <- paste0("gene-level adj. p-val > ",  geneSignifThresh, "; TAD-level adj. p-val <= ",  tadSignifThresh,
#                  "; abs(logFC) <= ", lowFC_thresh)
# 
# boxplot(as.formula(paste0(plot_var, "~ ", "rD_range")), ylab=plot_var, cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, main = plotDesc, xlab="ratioDown range",data = get(plot_dt))
# 
# plot_dt <- "min2_nMissedByTAD_dt"
# boxplot(as.formula(paste0(plot_var, "~ ", "rD_range")), ylab=plot_var, cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, main = plotDesc, xlab="ratioDown range",data = get(plot_dt))
# 
# multipleMissed_TADs <- min2_nMissedByTAD_dt$ds_id

###################################################################
# rD distribution of missed genes 
###################################################################

source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")

missed_dt <- all_dt[all_dt$missedGeneSignif,]



plotDesc <- "TAD logFC distribution of missed genes (missedLowFCsameDir vs. other missed)"
subTit <- paste0("gene-level adj. p-val > ",  geneSignifThresh, "; TAD-level adj. p-val <= ",  tadSignifThresh,
                 "; abs(logFC) <= ", lowFC_thresh)

plot_var <- "meanLogFC"
outFile <- file.path(outFolder, paste0(plot_var, "_missedGenes_lowFCsameDir_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(
  lapply(split(missed_dt, missed_dt$missedLowFCsameDir), function(x) x[["logFC"]]), plotTit = plotDesc
)
mtext(side=3, text = subTit)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



plotDesc <- "TAD ratioDown distribution of missed genes (missedLowFCsameDir vs. other missed)"
subTit <- paste0("gene-level adj. p-val > ",  geneSignifThresh, "; TAD-level adj. p-val <= ",  tadSignifThresh,
                 "; abs(logFC) <= ", lowFC_thresh)


plot_var <- "ratioDown"
outFile <- file.path(outFolder, paste0(plot_var, "_missedGenes_lowFCsameDir_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

plot_multiDens(
  lapply(split(missed_dt, missed_dt$missedLowFCsameDir), function(x) x[["TAD_rD"]]), plotTit = plotDesc
)
mtext(side=3, text = subTit)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


plotDesc <- "TAD varFC distribution of missed genes (missedLowFCsameDir vs. other missed)"
subTit <- paste0("gene-level adj. p-val > ",  geneSignifThresh, "; TAD-level adj. p-val <= ",  tadSignifThresh,
                 "; abs(logFC) <= ", lowFC_thresh)

plot_var <- "meanVarFC"
outFile <- file.path(outFolder, paste0(plot_var, "_missedGenes_lowFCsameDir_density.", plotType))

do.call(plotType, list(outFile, height=myHeight, width=myWidth))

plot_multiDens(
  lapply(split(missed_dt, missed_dt$missedLowFCsameDir), function(x) log10(x[["TAD_varFC"]])), plotTit = plotDesc
)
mtext(side=3, text = subTit)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

# ###################################################################
# # FC distribution of missed genes 
# ###################################################################
# 
# plotDesc <- "LogFC distribution of missed genes"
# subTit <- paste0("gene-level adj. p-val > ",  geneSignifThresh, "; TAD-level adj. p-val <= ",  tadSignifThresh,
#                  "; abs(logFC) <= ", lowFC_thresh)
# 
# plot(density(missed_dt$logFC))


###################################################################
# gene focus 
###################################################################

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

all_dt_2=all_dt
save(all_dt_2, file=file.path(outFolder, "all_dt_2.Rdata"), version=2)
load("MISSED_GENES_FEATURES/all_dt2.Rdata")

all_dt$focusGenes <- all_dt$missedLowFCsameDir & 
                     (all_dt$TAD_rD >= rD_highThresh | all_dt$TAD_rD <= rD_lowThresh)

focus_dt <- all_dt[all_dt$focusGenes,]

focus_dt$dataset <- file.path(focus_dt$hicds, focus_dt$exprds )

focus_dt$ds_id <- file.path(focus_dt$hicds,focus_dt$exprds, focus_dt$region)
# the # of times (# of datasets) the genes has been missed
nRec_dt <- aggregate(ds_id  ~ entrezID, data=focus_dt, FUN=function(x)length(unique(x)))
# the datasets in which the genes has been missed
dsRec_dt <- aggregate(ds_id  ~ entrezID, data=focus_dt, FUN=function(x)paste0(unique(x), collapse=","))

recFocus_dt <- merge(nRec_dt, dsRec_dt, by="entrezID", all=T)
colnames(recFocus_dt)[2:3] <- c("missed_nDS", "missed_DS")
stopifnot(!is.na(recFocus_dt))                     
recFocus_dt <- recFocus_dt[order(recFocus_dt$missed_nDS, decreasing = TRUE),] 

recFocus_dt$geneSymbol <- entrez2symb[paste0(recFocus_dt$entrezID)]
stopifnot(!is.na(recFocus_dt$geneSymbol))

# recFocus_dt => missedLowFCsameDir + belonging to rD >= and <= rD
save(recFocus_dt, file="recFocus_dt.Rdata", version=2)
nrow(recFocus_dt)
# 713
recFocus_dt <- recFocus_dt[,c("entrezID", "geneSymbol", "missed_nDS" ,"missed_DS" )]
outFile <- file.path(outFolder, "recurrently_missedLowFCsameDir_genes.txt")
write.table(recFocus_dt, file = outFile, sep="\t", col.names = T, row.names = FALSE, quote=F)
cat(paste0("... written: ", outFile, "\n"))


###################################################################
# gene focus gene multiple TADs
###################################################################

focus_symb_dt <-focus_dt
focus_symb_dt$geneSymbol <- entrez2symb[paste0(focus_symb_dt$entrezID)]
stopifnot(!is.na(focus_symb_dt$geneSymbol))
missedLowFCsameDir_genes_dt <- aggregate(geneSymbol~ds_id, data=focus_symb_dt, FUN=function(x) paste0(x, collapse=","))
colnames(missedLowFCsameDir_genes_dt)[colnames(missedLowFCsameDir_genes_dt)=="geneSymbol"] <- "missedLowFCsameDirGenes"

all_symb_dt <- all_dt
all_symb_dt$geneSymbol <- entrez2symb[paste0(all_symb_dt$entrezID)]
all_symb_dt$ds_id <- file.path(all_symb_dt$hicds, all_symb_dt$exprds, all_symb_dt$region)
stopifnot(!is.na(all_symb_dt$geneSymbol))
all_genes_dt <- aggregate(geneSymbol~ds_id, data=all_symb_dt, FUN=function(x) paste0(x, collapse=","))
# save(missedLowFCsameDir_genes_dt, file="missedLowFCsameDir_genes_dt.Rdata", version=2)
colnames(all_genes_dt)[colnames(all_genes_dt)=="geneSymbol"] <- "allGenes"

gene_dt <- merge(missedLowFCsameDir_genes_dt, all_genes_dt, by="ds_id", all=T)

nMissedByTAD_dt <- aggregate(missedLowFCsameDir ~ hicds + exprds + rD_range + region + ds_id, FUN=sum, data=focus_dt)
table(nMissedByTAD_dt$missedLowFCsameDir)
# 1   2   3   4   5   6   7   8   9  10  11  12  14 
# 286 134  95  54  27  16  10   4   3   3   2   1   1 

# => look at TAD with multiple missed genes but not multiple missed focus genes

min1_nMissedByTAD_focus_dt <- nMissedByTAD_dt[nMissedByTAD_dt$missedLowFCsameDir > 0,]
nrow(min1_nMissedByTAD_focus_dt)
min2_nMissedByTAD_focus_dt <- nMissedByTAD_dt[nMissedByTAD_dt$missedLowFCsameDir > 1,]
nrow(min2_nMissedByTAD_focus_dt)
min2_nMissedByTAD_focus_dt$ds_id <- file.path(min2_nMissedByTAD_focus_dt$hicds, min2_nMissedByTAD_focus_dt$exprds, min2_nMissedByTAD_focus_dt$region)

table(min2_nMissedByTAD_focus_dt$missedLowFCsameDir)

min2_nMissedByTAD_focus_dt[order(min2_nMissedByTAD_focus_dt$missedLowFCsameDir, decreasing = TRUE),]

focus_dt$ds_id <- file.path(focus_dt$hicds, focus_dt$exprds, focus_dt$region)

multiple_focus_dt <- focus_dt[focus_dt$ds_id %in% min2_nMissedByTAD_focus_dt$ds_id,]
multiple_focus_dt$geneSymbol <- entrez2symb[paste0(multiple_focus_dt$entrezID)]


out_dt <- nMissedByTAD_dt[nMissedByTAD_dt$missedLowFCsameDir > 2,]

save(out_dt, file="out_dt.Rdata", version=2)
save(gene_dt, file="gene_dt.Rdata", version=2)
     

out_dt <- merge(out_dt, gene_dt, all.x=T, all.y=F, by="ds_id")
out_dt <- out_dt[order(out_dt$missedLowFCsameDir, decreasing = TRUE),]
nrow(out_dt)
outFile <- file.path(outFolder, "tad_with_multiple_missedLowFCsameDir_genes.txt")
write.table(out_dt, file = outFile, sep="\t", col.names = T, row.names = FALSE, quote=F)
cat(paste0("... written: ", outFile, "\n"))


### min. 3 missed genes by TAD

# min3_nMissedByTAD_focus_dt <- nMissedByTAD_dt[nMissedByTAD_dt$missedLowFCdiffDir > 2,]
# min3_nMissedByTAD_focus_dt$ds_id <- file.path(min3_nMissedByTAD_focus_dt$hicds, min3_nMissedByTAD_focus_dt$exprds, min3_nMissedByTAD_focus_dt$region)
# min3_nMissedByTAD_focus_dt$ds_id
# 
# multiple_focus_dt <- focus_dt[focus_dt$ds_id %in% min3_nMissedByTAD_focus_dt$ds_id,]
# multiple_focus_dt$geneSymbol <- entrez2symb[paste0(multiple_focus_dt$entrezID)]
# stopifnot(!is.na(multiple_focus_dt$geneSymbol))

# # min3_nMissedByTAD_focus_dt$ds_id
# [1] "ENCSR444WCZ_A549_40kb/TCGAluad_nonsmoker_smoker/chr16_TAD3"     
# [2] "LG2_40kb/TCGAluad_nonsmoker_smoker/chr16_TAD3"                  
# [3] "PA2_40kb/TCGApaad_wt_mutKRAS/chr19_TAD202"                      
# [4] "ENCSR862OGI_RPMI-7951_40kb/TCGAskcm_lowInf_highInf/chr19_TAD215"
# [5] "K562_40kb/TCGAlaml_wt_mutFLT3/chr19_TAD3"                       
# [6] "ENCSR862OGI_RPMI-7951_40kb/TCGAskcm_wt_mutCTNNB1/chr6_TAD119"   
# Rscript look_TAD_expression_withRank_v2.R ENCSR444WCZ_A549_40kb TCGAluad_nonsmoker_smoker chr16_TAD3
# Rscript look_TAD_expression_withRank_v2.R LG2_40kb TCGAluad_nonsmoker_smoker chr16_TAD3
# Rscript look_TAD_expression_withRank_v2.R PA2_40kb TCGApaad_wt_mutKRAS chr19_TAD202
# Rscript look_TAD_expression_withRank_v2.R ENCSR862OGI_RPMI-7951_40kb TCGAskcm_lowInf_highInf chr19_TAD215
# Rscript look_TAD_expression_withRank_v2.R K562_40kb TCGAlaml_wt_mutFLT3 chr19_TAD3
# Rscript look_TAD_expression_withRank_v2.R ENCSR862OGI_RPMI-7951_40kb TCGAskcm_wt_mutCTNNB1 chr6_TAD119

stop("-ok\n")


# * TRASH *#




lowFC_thresh <- 0.5
rD_lowThresh <- 1/3
rD_highThresh <- 2/3
# genes of interest:
# - low FC
# - high RD
# - same sign as max
# - in signif. TADs


sum(all_dt$lowFC)
sum(all_dt$lowFC & all_dt$tadSignif)
sum(all_dt$lowSameDirFC & all_dt$tadSignif)

sum( abs(all_dt$logFC) <= lowFC_thresh &
  (all_dt$TAD_rD >= rD_highThresh | all_dt$TAD_rD <= rD_lowThresh) &
  sign(all_dt$TAD_maxFC) == sign(all_dt$logFC) & 
    all_dt$tadSignif & 
    !all_dt$geneSignif)
# 738

sum(  sign(all_dt$TAD_maxFC) == sign(all_dt$logFC) & 
        (all_dt$TAD_rD >= rD_highThresh | all_dt$TAD_rD <= rD_lowThresh) &
        abs(all_dt$logFC) <= lowFC_thresh )
# 131277



signif_dt <- all_dt[all_dt$tadSignif,]
sum(signif_dt$lowSameDirFC)

# how many low FC by TAD that are in the same dir as TADmeanFC

agg_signif_dt <- aggregate(lowSameDirFC ~ region+hicds+exprds+rD_range, data = signif_dt, FUN=mean )

agg_signif_dt[agg_signif_dt$region=="chr19_TAD5"& 
                agg_signif_dt$hicds=="Barutcu_MCF-10A_40kb"&
                agg_signif_dt$exprds=="TCGAbrca_lum_bas",]

boxplot(lowSameDirFC~rD_range,data=agg_signif_dt)

agg2_signif_dt <- aggregate(missGeneSignif ~ hicds+exprds+rD_range, data = all_dt, FUN=sum)
boxplot(missGeneSignif~rD_range,data=agg2_signif_dt)

agg2b_signif_dt <- aggregate(lowSameDirFC ~ hicds+exprds+rD_range, data = all_dt, FUN=sum)
boxplot(lowSameDirFC~rD_range,data=agg2b_signif_dt)




agg2c_signif_dt <- aggregate(lowSameDirFCmissed ~ hicds+exprds+rD_range, data = all_dt, FUN=sum)
boxplot(lowSameDirFCmissed~rD_range,data=agg2c_signif_dt)


agg2d_dt <- merge(agg2c_signif_dt, agg2_signif_dt, by=c("hicds", "exprds", "rD_range"), all=TRUE)
agg2d_dt$ratioLowInMissed <- agg2d_dt$lowSameDirFCmissed/agg2d_dt$missGeneSignif
boxplot(ratioLowInMissed~rD_range,data=agg2d_dt)


agg2e_signif_dt <- aggregate(lowSameDirFC ~ hicds+exprds+rD_range+geneSignif, data = all_dt, FUN=sum)
boxplot(lowSameDirFC+geneSignif~rD_range,data=agg2e_signif_dt)

agg2e_signif_dt <- aggregate(lowSameDirFC ~ hicds+exprds+rD_range+geneSignif, data = all_dt, FUN=sum)
boxplot(lowSameDirFC+geneSignif~rD_range,data=agg2e_signif_dt)



agg3_signif_dt <- aggregate(missGeneSignif ~ hicds+exprds+pvalCorr_range, data = all_dt[all_dt$TAD_meanCorrPval <= 0.05,], FUN=sum)
boxplot(missGeneSignif~pvalCorr_range,data=agg3_signif_dt)

agg3c_signif_dt <- aggregate(lowSameDirFCmissed ~ hicds+exprds+pvalCorr_range, data = all_dt[all_dt$TAD_meanCorrPval <= 0.05,],  FUN=sum)
boxplot(lowSameDirFCmissed~pvalCorr_range,data=agg3c_signif_dt)

agg3d_dt <- merge(agg3c_signif_dt, agg3_signif_dt, by=c("hicds", "exprds", "pvalCorr_range"), all=TRUE)
agg3d_dt$ratioLowInMissed <- agg3d_dt$lowSameDirFCmissed/agg3d_dt$missGeneSignif
boxplot(ratioLowInMissed~pvalCorr_range,data=agg3d_dt)

missedGenes_dt <- all_dt[all_dt$missGeneSignif,]
agg4a_signif_dt <- aggregate(lowSameDirFC ~ hicds+exprds+rD_range, data = missedGenes_dt, FUN=sum)
boxplot(lowSameDirFC~rD_range, ylab=plot_var, cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, main = plotDesc, xlab="ratioDown range",data=agg4a_signif_dt)
agg4b_signif_dt <- aggregate(lowDiffDirFC ~ hicds+exprds+rD_range, data = missedGenes_dt, FUN=sum)
boxplot(lowDiffDirFC~rD_range,data=agg4b_signif_dt)

agg4 <- merge(agg4a_signif_dt, agg4b_signif_dt, by=c("hicds", "exprds", "rD_range"), all=T)
agg4$sameOverDiff <- agg4$lowSameDirFC/agg4$lowDiffDirFC
boxplot(lowDiffDirFC~rD_range,data=agg4b_signif_dt)boxplot(lowDiffDirFC~rD_range,data=agg4b_signif_dt)


    
# dt$tadSignif <- dt$tad_adjCombPval <= tadSignifThresh
    # dt$geneSignif <- dt$adj.P.Val <= geneSignifThresh
    # 
    # dt$divmaxFC <- dt$logFC/dt$TAD_maxFC
    # 
    # gene_missed_dt <- dt[!dt$geneSignif,]
    # x_var <- "logFC"
    # for(y_var in c("TAD_rD", "divmaxFC", "TAD_varFC", "TAD_logFCpval_log10","TAD_meanCorrPval_log10")) {
    #   plot(
    #     y=missed_dt[,paste0(y_var)],
    #     x=missed_dt[,paste0(x_var)],
    #     xlab=paste0(x_var),
    #     ylab=paste0(y_var),
    #     col=as.numeric(missed_dt$tadSignif)+1,
    #     pch=16
    #   )
    #   
    #   pre_legend_title <- expression(paste(bold("TAD-level")))
    #   legend("topright",
    #          title=pre_legend_title,
    #          
    #          legend=c( "not signif.", "signif."),
    #          pch=c(16,16), 
    #          col = c(1,2), 
    #          # text.width = c(2,1,1),
    #          bty="n")
    #   
    # }
    # 
    # 
    # 
    # densplot(
    #   y=missed_dt$TAD_rD,
    #   x=missed_dt$logFC,
    #   pch=16
    # )
    # 
    # 
    
    
    
