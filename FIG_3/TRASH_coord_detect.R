



# Rscript coord_detect.R


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


outFolder <- file.path("COORD_DETECT")
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

stop("-ok")

all_dt$geneSignif <- all_dt$adj.P.Val <= geneSignifThresh
all_dt$tadSignif <- all_dt$tad_adjCombPval <= tadSignifThresh
all_dt$sameDirFC <- sign(all_dt$logFC) == sign(all_dt$TAD_meanFC)
all_dt$lowFC <- abs(all_dt$logFC) <= lowFC_thresh
all_dt$lowSameDirFC <- all_dt$sameDirFC & all_dt$lowFC
all_dt$lowDiffDirFC <- !all_dt$sameDirFC & all_dt$lowFC

all_dt$missGeneSignif <- !all_dt$geneSignif & all_dt$tadSignif

all_dt$lowSameDirFCmissed <- all_dt$lowSameDirFC & all_dt$missGeneSignif

all_dt$lowDiffDirFCmissed <- all_dt$lowDiffDirFC & all_dt$missGeneSignif


rd_ranges <- seq(0,1,by=0.1)
pval_ranges <- c(seq(0,0.05,by=0.01), 1)

all_dt$rD_range <- cut(all_dt$TAD_rD, breaks=rd_ranges, include.lowest=TRUE)
stopifnot(!is.na(all_dt$rD_range))

all_dt$pvalCorr_range <- cut(all_dt$TAD_meanCorrPval, breaks=pval_ranges, include.lowest=TRUE)
stopifnot(!is.na(all_dt$pvalCorr_range))


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
boxplot(lowSameDirFC~rD_range,data=agg4a_signif_dt)
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
    
    
    
