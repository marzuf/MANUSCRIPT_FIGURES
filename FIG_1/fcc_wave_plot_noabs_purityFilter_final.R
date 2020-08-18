#!/usr/bin/Rscript

# 28.04.2020 => do not take the abs before the cumsum




outFolder <- file.path("FCC_WAVE_PLOT_NOABS_PURITYFILTER_FINAL")
dir.create(outFolder, recursive=TRUE)

options(scipen=100)

startTime <- Sys.time()

# Rscript fcc_wave_plot_noabs_purityFilter_final.R

plotType <- "svg"
source("../settings.R")
source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")

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
# set colors:
polygonPermutCol <- rgb(188/255,188/255,188/255, 0.3)
meanPermutCol <- rgb(135/255,135/255,135/255, 0.3)
# qt95PermutCol <- rgb(140/255,140/255,140/255, 0.3)
qt95PermutCol <- rgb(140/255,140/255,140/255)
lwdObs <- 1.2
colObs <- "darkred"
pointObsCol <- colObs

# all_obs_hicds=all_obs_hicds[1]

#all_obs_hicds="LG1_40kb"

purity_ds <- "aran"
pm <- "CPE"
transfExpr <- "log10"
tad_signifThresh <- 0.01
corrPurityQtThresh <- 0.05
purityFtd_col <- "darkgreen"
corMet <- "pearson"
purity_plot_name <- "Aran - CPE"

purity_file <- file.path(runFolder, "ALLTADS_AND_PURITY_FINAL", purity_ds,pm, transfExpr, "all_ds_corrPurity_dt.Rdata") # here _final INPUT
purityData <- get(load(purity_file))
agg_purity <- aggregate(purityCorr~dataset+region, FUN=mean, data=purityData)

result_file <- file.path(runFolder, "CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)

merge_dt <- merge(agg_purity, resultData, by=c("dataset", "region"))
merge_dt$signif <- merge_dt$adjPvalComb <= tad_signifThresh
merge_dt$signif_lab <- ifelse(merge_dt$signif, paste0("adj. p-val <=", tad_signifThresh), paste0("adj. p-val >", tad_signifThresh) )
purityCorrThresh <- as.numeric(quantile(merge_dt$purityCorr[!merge_dt$signif], probs = corrPurityQtThresh ))

subTit <- paste0(corMet, "'s corr.", " - ", purity_plot_name, " data")
plotTit <- paste0("Purity corr. distribution")
myx_lab <- paste0(transfExpr, " expr. and purity correlation (meanTAD)")
outFile <- file.path(outFolder, paste0("exprPurityCorr_meanTAD_signif_notSignif_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.5))
plot_multiDens(split(merge_dt$purityCorr, merge_dt$signif_lab),
               plotTit = plotTit, my_xlab = myx_lab)
lines(density(merge_dt$purityCorr), col="green")
abline(v=purityCorrThresh, col="blue")
legend("topleft", lty=1, lwd=2, col=c("green", "blue"), bty="n", legend=c("all", paste0(corrPurityQtThresh, "-qt non-signif. TADs\n(=", round(purityCorrThresh, 2), ")")))
mtext(side=3, text = subTit)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

merge_dt$region_id <- file.path(merge_dt$dataset, merge_dt$region)
TADs_to_discard <- merge_dt$region_id[merge_dt$purityCorr <=  purityCorrThresh]
TADs_to_keep <- merge_dt$region_id[merge_dt$purityCorr >  purityCorrThresh]
stopifnot(length(TADs_to_keep) + length(TADs_to_discard) == nrow(merge_dt))


all_fcc_dt <- foreach(hicds = all_obs_hicds, .combine='rbind') %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    hicds_dt <- foreach(exprds = all_obs_exprds[[paste0(hicds)]], .combine='rbind') %do% {

      ds_TADs_to_keep <- TADs_to_keep[dirname(TADs_to_keep) ==  file.path(hicds,exprds) ]
      if(length(ds_TADs_to_keep) == 0) return(NULL)
      
      
      # if(! (hicds == "ENCSR444WCZ_A549_40kb" & exprds == "TCGAluad_mutKRAS_mutEGFR")) return(NULL)
#if(exprds != "TCGAluad_mutKRAS_mutEGFR") return(NULL)

		stopifnot(hicds %in% names(hicds_names))
		stopifnot(exprds %in% names(exprds_names))

		hicds_lab <- hicds_names[paste0(hicds)]
		exprds_lab <- exprds_names[paste0(exprds)]
			  
      cat(paste0("> START ", hicds, " - ", exprds, "\n"))
      
      settingF <- file.path(setDir, "/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA", "PIPELINE", "INPUT_FILES", hicds, paste0("run_settings_", exprds, ".R"))
      cat(paste0(settingF, "\n"))
      stopifnot(file.exists(settingF))
      
      pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")
      
      
      script0_name <- "0_prepGeneData"
      script1_name <- "1_runGeneDE"
      script8_name <- "8cOnlyFCC_runAllDown"
      script_name <- "fcc_wave_plot_noabs.R"
      cat(paste0("> START ", script_name,  "\n"))
      
      source(settingF)

      pipLogFile <- ""
      
      ################################****************************************************************************************
      ####################################################### PREPARE INPUT
      ################################****************************************************************************************
      
      # INPUT DATA
      gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
      gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
      
      cat(paste0(pipOutFold, "\n"))
      
      # UPDATE SELECT THE GENES ACCORDING TO THE SETTINGS PREPARED IN 0_PREPGENEDATA
      initList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_geneList.Rdata"))))
      geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))
      txt <- paste0(toupper(script_name), "> Start with # genes: ", length(geneList), "/", length(initList), "\n")
      cat(txt)
      
      stopifnot(!any(duplicated(names(geneList))))
      
      gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% geneList,]
      geneNbr <- setNames(as.numeric(table(gene2tadDT$region)), names(table(gene2tadDT$region)))
      
      ### SET OUTPUT
      
      # if plotSeparated == TRUE => 1 plot per ratio, otherwise all on the same figure (#x2 plots)
      plotSeparated <- F
      
      # "permThresh" the quantile of permutations to take is loaded from main_settings.R
      
      ###############
      ##### retrieve the direction of up/down
      ###############
      # retrieve the direction of up/down
      DE_topTable <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_topTable.Rdata"))))
      DE_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_geneList.Rdata"))))
      exprDT <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_rnaseqDT.Rdata"))))
      # samp1 <- eval(parse(text=load(paste0(setDir, "/", sample1_file))))
      # samp2 <- eval(parse(text=load(paste0(setDir, "/", sample2_file))))
      samp1 <- eval(parse(text=load(paste0(sample1_file))))
      samp2 <- eval(parse(text=load(paste0(sample2_file))))
      DE_topTable <- DE_topTable[DE_topTable$genes %in% names(DE_geneList),]
      stopifnot(all(DE_topTable$genes %in% names(DE_geneList)))
      stopifnot(!any(duplicated(names(DE_geneList))))
      stopifnot(all(colnames(exprDT) %in% c(samp1, samp2)))
      stopifnot(all(samp1 %in% colnames(exprDT)))
      stopifnot(all(samp2 %in% colnames(exprDT)))
      maxDownGene <- DE_topTable$genes[which.min(DE_topTable$logFC)]
      stopifnot(maxDownGene %in% rownames(exprDT))
      mean_expr1 <- mean(unlist(c(exprDT[maxDownGene, samp1])), na.rm=T)
      mean_expr2 <- mean(unlist(c(exprDT[maxDownGene, samp2])), na.rm=T)
      
      if(mean_expr1 > mean_expr2) {
        subtitDir <- paste0("down: ", toupper(cond1), " > ", toupper(cond2))
      } else{
        subtitDir <- paste0("down: ", toupper(cond2), " > ", toupper(cond1))
      }
      
      
      
      ################################****************************************************************************************
      ####################################################### ITERATE OVER RATIOS TO PLOT
      ################################****************************************************************************************
      
      curr_ratio_type <- "prodSignedRatio"
      
      #for(curr_ratio_type in allDown_limited) {
      cat(paste0("*** START ", curr_ratio_type, "\n"))
      
      obs_curr_down <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/all_obs_", curr_ratio_type, ".Rdata"))))
      permut_currDown <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/", curr_ratio_type, "_permDT.Rdata"))))
      
      # ensure I used the same set of TADs for the permutation and for the calculation
      # (NB: would also be possible to filter the obs_curr_down, but not the permut_currDown)
      stopifnot(all(names(obs_curr_down) %in% rownames(permut_currDown)))
      stopifnot(all(rownames(permut_currDown) %in% names(obs_curr_down)))
      
      interReg <- intersect(names(obs_curr_down),rownames(permut_currDown) )
      
      ############################################################
      ############################################################ # filter the TADs and sort
      ############################################################
      filter_obs_curr_down <- sort(obs_curr_down[interReg], decreasing = T)
      
      filterPurityFtd_obs_curr_down <- filter_obs_curr_down[names(filter_obs_curr_down) %in% basename(ds_TADs_to_keep) ]
      
      filter_permut_currDown_unsort <- permut_currDown[interReg,]
      stopifnot(length(filter_obs_curr_down) == nrow(filter_permut_currDown_unsort))
      
      filter_permut_currDown <- apply(filter_permut_currDown_unsort, 2, sort, decreasing=T)
      rownames(filter_permut_currDown) <- NULL
      stopifnot(length(filter_obs_curr_down) == nrow(filter_permut_currDown_unsort))
      
      # FOR ratioDown => plot ratioConcord, departure from 0.5
      
      my_stat_curr_ratio <- "prodSignedRatio"
      # => raw (departure 0)
      # prodSignedRatio -> does not need to be transformed
      filter_obs_curr_down_half <- filter_obs_curr_down
      filterPurityFtd_obs_curr_down_half <- filterPurityFtd_obs_curr_down
      filter_permut_currDown_half <- filter_permut_currDown
      
      # PLOT THE 1ST PLOT
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_genomeWide_FCC_cumsum_obs_permut.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      
      par(bty="l")
      
      my_main <- paste0("Genome-wide intra-TAD fold-change concordance")
      my_sub <- paste0(hicds_lab, " - " , exprds_lab)
      my_xlab <- paste0("TADs ranked by FCC")
      my_ylab <- paste0("FCC cumulative sum")
      
      observ_vect= filter_obs_curr_down_half
      purityFtd_observ_vect = filterPurityFtd_obs_curr_down_half
      permut_DT=filter_permut_currDown_half
      my_stat = my_stat_curr_ratio
      
      
      observ_vect <- sort(observ_vect, decreasing = T)
      purityFtd_observ_vect <- sort(purityFtd_observ_vect, decreasing = T)
      permut_DT <- apply(permut_DT, 2, sort, decreasing=T)
      
      x_val <- c(1:length(observ_vect))
      purityFtd_x_val <- c(1:length(purityFtd_observ_vect))
      
      # CHANGED HERE 28.04.2020 -> do not take abs !
      # diff_05_permut <- apply(permut_DT, 2, function(x) cumsum(abs(x-0)))
      diff_05_permut <- apply(permut_DT, 2, function(x) cumsum(x-0))
      
      meanPermut_cumsum <- apply(diff_05_permut, 1, mean)
      qt95Permut_cumsum <- apply(diff_05_permut, 1, quantile, probs=0.95)
      
      # CHANGED HERE 28.04.2020 -> do not take abs !
      # obs_cumsum <- cumsum(abs(observ_vect - 0))
      obs_cumsum <- cumsum(observ_vect - 0)
      purityFtd_obs_cumsum <- cumsum(purityFtd_observ_vect - 0)

		par(mar = par()$mar + c(0,1,0,0))
      
      plot(obs_cumsum ~ x_val,
           main= my_main,
           type = "l",
           pch = 16, cex = 0.7,
           xlab= my_xlab, 
           ylab= my_ylab,
           cex.main = plotCex,
           cex.lab = plotCex,
           cex.axis = plotCex,
           col = colObs,
           axes=FALSE,
           lwd = lwdObs,
           bty="l")
      box()
      axis(2, cex.axis=plotCex, cex.lab=plotCex)
      axis(1, cex.axis=plotCex, cex.lab=plotCex, at = seq(from=0, to=2000, by=200))
      
      # add here the purity filter curve
      lines(
        x=purityFtd_x_val,
        y=purityFtd_obs_cumsum,
        col=purityFtd_col
      )
      legend("bottomright", legend="purity-filtered TADs", col=purityFtd_col, bty="n", lty=1)

            mtext(my_sub, font=3)
      
      polygon(x = c(x_val, rev(x_val)), 
              y = c( apply(diff_05_permut, 1, function(x) min(x)), rev(apply(diff_05_permut, 1, function(x) max(x)))),
              border=NA,
              col = polygonPermutCol)
      # lines(
      #   x = x_val,
      #   y= meanPermut_cumsum,
      #   col = meanPermutCol
      # )
      
      lines(
        x = x_val,
        y= qt95Permut_cumsum,
        col = qt95PermutCol
      )
      
      auc_obs <- auc(x = x_val, y = obs_cumsum)
      auc_permut <- auc(x = x_val, y = meanPermut_cumsum)
      auc_permutQt <- auc(x = x_val, y = qt95Permut_cumsum)
      
      stopifnot(length(purityFtd_x_val) == length(purityFtd_obs_cumsum))
      stopifnot(length(purityFtd_x_val) <= length(meanPermut_cumsum))
      stopifnot(length(purityFtd_x_val) <= length(qt95Permut_cumsum))
      
      auc_obsPF <- auc(x = purityFtd_x_val, y = purityFtd_obs_cumsum)
      auc_permutPF <- auc(x = purityFtd_x_val, y = meanPermut_cumsum[1:length(purityFtd_x_val)])
      auc_permutQtPF <- auc(x = purityFtd_x_val, y = qt95Permut_cumsum[1:length(purityFtd_x_val)])
      
# 
# saveFile <- file.path(outFolder, paste0("fig1F_", hicds, "_", exprds, "_", "fcc_cumsum_dt.Rdata"))      
# fcc_cumsum_dt <- data.frame(
# 			hicds=hicds,
# 			exprds=exprds,
# 			x_rank = x_val,
# 			obs_FCC = obs_cumsum,
# 			qt95Permut_FCC=qt95Permut_cumsum,
# 			stringsAsFactors=FALSE
# )
# save(fcc_cumsum_dt, file= saveFile)
# cat(paste0("... written: ", outFile, "\n"))
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", "auc_obs.Rdata"))
      save(auc_obs, file= outFile)
      cat(paste0("... written: ", outFile, "\n"))
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", "auc_permutQt.Rdata"))
      save(auc_permutQt, file= outFile)
      cat(paste0("... written: ", outFile, "\n"))
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", "auc_obsPF.Rdata"))
      save(auc_obsPF, file= outFile)
      cat(paste0("... written: ", outFile, "\n"))
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", "auc_permutQtPF.Rdata"))
      save(auc_permutQtPF, file= outFile)
      cat(paste0("... written: ", outFile, "\n"))
      
      
      
      # save(auc_obs, file= "auc_obs.Rdata")
      # save(auc_permut, file="auc_permut.Rdata")
      # save(x_val,file= "x_val.Rdata")
      # save(obs_cumsum, file="obs_cumsum.Rdata")
      # save(meanPermut_cumsum, file ="meanPermut_cumsum.Rdata")
      
      pct_inc <- round(auc_obs/auc_permut,2)
      pct_inc_qt <- round(auc_obs/auc_permutQt,2)
      
      # 
      # legend("topleft",
      #      xjust=0.5, yjust=0,
      #      pch = c(16, 15), 
      #      legend = c(paste0("observed (n=", length(observ_vect), ")"), "min-max permut"), 
      #      pt.cex = c(0.7, 2),
      #      col = c(pointObsCol, polygonPermutCol),
      #      bty="n")
      
      # legend("topleft",
      #        xjust=0.5, yjust=0,
      #        pch = c(16, 15, -1, -1 ),
      #        lty=c(-1,-1,1,1),
      #        legend = c(paste0("observed (n=", length(observ_vect), ")"), "min-max permut.", "avg. permut.", "0.95-qt permut."), 
      #        pt.cex = c(0.7, 2),
      #        col = c(pointObsCol, polygonPermutCol, meanPermutCol, qt95PermutCol),
      #        bty="n")
      # legtxt1 <- as.expression(bquote(frac(AUC[obs.], AUC[avg.permut.]) == .(pct_inc)))
      # legtxt2 <- as.expression(bquote(frac(AUC[obs.], AUC[qt.permut.]) == .(pct_inc_qt)))
      # legend("bottomright", legend=c(legtxt1, legtxt2), bty="n")
      legend("topleft",
             xjust=0.5, yjust=0,
             pch = c(16, 15, -1 ),
             lty=c(-1,-1,1),
             legend = c(paste0("observed (n=", length(observ_vect), ")"), "min-max permut.", "0.95-qt permut."), 
             pt.cex = c(0.7, 2),
             col = c(pointObsCol, polygonPermutCol, qt95PermutCol),
             bty="n")
      legtxt <- as.expression(bquote(frac(AUC[obs.], AUC[permut.]) == .(pct_inc_qt)))
      # legend("bottomright", legend=c(legtxt), bty="n")
      legend("right", legend=c(legtxt), bty="n")
      
      # text(x = length(observ_vect), y = 0+50,
      #      pos=2,
      #      bquote(frac(AUC[obs.], AUC[permut.]) == .(pct_inc)))
      
      
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      data.frame(
        hicds = hicds,
        exprds = exprds,
        fcc_auc = auc_obs/auc_permutQt ,
		fcc_auc_mean = auc_obs/auc_permut,
		
		fcc_auc_PF = auc_obsPF/auc_permutQtPF ,
		fcc_auc_mean_PF = auc_obsPF/auc_permutPF,
        stringsAsFactors=FALSE
      )
      
      
    }
    hicds_dt
}

outFile <- file.path(outFolder, "all_fcc_dt.Rdata")
save(all_fcc_dt, file = outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))



txt <- paste0(startTime, "\n", Sys.time(), "\n")
#printAndLog(txt, pipLogFile)
cat(txt)

cat(paste0("*** DONE: ", script_name, "\n"))


