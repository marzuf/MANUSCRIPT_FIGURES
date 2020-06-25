
# Rscript ratioDown_signif.R

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
require(reshape2)
require(ggpubr)
require(ggsci)

library(ggridges)

ggsci_pal <- "aaas"
ggsci_subpal <- ""


geneSignifThresh <- 0.01
tadSignifThresh <- 0.01

all_ranks_dt <- get(load(file.path(runFolder, "GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata")))

mycols <- pal_uchicago()(5)[4:5]

outFolder <- file.path("RATIODOWN_SIGNIF")
dir.create(outFolder, recursive=TRUE)

options(scipen=100)

startTime <- Sys.time()

myHeightGG <- 7
myWidthGG <- 9

buildData <- F

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
  # load("RATIODOWN_SIGNIF/all_dt.Rdata")
}

rdCoord_threshLow <- 0.2
rdCoord_threshHigh <- 1-rdCoord_threshLow

all_dt$ratioDownCoord0a <- all_dt$TAD_rD >= 0.8 | all_dt$TAD_rD <= 0.2
all_dt$ratioDownCoord0 <- all_dt$TAD_rD >= rdCoord_threshHigh | all_dt$TAD_rD <= rdCoord_threshLow
all_dt$ratioDownCoord <- abs(all_dt$TAD_rD - 0.5) >= (0.5-rdCoord_threshLow)
stopifnot(all_dt$ratioDownCoord == all_dt$ratioDownCoord0)
stopifnot(all_dt$ratioDownCoord == all_dt$ratioDownCoord0a)
all_dt$ratioDownCoord0 <- NULL
all_dt$ratioDownCoord0a <- NULL
all_dt$ratioDownCoord_lab <- ifelse(all_dt$ratioDownCoord, paste0("<= ", rdCoord_threshLow, " | >= ", rdCoord_threshHigh), "other")
all_dt$geneSignif <- all_dt$adj.P.Val <= geneSignifThresh
all_dt$tadSignif <- all_dt$tad_adjCombPval <= tadSignifThresh

all_dt$signif_lab <- ifelse(all_dt$geneSignif & all_dt$tadSignif, "both", 
                            ifelse(all_dt$geneSignif, "gene-level only",
                                   ifelse(all_dt$tadSignif, "TAD-level only", "never")))
signif_dt <- all_dt[all_dt$geneSignif | all_dt$tadSignif,]
stopifnot(!"never" %in% signif_dt$signif_lab)


totNbr_signif_dt <- aggregate(entrezID ~ ratioDownCoord_lab + signif_lab, data=signif_dt, FUN=length)
colnames(totNbr_signif_dt)[colnames(totNbr_signif_dt) == "entrezID"] <- "nbrGenes"

totNbr_signif_dt$nbrGenes_log10 <- log10(totNbr_signif_dt$nbrGenes)


nDS <- length(unique(file.path(all_dt$hicds, all_dt$exprds)))

plotTit <- paste0("Signif. genes and TAD ratio of down-reg. genes")
subTit <- paste0("TAD adj. p-val <= ", tadSignifThresh, "; gene adj. p-val <= ", geneSignifThresh, "; all datasets - n=", nDS)

nbr_p <- ggplot(totNbr_signif_dt, aes(x = signif_lab, fill = ratioDownCoord_lab, y = nbrGenes_log10)) +
  ggtitle(plotTit, subtitle = subTit)+
  geom_bar(stat="identity") +
  scale_y_continuous(name="# of genes [log10]", breaks = scales::pretty_breaks(n = 8), expand=c(0,0))+
  # annotation_logticks()+
  labs(fill = "TAD ratio\ndown-reg.\ngenes", x = "Genes signif.")+
  scale_fill_manual(values=mycols)+
  my_box_theme+ 
  theme(
    legend.text = element_text(size=12),
    axis.text.x = element_text(size=16),
    axis.title.x = element_text(size=16),
    axis.text.y = element_text(size=14),
    axis.title.y = element_text(size=16),
    axis.line = element_line()
  )
  
outFile <- file.path(outFolder, paste0("nbrGeneSignif_ratioDown_barplot.", plotType))
ggsave(plot = nbr_p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

ratio_totNbr_signif_dt <- do.call(rbind, by(totNbr_signif_dt, totNbr_signif_dt$signif_lab, function(x){
  x$nbrGenesRatio <- x$nbrGenes/sum(x$nbrGenes)
  x
}))

ratio_p <- ggplot(ratio_totNbr_signif_dt, aes(x = signif_lab, fill = ratioDownCoord_lab, y = nbrGenesRatio)) +
  ggtitle(plotTit, subtitle = subTit)+
  geom_bar(stat="identity") +
  scale_y_continuous(name="Ratio of genes", breaks = scales::pretty_breaks(n = 8), expand=c(0,0))+
  # annotation_logticks()+
  labs(fill = "TAD ratio\ndown-reg.\ngenes", x = "Genes signif.")+
  scale_fill_manual(values=mycols)+
  my_box_theme+ 
  theme(
    legend.text = element_text(size=12),
    axis.text.x = element_text(size=16),
    axis.title.x = element_text(size=16),
    axis.text.y = element_text(size=14),
    axis.title.y = element_text(size=16),
    axis.line = element_line()
  )



outFile <- file.path(outFolder, paste0("ratioGeneSignif_ratioDown_barplot.", plotType))
ggsave(plot = ratio_p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

ratioDown_ratioSignif_geneLevel_tadLevel_dt <- ratio_totNbr_signif_dt
saveFile <- file.path(outFolder, "supp_fig2A_ratioDown_ratioSignif_geneLevel_tadLevel_dt.Rdata")
save(ratioDown_ratioSignif_geneLevel_tadLevel_dt, file=saveFile, version=2)
cat(paste0("... written:" , saveFile, "\n"))



###### abs. FC distribution

all_dt$logFC_abs <- abs(all_dt$logFC)
fc_qt <- quantile(x=all_dt$logFC_abs, probs=0.95)

cut_dt <- all_dt[all_dt$logFC_abs <= fc_qt,]

cut_dt_signif <- cut_dt[cut_dt$signif_lab != "never",]

plotTit <- "Distribution of |logFC|"
subTit <- paste0("|logFC| <= 0.95 qt. (<=", round(fc_qt,2), ")")

dist_p <- ggplot(cut_dt_signif, aes(x=logFC_abs, fill=signif_lab, col=signif_lab)) + 
  ggtitle(plotTit, subtitle = subTit)+
  labs(fill="Gene signif.", x = paste0("|logFC|"), y ="Density") + 
  guides(color=FALSE)+
  # geom_histogram(aes(y=..density..), alpha=0.5, 
  #                position="identity")+
  geom_density(alpha=.2) +
  # geom_histogram(aes(y=..density..))+
  # geom_density(alpha=.3)+
  eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")"))) + 
  eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")"))) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8))+
  my_box_theme+
  theme(
    axis.line = element_line(),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    legend.text =element_text(size=12)
  )

outFile <- file.path(outFolder, paste0("absLogFC_signifType_density.", plotType))
ggsave(dist_p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

dist_hist_p <- ggplot(cut_dt_signif, aes(x=logFC_abs, fill=signif_lab, col=signif_lab)) + 
  ggtitle(plotTit, subtitle = subTit)+
  labs(fill="Gene signif.", x = paste0("|logFC|"), y ="Density") + 
  guides(color=FALSE)+
  geom_histogram(aes(y=..density..), alpha=0.5, 
                 position="identity")+
  geom_density(alpha=.2) +
  # geom_histogram(aes(y=..density..))+
  # geom_density(alpha=.3)+
  eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")"))) + 
  eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")"))) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  my_box_theme+
  theme(
    axis.line = element_line(),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    legend.text =element_text(size=12)
  )

outFile <- file.path(outFolder, paste0("absLogFC_signifType_density_withHist.", plotType))
ggsave(dist_hist_p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


ridge1_p <- ggplot(cut_dt_signif, aes(x = logFC_abs, y = signif_lab, fill = signif_lab)) +
  ggtitle(plotTit, subtitle = subTit)+
  labs(fill="Gene signif.", x = paste0("|logFC|"), y="") + 
  geom_density_ridges_gradient(scale = 0.95, size = 0.1, rel_min_height = 0.001) +
  eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")"))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  my_box_theme+
  theme(
    axis.line = element_line(),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    legend.text =element_text(size=12)
  )
outFile <- file.path(outFolder, paste0("absLogFC_signifType_density_ggridge_v1.", plotType))
ggsave(ridge1_p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


ridge2_p <- ggplot(cut_dt_signif, aes(x = logFC_abs, y = signif_lab, fill = signif_lab)) +
  ggtitle(plotTit, subtitle = subTit)+
  labs(fill="Gene signif.", x = paste0("|logFC|"), y="") + 
  geom_density_ridges_gradient(scale = 1.5, size = 0.1, rel_min_height = 0.001) +
  eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ", alpha=0.7)"))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  my_box_theme+
  theme(
    axis.line = element_line(),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    legend.text =element_text(size=12)
  )


outFile <- file.path(outFolder, paste0("absLogFC_signifType_density_ggridge_v2.", plotType))
ggsave(ridge2_p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

