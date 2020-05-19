# Rscript cmp_nSignif_sampleSize.R


setDir <- "/media/electron"
setDir <- ""

require(foreach)
require(doMC)
registerDoMC(40)
require(ggplot2)
require(ggpubr)
require(ggsci)
require(reshape2)

plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 9

sampleSize_datasets <- c(file.path("GSE105381_HepG2_40kb", "TCGAlihc_wt_mutCTNNB1"),
                         file.path("ENCSR312KHQ_SK-MEL-5_40kb", "TCGAskcm_lowInf_highInf"),
                         file.path("ENCSR079VIJ_G401_40kb","TCGAkich_norm_kich"))


nSubsamp <- 5

ggsci_pal <- "aaas"
ggsci_subpal <- ""


tadSignifThresh <- 0.01
geneSignifThresh <- 0.01

init_hicds <- "GSE105381_HepG2_40kb"
init_exprds <- "TCGAlihc_wt_mutCTNNB1"

source("../settings.R")

myWidthGG <- myWidthGG*1.2
myWidthGG <- 7
myHeightGG <- 5

myWidthGG <- myWidthGG*1.2 


outFolder <- file.path("CMP_NSIGNIF_SAMPLESIZE")
dir.create(outFolder, recursive = TRUE)


script1_name <- "1_runGeneDE"
script11_name <- "11sameNbr_runEmpPvalCombined"


all_tadSignif_dt <- foreach(ds = sampleSize_datasets, .combine='rbind') %do% {
  init_hicds <- dirname(ds)
  init_exprds <- basename(ds)
  stopifnot(init_hicds %in% names(hicds_names))
  stopifnot(init_exprds %in% names(exprds_names))
  stopifnot(init_hicds %in% all_obs_hicds)
  randomsub_hicds <- all_hicds[grepl("RANDOMSUB", all_hicds) & grepl(gsub(paste0("_40kb"), "", init_hicds), all_hicds)]
  my_hicds <- c(init_hicds, randomsub_hicds)
  my_exprds <- get_exprds(my_hicds)
  ds_tadSignif_dt <- foreach(hicds = my_hicds, .combine='rbind') %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    hicds_dt <- foreach(exprds = my_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      if(init_exprds != exprds) return(NULL)
      settingF <- file.path(runFolder, "PIPELINE", "INPUT_FILES", hicds, paste0("run_settings_", exprds, ".R"))
      stopifnot(file.exists(settingF))
      source(settingF)  
      s1 <- get(load(file.path(setDir, sample1_file)))
      s2 <- get(load(file.path(setDir, sample2_file)))
      # stopifnot(length(s1) == length(s2))
      # if(nsub != "") stopifnot(length(s1) == nsub)
      pval_file <- file.path(pipFolder,  hicds, exprds, script11_name, "emp_pval_combined.Rdata")
      #        if(!file.exists(pval_file)) return(NULL)
      stopifnot(file.exists(pval_file))
      pval <- get(load(pval_file))
      adj_pval <- p.adjust(pval, method="BH")
      data.frame(
        ds = ds,
        hicds = hicds,
        exprds = exprds,
        nSamp1 = length(s1),
        nSamp2 = length(s2),
        region = names(adj_pval),
        adjCombPval = as.numeric(adj_pval),
        stringsAsFactors = FALSE
      )
    }
    hicds_dt
  }
  ds_tadSignif_dt
}

outFile <- file.path(outFolder, "all_tadSignif_dt.Rdata")
save(all_tadSignif_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

all_tadSignif_dt <- do.call(rbind, by(all_tadSignif_dt, all_tadSignif_dt$ds, function(x){
  x$nSamp1_ratio <- x$nSamp1/max(x$nSamp1)
  x$nSamp2_ratio <- x$nSamp2/max(x$nSamp2)
  x
  }))
stopifnot(all_tadSignif_dt$nSamp1_ratio == all_tadSignif_dt$nSamp2_ratio)
all_tadSignif_dt$nSamp_ratio <- all_tadSignif_dt$nSamp1_ratio

tadSignif_agg_dt <- aggregate(adjCombPval~ds+hicds+exprds+nSamp_ratio,data=all_tadSignif_dt,FUN=function(x) mean(x <= tadSignifThresh))
colnames(tadSignif_agg_dt)[colnames(tadSignif_agg_dt) == "adjCombPval"] <- "ratioSignifTADs"

outFile <- file.path(outFolder, "tadSignif_agg_dt.Rdata")
save(tadSignif_agg_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

all_geneSignif_dt <- foreach(ds = sampleSize_datasets, .combine='rbind') %do% {
  init_hicds <- dirname(ds)
  init_exprds <- basename(ds)
  stopifnot(init_hicds %in% names(hicds_names))
  stopifnot(init_exprds %in% names(exprds_names))
  stopifnot(init_hicds %in% all_obs_hicds)
  randomsub_hicds <- all_hicds[grepl("RANDOMSUB", all_hicds) & grepl(gsub(paste0("_40kb"), "", init_hicds), all_hicds)]
  my_hicds <- c(init_hicds, randomsub_hicds)
  my_exprds <- get_exprds(my_hicds)
  ds_geneSignif_dt <- foreach(hicds = my_hicds, .combine='rbind') %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    hicds_dt <- foreach(exprds = my_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      #
      if(init_exprds != exprds) return(NULL)
      settingF <- file.path(runFolder, "PIPELINE", "INPUT_FILES", hicds, paste0("run_settings_", exprds, ".R"))
      stopifnot(file.exists(settingF))
      source(settingF)  
      s1 <- get(load(file.path(setDir, sample1_file)))
      s2 <- get(load(file.path(setDir, sample2_file)))
      pval_file <- file.path(pipFolder,  hicds, exprds, script1_name, "DE_topTable.Rdata")
  #    if(!file.exists(pval_file)) return(NULL)
      stopifnot(file.exists(pval_file))
      de_dt <- get(load(pval_file))
      data.frame(
        ds=ds,
        hicds = hicds,
        exprds = exprds,
        nSamp1 = length(s1),
        nSamp2 = length(s2),
        gene = as.character(de_dt$genes) ,
        adjPval = de_dt$adj.P.Val,
        stringsAsFactors = FALSE
      )
    }
    hicds_dt
  }
  ds_geneSignif_dt
}

outFile <- file.path(outFolder, "all_geneSignif_dt.Rdata")
save(all_geneSignif_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))
# all_geneSignif_dt <- get(load(file.path(outFolder, "all_geneSignif_dt.Rdata")))

all_geneSignif_dt <- do.call(rbind, by(all_geneSignif_dt, all_geneSignif_dt$ds, function(x){
  x$nSamp1_ratio <- x$nSamp1/max(x$nSamp1)
  x$nSamp2_ratio <- x$nSamp2/max(x$nSamp2)
  x
}))
stopifnot(all_geneSignif_dt$nSamp1_ratio == all_geneSignif_dt$nSamp2_ratio)
all_geneSignif_dt$nSamp_ratio <- all_geneSignif_dt$nSamp1_ratio

geneSignif_agg_dt <- aggregate(adjPval~ds+hicds+exprds+nSamp_ratio,data=all_geneSignif_dt,FUN=function(x) mean(x <= geneSignifThresh))
colnames(geneSignif_agg_dt)[colnames(geneSignif_agg_dt) == "adjPval"] <- "ratioSignifGenes"

outFile <- file.path(outFolder, "geneSignif_agg_dt.Rdata")
save(geneSignif_agg_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

allSignif_agg_dt <- merge(geneSignif_agg_dt, tadSignif_agg_dt, by=c("ds", "hicds", "exprds", "nSamp_ratio"), all=TRUE)
stopifnot(!is.na(allSignif_agg_dt))

save(allSignif_agg_dt, file=file.path(outFolder, "allSignif_agg_dt.Rdata"), version=2)
# load("CMP_NSIGNIF_SAMPLESIZE/allSignif_agg_dt.Rdata")

# just to the check the data correctly aggregated
tmp_dt <- allSignif_agg_dt[,c("ds", "hicds", "nSamp_ratio")]
check_dt <- aggregate(hicds~ds+nSamp_ratio, data=tmp_dt, FUN=length)
stopifnot(sum(check_dt$hicds) == length(sampleSize_datasets) + nSubsamp*(nrow(check_dt)-3))

allSignif_agg_dt_m <- melt(allSignif_agg_dt, id=c("ds", "hicds", "exprds", "nSamp_ratio" ))
mean_plot_dt <- aggregate(value~ds + nSamp_ratio + variable, data=allSignif_agg_dt_m, FUN=mean)
mean_plot_dt$variable <- as.character(mean_plot_dt$variable)
mean_plot_dt$variable[mean_plot_dt$variable == "ratioSignifGenes"] <- paste0("genes (adj. p-val <= ", geneSignifThresh, ")") 
mean_plot_dt$variable[mean_plot_dt$variable == "ratioSignifTADs"] <- paste0("TADs (adj. p-val <= ", tadSignifThresh, ")") 

mean_plot_dt$ds_lab <- paste0(hicds_names[dirname(mean_plot_dt$ds)], "\n", 
                              exprds_names[basename(mean_plot_dt$ds)])

stopifnot(!is.na(mean_plot_dt))

customSignif_p <- function(p) {
  p <- p + 
    geom_line() +  
    geom_point()+
    eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")"))) +
    my_box_theme +
    guides(color=guide_legend(keyheight=0.5,
                             default.unit="inch"))+
    theme(
      axis.line = element_line(),
      legend.text=element_text(size=12),
      legend.title=element_text(size=14),
      legend.key = element_rect(fill = NA))+
    labs(color="Datasets", linetype="Ratio signif.", x="Subsampling ratio")
  return(p)
}


##################################
## Plot both TADs and genes 
##################################

subTit <- paste0("average of ", nSubsamp, " subsamplings")
plotTit <- paste0("Ratio of features detected signif.")

both_p <- customSignif_p(
  ggplot(mean_plot_dt, aes(x=nSamp_ratio, y =value, color=ds_lab, linetype=variable)) + 
    scale_y_continuous(name="Ratio of signif. features", breaks = scales::pretty_breaks(n = 8))+
  ggtitle(plotTit, subtitle=subTit))

outFile <- file.path(outFolder, paste0("allDS_ratioSignif_genesAndTADs_linePlot.", plotType))
ggsave(plot = both_p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

##################################
## Plot only TADs  
##################################

subTit <- paste0("adj. p-val <= ", tadSignifThresh, "; average of ", nSubsamp, " subsamplings")
plotTit <- paste0("Ratio of TADs detected signif.")

onlyTADs_p <- customSignif_p(ggplot(mean_plot_dt[grepl("TADs", mean_plot_dt$variable),], aes(x=nSamp_ratio, y =value, color=ds_lab)) +
                               scale_y_continuous(name="Ratio of signif. TADs", breaks = scales::pretty_breaks(n = 8))+
                               ggtitle(plotTit, subtitle=subTit))

outFile <- file.path(outFolder, paste0("allDS_ratioSignif_onlyTADs_linePlot.", plotType))
ggsave(plot = onlyTADs_p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

##################################
## Plot only genes 
##################################

subTit <- paste0("adj. p-val <= ", geneSignifThresh, "; average of ", nSubsamp, " subsamplings")
plotTit <- "Ratio of genes detected signif."

onlyGenes_p <- customSignif_p(ggplot(mean_plot_dt[grepl("genes", mean_plot_dt$variable),], aes(x=nSamp_ratio, y =value, color=ds_lab)) + 
                                scale_y_continuous(name="Ratio of signif. genes", breaks = scales::pretty_breaks(n = 8))+
                                ggtitle(plotTit, subtitle=subTit))

outFile <- file.path(outFolder, paste0("allDS_ratioSignif_onlyGenes_linePlot.", plotType))
ggsave(plot = onlyGenes_p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


# wrap plot
wrap_both <- ggarrange(onlyTADs_p, onlyGenes_p, common.legend = TRUE, legend="right")

outFile <- file.path(outFolder, paste0("allDS_ratioSignif_wrapBoth_linePlot.", plotType))
ggsave(plot = wrap_both, filename = outFile, height=myHeightGG, width = myWidthGG*1.8)
cat(paste0("... written: ", outFile, "\n"))




