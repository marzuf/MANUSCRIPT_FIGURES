# Rscript fcc_ratioNbrTADs_heatmap_permG2T.R 

require(doMC)
require(foreach)
require(ggplot2)
require(ggsci)
require(colorRamps)

require(ggpubr)

registerDoMC(40)

plotType <- "svg"

source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")

a_t <- "permg2t"

myWidthGG <- 10
myHeightGG <- 6


outFolder <- "FCC_RATIONBRTADS_HEATMAP_PERMG2T"
dir.create(outFolder, recursive = TRUE)

buildData <- FALSE

fcc_fract <- seq(from=-1, to=1, by=0.1)
# fcc_fract_names <- paste0("FCC > ", fcc_fract[1:(length(fcc_fract)-1)], " and FCC <= ",fcc_fract[2:length(fcc_fract)])
fcc_fract_names <- paste0("FCC \u2208 ]", fcc_fract[1:(length(fcc_fract)-1)], ", ",fcc_fract[2:length(fcc_fract)], "]")
fcc_fract_names <- paste0("]", fcc_fract[1:(length(fcc_fract)-1)], ", ",fcc_fract[2:length(fcc_fract)], "]")

fcc_fract_names[fcc_fract_names == "]-1, -0.9]"] <- "[-1, -0.9]"


# fract_sort <- "FCC > 0.75 and FCC <= 1"
fract_sort <- fcc_fract_names[length(fcc_fract_names)]

ggsci_pal <- "lancet"
ggsci_subpal <- ""

legTitle <- "FCC ranges:"
fractBarSubTitle <- "AUC ratios:\n"
fractBarTitle <- "Fold-change concordance scores"

plotMargin <- c(1,2,1,1)

auc_ratio_file <- file.path("FCC_WAVE_PLOT_NOABS/all_fcc_dt.Rdata")
stopifnot(file.exists(auc_ratio_file))
x <- get(load(auc_ratio_file))
x$dataset <- file.path(x$hicds,  x$exprds)
x <- x[order(x$fcc_auc, decreasing=TRUE),]
fcc_ds_order <- x$dataset

dsCols <- all_cols[all_cmps[basename(fcc_ds_order)]]
stopifnot(!is.na(dsCols))

signif_file <- file.path(runFolder, "CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(signif_file))

all_hicds <- all_obs_hicds
all_exprds <- all_obs_exprds

# all_hicds=all_hicds[1]

keepPermut <- 1000

if(buildData){
  all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
    cat(paste0("... start: ", hicds, "\n"))
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      cat(paste0("... start: ", hicds," - ", exprds,  "\n"))

      fcc_file <- file.path(pipFolder, hicds, exprds, "8cOnlyFCC_runAllDown", "prodSignedRatio_permDT.Rdata")
      if(!file.exists(fcc_file)) return(NULL)
      stopifnot(file.exists(fcc_file))
      fcc_perm_dt <- get(load(fcc_file))
      stopifnot(ncol(fcc_perm_dt) >= keepPermut)
      
      keepCols <- sample(x=1:ncol(fcc_perm_dt), size = keepPermut)
      
      stopifnot(length(keepCols) == keepPermut)
      
      fcc_perm_dt <- fcc_perm_dt[,keepCols ]
      all_fcc <- as.numeric(fcc_perm_dt)
      
      # [1] -1.00 -0.75   -0.50       -0.25   0.00    0.25  0.50  0.75  1.00
      # [1]   0   0        10        411      571      355  97 350
      # > sum(all_fcc > 0.75 & all_fcc <= 1)
      # [1] 350
      # > sum(all_fcc > 0.5 & all_fcc <= 0.75)
      # [1] 97
      # sum(all_fcc > -0.5 & all_fcc <= -0.25)
      # [1] 10
      fcc_hist <- hist(all_fcc, breaks=fcc_fract)$counts
      names(fcc_hist) <- fcc_fract_names
      fcc_hist_nbr <- fcc_hist
      fcc_hist <- fcc_hist/length(all_fcc)
      stopifnot(sum(fcc_hist) == 1)
      
      data.frame(
        hicds = hicds,
        exprds = exprds,
        intervalFCC = names(fcc_hist),
        ratioFCC = as.numeric(fcc_hist),
        nFCC = as.numeric(fcc_hist_nbr),
        stringsAsFactors = FALSE
      )
    }
    exprds_dt
  }
  outFile <- file.path(outFolder, paste0("all_dt.Rdata"))
  auc_fract_file <- outFile
  save(all_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  inFile <- file.path(outFolder, paste0("all_dt.Rdata"))
  auc_fract_file <- inFile
  all_dt <- get(load(inFile))
 # load("FCC_RATIONBRTADS_HEATMAP_PERMG2T//all_dt.Rdata")
}

all_dt$dataset <- file.path(all_dt$hicds, all_dt$exprds)
all_dt$dataset <- factor(all_dt$dataset, levels = fcc_ds_order)
stopifnot(!is.na(all_dt$dataset))

nDS <- length(fcc_ds_order)

all_dt$intervalFCC <- factor(all_dt$intervalFCC, levels=fcc_fract_names)
all_dt$nFCC_log10 <- log10(all_dt$nFCC + 1)

curr_heat_theme <- theme(
  text = element_text(family=fontFamily),
  legend.title = element_text(face="bold"),
  axis.text.x = element_text(colour = dsCols, size=12),
  axis.text.y= element_text(colour = "black", size=12),
  axis.title.x = element_text(colour = "black", size=14, face="bold"),
  axis.title.y = element_text(colour = "black", size=14, face="bold"),
  plot.title = element_text(hjust=0.5, size=16, face="bold"),
  plot.subtitle = element_text(hjust=0.5, size=14, face="italic"),
  panel.background = element_rect(fill = "transparent")
  # legend.background =  element_rect()
) 

subTit <- "PERMG2T data"

densityRatio_plot <- ggplot(all_dt, aes(x = dataset, y = intervalFCC, fill = ratioFCC))+
  geom_tile() +
  ggtitle(paste0("FCC score distribution"),
          subtitle = subTit) +
  scale_x_discrete(name="Datasets ranked by decreasing AUC FCC ratio", labels = rep(labsymbol, nDS ), expand = c(0, 0))  +
  scale_y_discrete(name="FCC score",  expand = c(0, 0))+
  labs(fill = "Ratio of TADs")+
  # geom_vline(xintercept=seq(from=1.5, by=1, length.out = nDS-1), linetype=3) + 
  geom_vline(xintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  geom_hline(yintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  theme(axis.line=element_line())+
  scale_fill_gradientn(colours=colorRamps::matlab.like2(20))+
  curr_heat_theme


outFile <- file.path(outFolder, paste0("FCC_score_dist_allDS_", a_t, "_ratioTADs_heatmap.", plotType))
ggsave(densityRatio_plot, filename = outFile,  height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


cut_dt <- all_dt[as.character(all_dt$intervalFCC) %in% unique(as.character(all_dt$intervalFCC[all_dt$ratioFCC > 0])),]

densityRatio_plot_cut <- ggplot(cut_dt, aes(x = dataset, y = intervalFCC, fill = ratioFCC))+
  geom_tile() +
  ggtitle(paste0("FCC score distribution"),
          subtitle = subTit) +
  scale_x_discrete(name="Datasets ranked by decreasing AUC FCC ratio", labels = rep(labsymbol, nDS ), expand = c(0, 0))  +
  scale_y_discrete(name="FCC score",  expand = c(0, 0))+
  labs(fill = "Ratio of TADs")+
  curr_heat_theme+
  # geom_vline(xintercept=seq(from=1.5, by=1, length.out = nDS-1), linetype=3) +
  geom_vline(xintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  geom_hline(yintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  theme(axis.line=element_line())+
  scale_fill_gradientn(colours=colorRamps::matlab.like2(20))

outFile <- file.path(outFolder, paste0("FCC_score_dist_allDS_", a_t, "_ratioTADs_heatmap_cut.", plotType))
ggsave(densityRatio_plot_cut, filename = outFile,  height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


densityNbr_plot <- ggplot(all_dt, aes(x = dataset, y = intervalFCC, fill = nFCC))+
  geom_tile() +
  ggtitle(paste0("FCC score distribution"),
          subtitle = subTit) +
  scale_x_discrete(name="Datasets ranked by decreasing AUC FCC ratio", labels = rep(labsymbol, nDS ), expand = c(0, 0))  +
  scale_y_discrete(name="FCC score",  expand = c(0, 0))+
  labs(fill = "# TADs")+
  # geom_vline(xintercept=seq(from=1.5, by=1, length.out = nDS-1), linetype=3) + 
  geom_vline(xintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  geom_hline(yintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  theme(axis.line=element_line())+
  scale_fill_gradientn(colours=colorRamps::matlab.like2(20))+
  curr_heat_theme

outFile <- file.path(outFolder, paste0("FCC_score_dist_allDS_", a_t, "_nbrTADs_heatmap.", plotType))
ggsave(densityNbr_plot, filename = outFile,  height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



subTit2 <- paste0(subTit, " [log10(#TADs+1)]")

densityNbrLog_plot <- ggplot(all_dt, aes(x = dataset, y = intervalFCC, fill = nFCC_log10))+
  geom_tile() +
  ggtitle(paste0("FCC score distribution"),
          subtitle = subTit2) +
  scale_x_discrete(name="Datasets ranked by decreasing AUC FCC ratio", labels = rep(labsymbol, nDS ), expand = c(0, 0))  +
  scale_y_discrete(name="FCC score",  expand = c(0, 0))+
  labs(fill = "# TADs\n[log10]")+
  # geom_vline(xintercept=seq(from=1.5, by=1, length.out = nDS-1), linetype=3) + 
  geom_vline(xintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  geom_hline(yintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  theme(axis.line=element_line())+
  scale_fill_gradientn(colours=colorRamps::matlab.like2(20))+
  curr_heat_theme

outFile <- file.path(outFolder, paste0("FCC_score_dist_allDS_", a_t, "_nbrTADsLog_heatmap.", plotType))
ggsave(densityNbrLog_plot, filename = outFile,  height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

densityNbr_plot_cut <- ggplot(cut_dt, aes(x = dataset, y = intervalFCC, fill = nFCC))+
  geom_tile() +
  ggtitle(paste0("FCC score distribution"),
          subtitle = subTit) +
  scale_x_discrete(name="Datasets ranked by decreasing AUC FCC ratio", labels = rep(labsymbol, nDS ), expand = c(0, 0))  +
  scale_y_discrete(name="FCC score",  expand = c(0, 0))+
  labs(fill = "# TADs")+
  # geom_vline(xintercept=seq(from=1.5, by=1, length.out = nDS-1), linetype=3) + 
  geom_vline(xintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  geom_hline(yintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  theme(axis.line=element_line())+
  scale_fill_gradientn(colours=colorRamps::matlab.like2(20))+
  curr_heat_theme

outFile <- file.path(outFolder, paste0("FCC_score_dist_allDS_", a_t, "_nbrTADs_heatmap_cut.", plotType))
ggsave(densityNbr_plot_cut, filename = outFile,  height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



subTit2 <- paste0(subTit, " [log10(#TADs+1)]")

densityNbrLog_plot_cut <- ggplot(cut_dt, aes(x = dataset, y = intervalFCC, fill = nFCC_log10))+
  geom_tile() +
  ggtitle(paste0("FCC score distribution"),
          subtitle = subTit2) +
  scale_x_discrete(name="Datasets ranked by decreasing AUC FCC ratio", labels = rep(labsymbol, nDS ), expand = c(0, 0))  +
  scale_y_discrete(name="FCC score",  expand = c(0, 0))+
  labs(fill = "# TADs\n[log10]")+
  # geom_vline(xintercept=seq(from=1.5, by=1, length.out = nDS-1), linetype=3) + 
  geom_vline(xintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  geom_hline(yintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  theme(axis.line=element_line())+
  scale_fill_gradientn(colours=colorRamps::matlab.like2(20))+
  curr_heat_theme

outFile <- file.path(outFolder, paste0("FCC_score_dist_allDS_", a_t, "_nbrTADsLog_heatmap_cut.", plotType))
ggsave(densityNbrLog_plot_cut, filename = outFile,  height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))




