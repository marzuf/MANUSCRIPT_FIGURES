# Rscript fcc_ratioTADs_heatmap.R 

require(doMC)
require(foreach)
require(ggplot2)
require(ggsci)

require(ggpubr)

registerDoMC(40)

plotType <- "svg"

source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")


outFolder <- "FCC_RATIOTADS_HEATMAP"
dir.create(outFolder, recursive = TRUE)

buildData <- TRUE

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
x$dataset <- paste0(x$hicds, "\n", x$exprds)
x <- x[order(x$fcc_auc, decreasing=TRUE),]
fcc_ds_order <- x$dataset
fcc_auc_sort <- TRUE

signif_file <- file.path(runFolder, "CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(signif_file))

all_hicds <- all_obs_hicds
all_exprds <- all_obs_exprds

if(buildData){
  all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
    cat(paste0("... start: ", hicds, "\n"))
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      cat(paste0("... start: ", hicds," - ", exprds,  "\n"))
      fcc_file <- file.path(pipFolder, hicds, exprds, step8fcc_folder, "all_obs_prodSignedRatio.Rdata")
      stopifnot(file.exists(fcc_file))
      all_fcc <- get(load(fcc_file))
      
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
  # load("FCC_RATIOTADS_HEATMAP/all_dt.Rdata")
}

all_dt$dataset <- file.path(all_dt$hicds, all_dt$exprds)
all_dt$intervalFCC <- factor(all_dt$intervalFCC, levels=fcc_fract_names)

density_plot <- ggplot(all_dt, aes(x = dataset, y = intervalFCC, fill = ratioFCC))+
  geom_tile() +
  ggtitle(paste0("FCC score distribution - ", a_t),
          subtitle = subTit) +
  scale_x_discrete(name="Datasets ranked by decreasing AUC FCC ratio", labels = rep(labsymbol, nDS ), expand = c(0, 0))  +
  scale_y_discrete(name="FCC score",  expand = c(0, 0))+
  labs(fill = "Density")+
  theme(
    text = element_text(family=fontFamily),
    axis.text.x = element_text(colour = dsCols, size=12),
    axis.text.y= element_text(colour = "black", size=12),
    axis.title.x = element_text(colour = "black", size=14, face="bold"),
    axis.title.y = element_text(colour = "black", size=14, face="bold"),
    plot.title = element_text(hjust=0.5, size=16, face="bold"),
    plot.subtitle = element_text(hjust=0.5, size=14, face="italic"),
    panel.background = element_rect(fill = "transparent")
    # legend.background =  element_rect()
  ) + 
  geom_vline(xintercept=seq(from=1.5, by=1, length.out = nDS-1), linetype=3) +
    scale_fill_viridis_c(option="A")  

density_plot1 <- density_plot +
  # scale_fill_gradient( high="red", low="blue", na.value = "white" )  +
  scale_fill_gradient2( high="red", low="blue", na.value = "grey", mid ="white" , midpoint=mean(all_dt$ratioFCC))  


density_plot2 <- density_plot +
  scale_fill_viridis_c(option="A")  

outFile <- file.path(outFolder, paste0("FCC_score_dist_allDS_", a_t, "_densityheatmap.", plotType))
ggsave(density_plot1, filename = outFile,  height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



