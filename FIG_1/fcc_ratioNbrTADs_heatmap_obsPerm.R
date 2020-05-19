# Rscript fcc_ratioNbrTADs_heatmap_obsPerm.R 

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

outFolder <- "FCC_RATIONBRTADS_HEATMAP_OBSPERM"
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

myWidthGG <- 10
myHeightGG <- 6


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

obs_dt <- get(load(file.path("FCC_RATIONBRTADS_HEATMAP/all_dt.Rdata")))
# colnames(obs_dt)[4:5] <- paste0(colnames(obs_dt)[4:5], "_obs")
perm_dt <- get(load(file.path("FCC_RATIONBRTADS_HEATMAP_PERMG2T//all_dt.Rdata")))
# colnames(perm_dt)[4:5] <- paste0(colnames(perm_dt)[4:5], "_perm")

obsPerm_dt <- merge(obs_dt, perm_dt, by=c("hicds", "exprds", "intervalFCC"), suffixes = c("_obs", "_perm"), all=TRUE)
stopifnot(!is.na(obsPerm_dt))

obsPerm_dt$intervalFCC <- factor(obsPerm_dt$intervalFCC, levels=fcc_fract_names)


obsPerm_dt$nFCC_perm_mean <- obsPerm_dt$nFCC_perm/keepPermut


obsPerm_dt$ratioFCC_obs_perm <- obsPerm_dt$ratioFCC_obs/obsPerm_dt$ratioFCC_perm

obsPerm_dt$nFCC_obs_perm <- obsPerm_dt$nFCC_obs/obsPerm_dt$nFCC_perm_mean

stopifnot(round(na.omit(obsPerm_dt$nFCC_obs_perm),4) ==  round(na.omit(obsPerm_dt$ratioFCC_obs_perm),4))

obsPerm_dt$dataset <- file.path(obsPerm_dt$hicds, obsPerm_dt$exprds)
nDS <- length(unique(obsPerm_dt$dataset))

obsPerm_dt$ratioFCC_obs_perm_log2 <- log2(obsPerm_dt$ratioFCC_obs_perm)

cut_dt <- obsPerm_dt[as.character(obsPerm_dt$intervalFCC) %in% unique(as.character(obsPerm_dt$intervalFCC[obsPerm_dt$ratioFCC_obs > 0])),]

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

subTit <- "OBS/PERMG2T"

densityRatio_plot <- ggplot(obsPerm_dt, aes(x = dataset, y = intervalFCC, fill = ratioFCC_obs_perm_log2))+
  geom_tile() +
  ggtitle(paste0("FCC score distribution"),
          subtitle = subTit) +
  scale_x_discrete(name="Datasets ranked by decreasing AUC FCC ratio", labels = rep(labsymbol, nDS ), expand = c(0, 0))  +
  scale_y_discrete(name="FCC score",  expand = c(0, 0))+
  labs(fill = "Log2 obs./perm.\n# TADs")+
  # geom_vline(xintercept=seq(from=1.5, by=1, length.out = nDS-1), linetype=3) + 
  geom_vline(xintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  geom_hline(yintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  theme(axis.line=element_line())+
  scale_fill_gradientn(colours=colorRamps::matlab.like2(20))+
  curr_heat_theme


outFile <- file.path(outFolder, paste0("FCC_score_dist_allDS_obs_perm_ratio_heatmap.", plotType))
ggsave(densityRatio_plot, filename = outFile,  height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



densityRatio_plot_cut <- ggplot(cut_dt, aes(x = dataset, y = intervalFCC, fill = ratioFCC_obs_perm_log2))+
  geom_tile() +
  ggtitle(paste0("FCC score distribution"),
          subtitle = subTit) +
  scale_x_discrete(name="Datasets ranked by decreasing AUC FCC ratio", labels = rep(labsymbol, nDS ), expand = c(0, 0))  +
  scale_y_discrete(name="FCC score",  expand = c(0, 0))+
  labs(fill = "Log2 obs./perm.\n# TADs")+
  curr_heat_theme+
  # geom_vline(xintercept=seq(from=1.5, by=1, length.out = nDS-1), linetype=3) +
  geom_vline(xintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  geom_hline(yintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  theme(axis.line=element_line())+
  scale_fill_gradientn(colours=colorRamps::matlab.like2(20))

outFile <- file.path(outFolder, paste0("FCC_score_dist_allDS_obs_perm_ratio_heatmap_cut.", plotType))
ggsave(densityRatio_plot_cut, filename = outFile,  height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

