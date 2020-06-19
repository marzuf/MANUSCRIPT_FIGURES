
# Rscript dbchart_fcc_obsPerm.R 

require(doMC)
require(foreach)
require(ggplot2)
require(ggsci)
require(colorRamps)
require(reshape2)

require(ggpubr)

registerDoMC(50)

plotType <- "svg"

source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")

myWidthGG <- 9
myHeightGG <- 6

outFolder <- "DBCHART_FCC_OBSPERM_1"
dir.create(outFolder, recursive = TRUE)

buildData <- FALSE



all_hicds <- all_obs_hicds
all_exprds <- all_obs_exprds

# all_hicds=all_hicds[1:2]

fcc_thresh <- 1

keepPermut <- 1000

if(buildData){
  all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
    cat(paste0("... start: ", hicds, "\n"))
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      cat(paste0("... start: ", hicds," - ", exprds,  "\n"))
      
      gene_file <- file.path(pipFolder, hicds, exprds, step0_folder, "pipeline_geneList.Rdata")
      geneList <- get(load(gene_file))
      g2t_dt <- read.delim(file.path(runFolder, hicds, "genes2tad", "all_genes_positions.txt"),
                           header=FALSE,col.names=c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
      g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
      stopifnot(geneList %in% g2t_dt$entrezID)
      g2t_dt <- g2t_dt[g2t_dt$entrezID %in% geneList,]
      tad_size <- setNames(as.numeric(table(g2t_dt$region)), names(table(g2t_dt$region)))
      
      fcc_file <- file.path(pipFolder, hicds, exprds, step8fcc_folder, "all_obs_prodSignedRatio.Rdata")
      stopifnot(file.exists(fcc_file))
      all_obs_fcc <- get(load(fcc_file))
      
      obs_aboveThresh <- all_obs_fcc[all_obs_fcc >= fcc_thresh]
      nObs_aboveThresh <- sum(all_obs_fcc >= fcc_thresh)
      ratioObs_aboveThresh <- mean(all_obs_fcc >= fcc_thresh)
      
      stopifnot(ratioObs_aboveThresh >= 0 & ratioObs_aboveThresh <= 1)
      
      fcc_file <- file.path(pipFolder, hicds, exprds, "8cOnlyFCC_runAllDown", "prodSignedRatio_permDT.Rdata")
      if(!file.exists(fcc_file)) return(NULL)
      stopifnot(file.exists(fcc_file))
      fcc_perm_dt <- get(load(fcc_file))
      stopifnot(ncol(fcc_perm_dt) >= keepPermut)
      
      keepCols <- sample(x=1:ncol(fcc_perm_dt), size = keepPermut)
      stopifnot(length(keepCols) == keepPermut)
      fcc_perm_dt <- fcc_perm_dt[,keepCols ]

      all_perm_nAboveThresh <- apply(fcc_perm_dt, 2, function(x) sum(x >= fcc_thresh))
      all_perm_ratioAboveThresh <- apply(fcc_perm_dt, 2, function(x) mean(x >= fcc_thresh))
      
      mean_perm_nAboveThresh <- mean(all_perm_nAboveThresh)
      mean_perm_ratioAboveThresh <- mean(all_perm_ratioAboveThresh)
      
      stopifnot(mean_perm_ratioAboveThresh >= 0 & mean_perm_ratioAboveThresh <= 1)
      
      data.frame(
        hicds = hicds,
        exprds = exprds,
        nObs_aboveThresh=nObs_aboveThresh,
        ratioObs_aboveThresh=ratioObs_aboveThresh,
        nMeanPerm_aboveThresh=mean_perm_nAboveThresh,
        ratioMeanPerm_aboveThresh=mean_perm_ratioAboveThresh,
        stringsAsFactors = FALSE
      )
    }
    exprds_dt
  }
  outFile <- file.path(outFolder, paste0("all_dt.Rdata"))
  # auc_fract_file <- outFile
  save(all_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  inFile <- file.path(outFolder, paste0("all_dt.Rdata"))
  # auc_fract_file <- inFile
  all_dt <- get(load(inFile))
  # load("DBCHART_FCC_OBSPERM/all_dt.Rdata")
}

library(ggcharts)

auc_ratio_file <- file.path("FCC_WAVE_PLOT_NOABS/all_fcc_dt.Rdata")
stopifnot(file.exists(auc_ratio_file))
x <- get(load(auc_ratio_file))
x$dataset <- file.path(x$hicds,  x$exprds)
x <- x[order(x$fcc_auc, decreasing=TRUE),]
fcc_ds_order <- x$dataset
all_dt$dataset <- file.path(all_dt$hicds, all_dt$exprds)
all_dt$dataset <- factor(all_dt$dataset, levels=fcc_ds_order)
stopifnot(!is.na(all_dt$dataset))

all_dt <- all_dt[order(as.numeric(all_dt$dataset)),]

# db_p <- 
#   dumbbell_chart(all_dt, dataset, ratioObs_aboveThresh, ratioMeanPerm_aboveThresh, sort=FALSE) + 
#   theme(
#     axis.text.y = element_blank(),
#     axis.line = element_line()
#   ) + 
#   labs(x="Datasets ranked by decreasing\nFCC AUC ratio (top=higher) ")

ggsci_pal <- "lancet"
ggsci_subpal <- ""

ggsci_pal <- "lancet"
ggsci_subpal <- ""

observ_col <- pal_lancet()(3)[1]
permut_col <- pal_lancet()(3)[2]
mycols <- setNames(c(observ_col, permut_col), c("observed", "permut."))

# bar_colors <-  c("steelblue3", "orangered")
bar_colors <-  pal_lancet()(2)
bar_colors <- c(observ_col, permut_col)

line_color <- "darkgrey"
line_size <- 1
bar_names <- setNames(c("observed", "permut."), bar_colors)

point_size <- 4
horizontal <- FALSE

if(fcc_thresh == 1 | fcc_thresh == -1) {
  plotTit <- paste0("Ratio of TADs with FCC = ", fcc_thresh)
} else {
  plotTit <- paste0("Ratio of TADs with FCC >= ", fcc_thresh)
}
subTit <- paste0("(mean ", keepPermut, " permut)")

nDS <- length(unique(all_dt$dataset))

myy_lab <- "Ratio of TADs"
myx_lab <- paste0("Datasets ranked by decreasing FCC AUC ratio (n=", nDS, ")")


dumbbell_p <- ggplot(all_dt, aes(x = dataset)) + 
  geom_segment(mapping = aes(xend = dataset, y = ratioMeanPerm_aboveThresh, yend = ratioObs_aboveThresh), 
               color = line_color, size = line_size) +
  geom_point(aes(y = ratioObs_aboveThresh, color = mycols["observed"]),
                                             size = point_size) + 
  geom_point(aes(y = ratioMeanPerm_aboveThresh , color = mycols["permut."]),
             size = point_size) + 
            scale_color_manual(values = bar_colors, labels=bar_names) + 
  
  # ggcharts:::ggcharts_current_theme(grid = ifelse(horizontal, "Y", "X"))+
  
      theme(legend.position = "top") + 
  labs(x=myx_lab,y=myy_lab, color = "Data:")+
  # guides(color = guide_legend(title = NULL)) +
  ggtitle(plotTit, subtitle = subTit)+
  expand_limits(x= c(-1, length(fcc_ds_order) + 1))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8))+
  # eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")"))) + 
  # eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")"))) + 
  # 
  my_box_theme+
  theme(
    text = element_text(family=fontFamily, color = "black"),
    # axis.text.x=element_text(angle=90)
    plot.title = element_text(size=18, face="bold", family=fontFamily, hjust=0.5),
    plot.subtitle = element_text(size=14, face="italic", family=fontFamily, hjust=0.5),
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    axis.text.x=element_blank(),
    axis.ticks.x = element_blank(),
    axis.title=element_text(size=14),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14),
    legend.key = element_rect(fill = NA)
  ) + 
  theme(
    panel.background = element_blank(),
    plot.background = element_blank())


outFile <- file.path(outFolder, paste0("ratio_", fcc_thresh, "FCC_obsPerm_dumbbell_plot.", plotType))
ggsave(dumbbell_p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))
# 




# if (horizontal) {
#   plot <- plot + coord_flip()
# }




# 
# # load("CHECK_FCC1/all_dt.Rdata")
# all_dt$obs_ratioSize3 <- all_dt$obs_nFCC1size3/all_dt$obs_nFCC1
# all_dt$permMean_ratioSize3 <- all_dt$permMean_nFCC1size3/all_dt$permMean_nFCC1
# 
# m_all_dt <- melt(all_dt, id=c("hicds", "exprds"))
# m_all_dt$varType <- gsub("(.+?)_.+", "\\1", m_all_dt$variable)
# m_all_dt$varLab <- gsub("(.+?)_(.+)", "\\2", m_all_dt$variable)
# 
# m_all_dt$cmpType <- all_cmps[paste0(m_all_dt$exprds)]
# 
# nDS <- length(unique(file.path(m_all_dt$hicds, m_all_dt$exprds)))
# 
# plotTit <- "# TADs with FCC=1"
# subTit <- paste0("all datasets (n=", nDS, "); mean permut (# perm. = ", keepPermut, ")")
# 
# sub1_dt <- m_all_dt[grepl("nFCC", m_all_dt$varLab),]
# 
# plot_myBox <- function(p){
#   p <- p +   
#     geom_boxplot(notch = TRUE, outlier.shape=NA)+
#     geom_point(position=position_jitterdodge(),  alpha=0.5) +
#     eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")"))) + 
#     eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")"))) + 
#     my_box_theme+
#     theme(
#       legend.text=element_text(size=12),
#       axis.line=element_line(),
#       axis.text.x = element_text(size=14)
#     )
#   return(p)
# }
# # ggplot(sub1_dt, aes(x=varLab, color=varType, y=value)) +
# #   ggtitle(plotTit, subtitle = subTit)+
# #   labs(fill ="", color="", x="", y="# TADs" )+
# #   geom_boxplot(notch = TRUE, outlier.shape=NA) +
# #   geom_point(aes(color=cmpType),position=position_jitterdodge(),  alpha=0.5) +
# #   eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")"))) +
# #   eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")"))) +
# #   my_box_theme+
# #   theme(
# #     legend.text=element_text(size=12),
# #     axis.line=element_line(),
# #     axis.text.x = element_text(size=14)
# 
# 
# p1 <- plot_myBox(ggplot(sub1_dt, aes(x=varLab, color=varType, y=value)) + 
#   ggtitle(plotTit, subtitle = subTit)+
#   labs(fill ="", color="", x="", y="# TADs" ))
#   
# outFile <- file.path(outFolder, paste0("nFCC1_nFCC1size3_obsPerm_boxplot.", plotType))
# ggsave(p1, filename = outFile, height=myHeightGG, width=myWidthGG)
# cat(paste0("... written: ", outFile, "\n"))
# 
# 
# sub2_dt <- m_all_dt[!grepl("nFCC", m_all_dt$varLab) & grepl("ratio", m_all_dt$varLab),]
# 
# plotTit <- "Ratio of FCC=1 TADs of size # genes=3"
# 
# p2 <- plot_myBox(ggplot(sub2_dt, aes(x=varLab, color=varType, y=value)) + 
#   ggtitle(plotTit, subtitle = subTit)+
#     labs(fill ="", color="", x="", y="Ratio of FCC=1 TADs" ))
# 
# outFile <- file.path(outFolder, paste0("FCC1ratioSize3_obsPerm_boxplot.", plotType))
# ggsave(p2, filename = outFile, height=myHeightGG, width=myWidthGG)
# cat(paste0("... written: ", outFile, "\n"))
# 
# 
# sub3_dt <- m_all_dt[!grepl("nFCC", m_all_dt$varLab) & ! grepl("ratio", m_all_dt$varLab),]
# 
# plotTit <- "Mean and median # genes of FCC=1 TADs"
# 
# p3 <- plot_myBox(ggplot(sub3_dt, aes(x=varLab, color=varType, y=value)) + 
#   ggtitle(plotTit, subtitle = subTit)+
#     labs(fill ="", color="", x="", y="Ratio of FCC=1 TADs" ))
# 
# outFile <- file.path(outFolder, paste0("FCC1meanMedianSize_obsPerm_boxplot.", plotType))
# ggsave(p3, filename = outFile, height=myHeightGG, width=myWidthGG)
# cat(paste0("... written: ", outFile, "\n"))
# 
# 
# 
