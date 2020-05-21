
# Rscript check_fcc_geneSize_randommidposstrict.R 

require(doMC)
require(foreach)
require(ggplot2)
require(ggsci)
require(colorRamps)
require(reshape2)

require(ggpubr)

registerDoMC(40)

plotType <- "svg"

source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")

myWidthGG <- 7
myHeightGG <- 5

outFolder <- "CHECK_FCC_GENESIZE_RANDOMMIDPOSSTRICT"
dir.create(outFolder, recursive = TRUE)

buildData <- F

ggsci_pal <- "lancet"
ggsci_subpal <- ""

all_hicds <- all_obs_hicds
all_exprds <- all_obs_exprds


size_levels <- c("3",">3 & <=6",">6 & <=9",">9 & <=12",">12")
fcc_levels <- rev(c("1", ">=0.75 & <1",">=0.5 & <0.75",">=0.25 & <0.5","<0.25"))

# all_hicds=all_hicds[1]


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
      
      tad_sizeLabels <- ifelse(tad_size == 3, "3", 
                               ifelse(tad_size > 3 & tad_size <= 6, ">3 & <=6",
                                      ifelse(tad_size > 6 & tad_size <= 9, ">6 & <=9",
                                             ifelse(tad_size > 9 & tad_size <= 12, ">9 & <=12",
                                                    ifelse(tad_size > 12, ">12", NA)))))
      stopifnot(!is.na(tad_sizeLabels))
      stopifnot(length(tad_sizeLabels) == length(tad_size))
      
      
      fcc_file <- file.path(pipFolder, hicds, exprds, step8fcc_folder, "all_obs_prodSignedRatio.Rdata")
      stopifnot(file.exists(fcc_file))
      all_obs_fcc <- get(load(fcc_file))
      
      
      tad_fccLabels <- ifelse(all_obs_fcc == 1, "1", 
                               ifelse(all_obs_fcc >= 0.75 & all_obs_fcc <1, ">=0.75 & <1",
                                      ifelse(all_obs_fcc >= 0.5 & all_obs_fcc <0.75, ">=0.5 & <0.75",
                                             ifelse(all_obs_fcc >=0.25 & all_obs_fcc <0.5, ">=0.25 & <0.5",
                                                    ifelse(all_obs_fcc < 0.25, "<0.25", NA)))))
      stopifnot(!is.na(tad_fccLabels))
      stopifnot(length(tad_fccLabels) == length(all_obs_fcc))
      
      stopifnot(setequal(names(tad_fccLabels), names(tad_sizeLabels)))
      
      obsFCC_labels_dt <- data.frame(
        region = names(tad_fccLabels),
        tad_fccLabels = tad_fccLabels[names(tad_fccLabels)],
        tad_sizeLabels = tad_sizeLabels[names(tad_fccLabels)],
        stringsAsFactors = FALSE
      ) 
      stopifnot(!is.na(obsFCC_labels_dt))
      
      count_obsFCC_labels_dt <- aggregate(region ~ tad_fccLabels+tad_sizeLabels, data=obsFCC_labels_dt, FUN=length)
      colnames(count_obsFCC_labels_dt)[colnames(count_obsFCC_labels_dt)=="region"] <- "nTADs_obs"
      count_obsFCC_labels_dt$ratioTADs_obs <- count_obsFCC_labels_dt$nTADs_obs/length(tad_fccLabels)
      
      rd_hicds <- gsub("_40kb", "_RANDOMMIDPOSSTRICT_40kb",hicds)
      gene_file <- file.path(pipFolder, rd_hicds, exprds, step0_folder, "pipeline_geneList.Rdata")
      geneList <- get(load(gene_file))
      g2t_dt <- read.delim(file.path(runFolder, rd_hicds, "genes2tad", "all_genes_positions.txt"),
                           header=FALSE,col.names=c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
      g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
      stopifnot(geneList %in% g2t_dt$entrezID)
      g2t_dt <- g2t_dt[g2t_dt$entrezID %in% geneList,]
      tad_size <- setNames(as.numeric(table(g2t_dt$region)), names(table(g2t_dt$region)))
      
      tad_sizeLabels <- ifelse(tad_size == 3, "3", 
                               ifelse(tad_size > 3 & tad_size <= 6, ">3 & <=6",
                                      ifelse(tad_size > 6 & tad_size <= 9, ">6 & <=9",
                                             ifelse(tad_size > 9 & tad_size <= 12, ">9 & <=12",
                                                    ifelse(tad_size > 12, ">12", NA)))))
      stopifnot(!is.na(tad_sizeLabels))
      stopifnot(length(tad_sizeLabels) == length(tad_size))
      
      
      fcc_file <- file.path(pipFolder, rd_hicds, exprds, step8fcc_folder, "all_obs_prodSignedRatio.Rdata")
      stopifnot(file.exists(fcc_file))
      all_obs_fcc <- get(load(fcc_file))
      
      
      tad_fccLabels <- ifelse(all_obs_fcc == 1, "1", 
                              ifelse(all_obs_fcc >= 0.75 & all_obs_fcc <1, ">=0.75 & <1",
                                     ifelse(all_obs_fcc >= 0.5 & all_obs_fcc <0.75, ">=0.5 & <0.75",
                                            ifelse(all_obs_fcc >=0.25 & all_obs_fcc <0.5, ">=0.25 & <0.5",
                                                   ifelse(all_obs_fcc < 0.25, "<0.25", NA)))))
      stopifnot(!is.na(tad_fccLabels))
      stopifnot(length(tad_fccLabels) == length(all_obs_fcc))
      
      stopifnot(setequal(names(tad_fccLabels), names(tad_sizeLabels)))
      
      rdFCC_labels_dt <- data.frame(
        region = names(tad_fccLabels),
        tad_fccLabels = tad_fccLabels[names(tad_fccLabels)],
        tad_sizeLabels = tad_sizeLabels[names(tad_fccLabels)],
        stringsAsFactors = FALSE
      ) 
      stopifnot(!is.na(rdFCC_labels_dt))
      
      count_rdFCC_labels_dt <- aggregate(region ~ tad_fccLabels+tad_sizeLabels, data=rdFCC_labels_dt, FUN=length)
      colnames(count_rdFCC_labels_dt)[colnames(count_rdFCC_labels_dt)=="region"] <- "nTADs_rd"
      count_rdFCC_labels_dt$ratioTADs_rd <- count_rdFCC_labels_dt$nTADs_rd/length(tad_fccLabels)
   
              
      out_dt <- merge(count_obsFCC_labels_dt, count_rdFCC_labels_dt, by=c("tad_fccLabels", "tad_sizeLabels"), all=TRUE)  
        
      out_dt$hicds <- hicds
      out_dt$exprds <- exprds
        
        out_dt
        
        
   
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
  # load("FCC_RATIONBRTADS_HEATMAP_PERMG2T//all_dt.Rdata")
}



# obs_plot_dt <- aggregate(nTADs_obs~tad_fccLabels + tad_sizeLabels, data = all_dt, FUN=mean)
# colnames(obs_plot_dt)[colnames(obs_plot_dt) == "nTADs_obs"] <- "nTADs_mean"
# obs_plot_dt$dataType <- "obs."
# perm_plot_dt <- aggregate(nTADs_meanPerm~tad_fccLabels + tad_sizeLabels, data = all_dt, FUN=mean)
# colnames(perm_plot_dt)[colnames(perm_plot_dt) == "nTADs_meanPerm"] <- "nTADs_mean"
# perm_plot_dt$dataType <- "permut."
# plot_dt <- rbind(obs_plot_dt, perm_plot_dt)

obs_plot_dt <- aggregate(ratioTADs_obs~tad_fccLabels + tad_sizeLabels, data = all_dt, FUN=mean)
colnames(obs_plot_dt)[colnames(obs_plot_dt) == "ratioTADs_obs"] <- "ratioTADs_meanObs"
perm_plot_dt <- aggregate(ratioTADs_rd~tad_fccLabels + tad_sizeLabels, data = all_dt, FUN=mean)
colnames(perm_plot_dt)[colnames(perm_plot_dt) == "ratioTADs_rd"] <- "ratioTADs_meanRd"
plot_dt <- merge(obs_plot_dt, perm_plot_dt, by=c("tad_fccLabels", "tad_sizeLabels"))

plot_dt$tad_fccLabels <- factor(plot_dt$tad_fccLabels, levels = fcc_levels)
stopifnot(!is.na(plot_dt$tad_fccLabels))

plot_dt$tad_sizeLabels <- factor(plot_dt$tad_sizeLabels, levels = size_levels)
stopifnot(!is.na(plot_dt$tad_sizeLabels))

curr_heat_theme <- theme(
  text = element_text(family=fontFamily),
  legend.title = element_text(face="bold"),
  axis.text.x = element_text(colour = "black", size=12),
  axis.text.y= element_text(colour = "black", size=12),
  axis.title.x = element_text(colour = "black", size=14, face="bold"),
  axis.title.y = element_text(colour = "black", size=14, face="bold"),
  plot.title = element_text(hjust=0.5, size=16, face="bold"),
  plot.subtitle = element_text(hjust=0.5, size=14, face="italic"),
  panel.background = element_rect(fill = "transparent")
  # legend.background =  element_rect()
) 

myWidthGG <- myWidthGG*1.2

nDS <- length(unique(file.path(all_dt$exprds, all_dt$hicds)))

subTit <- paste0("mean all obs. datasets (n=", nDS, ")")

# plot_dt$ratioTADs_meanObs_lab <- ifelse(is.null(all_dt$ratioTADs_meanObs), "", round(plot_dt$ratioTADs_meanObs,2))
plot_dt$ratioTADs_meanObs_lab <- round(plot_dt$ratioTADs_meanObs,4)
plot_dt$ratioTADs_meanRd_lab <- round(plot_dt$ratioTADs_meanRd,4)

obs_p <- ggplot(plot_dt, aes(x = tad_fccLabels, y = tad_sizeLabels, fill = ratioTADs_meanObs, label=ratioTADs_meanObs_lab))+
  geom_tile() +
  ggtitle(paste0("TAD FCC and # genes"),
          subtitle = subTit) +
  scale_x_discrete(name="FCC range", expand = c(0, 0))  +
  scale_y_discrete(name="# genes/TAD", expand = c(0, 0))+
  labs(fill = "ratio TADs", x ="", y="")+
  # geom_vline(xintercept=seq(from=1.5, by=1, length.out = nDS-1), linetype=3) + 
  geom_vline(xintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  geom_hline(yintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  geom_text(fontface="bold")+
  theme(axis.line=element_line(),
        axis.text.x = element_text(size=10, hjust=0.5, vjust=1),
        axis.text.y = element_text(size=10, hjust=1, vjust=1)
        )+
  scale_fill_gradientn(colours=colorRamps::matlab.like2(20))+
  curr_heat_theme

outFile <- file.path(outFolder, paste0("meanObs_ratioTADs_nbrGenes_FCC_ranges_heatmap.", plotType))
ggsave(obs_p, filename=outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


subTit <- paste0("mean all randommidposstrict datasets")

permut_p <- ggplot(plot_dt, aes(x = tad_fccLabels, y = tad_sizeLabels, fill = ratioTADs_meanRd, label=ratioTADs_meanRd_lab))+
  geom_tile() +
  ggtitle(paste0("TAD FCC and # genes"),
          subtitle = subTit) +
  scale_x_discrete(name="FCC range", expand = c(0, 0))  +
  scale_y_discrete(name="# genes/TAD", expand = c(0, 0))+
  labs(fill = "ratio TADs", x ="", y="")+
  # geom_vline(xintercept=seq(from=1.5, by=1, length.out = nDS-1), linetype=3) + 
  geom_vline(xintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  geom_hline(yintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  geom_text(fontface="bold")+
  theme(axis.line=element_line(),
        axis.text.x = element_text(size=10, hjust=0.5, vjust=1),
        axis.text.y = element_text(size=10, hjust=1, vjust=1)
  )+
  scale_fill_gradientn(colours=colorRamps::matlab.like2(20))+
  curr_heat_theme

outFile <- file.path(outFolder, paste0("meanPermut_ratioTADs_nbrGenes_FCC_ranges_heatmap.", plotType))
ggsave(permut_p, filename=outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



subTit <- paste0("ratio mean obs./mean randommidposstrict")
plot_dt$obs_perm_ratio <- plot_dt$ratioTADs_meanObs/plot_dt$ratioTADs_meanRd
plot_dt$obs_perm_ratio_lab <- round(plot_dt$obs_perm_ratio, 4)

ratio_p <- ggplot(plot_dt, aes(x = tad_fccLabels, y = tad_sizeLabels, fill = obs_perm_ratio, label=obs_perm_ratio_lab))+
  geom_tile() +
  ggtitle(paste0("TAD FCC and # genes"),
          subtitle = subTit) +
  scale_x_discrete(name="FCC range", expand = c(0, 0))  +
  scale_y_discrete(name="# genes/TAD", expand = c(0, 0))+
  labs(fill = "ratio TADs", x ="", y="")+
  # geom_vline(xintercept=seq(from=1.5, by=1, length.out = nDS-1), linetype=3) + 
  geom_vline(xintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  geom_hline(yintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  geom_text(fontface="bold")+
  theme(axis.line=element_line(),
        axis.text.x = element_text(size=10, hjust=0.5, vjust=1),
        axis.text.y = element_text(size=10, hjust=1, vjust=1)
  )+
  scale_fill_gradientn(colours=colorRamps::matlab.like2(20))+
  curr_heat_theme

outFile <- file.path(outFolder, paste0("meanObsPermutRatio_ratioTADs_nbrGenes_FCC_ranges_heatmap.", plotType))
ggsave(ratio_p, filename=outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


subTit <- paste0("ratio mean obs./mean randommidposstrict [log2]")
plot_dt$obs_perm_ratioLog2 <- log2(plot_dt$ratioTADs_meanObs/plot_dt$ratioTADs_meanRd)
plot_dt$obs_perm_ratioLog2_lab <- round(plot_dt$obs_perm_ratioLog2, 4)

ratio_p <- ggplot(plot_dt, aes(x = tad_fccLabels, y = tad_sizeLabels, fill = obs_perm_ratioLog2, label=obs_perm_ratioLog2_lab))+
  geom_tile() +
  ggtitle(paste0("TAD FCC and # genes"),
          subtitle = subTit) +
  scale_x_discrete(name="FCC range", expand = c(0, 0))  +
  scale_y_discrete(name="# genes/TAD", expand = c(0, 0))+
  labs(fill = "ratio TADs", x ="", y="")+
  # geom_vline(xintercept=seq(from=1.5, by=1, length.out = nDS-1), linetype=3) + 
  geom_vline(xintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  geom_hline(yintercept=seq(from=0.5, by=1, length.out = nDS+1), linetype=1) + 
  geom_text(fontface="bold")+
  theme(axis.line=element_line(),
        axis.text.x = element_text(size=10, hjust=0.5, vjust=1),
        axis.text.y = element_text(size=10, hjust=1, vjust=1)
  )+
  scale_fill_gradientn(colours=colorRamps::matlab.like2(20))+
  curr_heat_theme

outFile <- file.path(outFolder, paste0("meanObsPermutRatioLog2_ratioTADs_nbrGenes_FCC_ranges_heatmap.", plotType))
ggsave(ratio_p, filename=outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



