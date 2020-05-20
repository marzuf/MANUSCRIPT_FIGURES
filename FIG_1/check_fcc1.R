
# Rscript check_fcc1.R 

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

outFolder <- "CHECK_FCC1"
dir.create(outFolder, recursive = TRUE)

buildData <- FALSE

ggsci_pal <- "lancet"
ggsci_subpal <- ""

all_hicds <- all_obs_hicds
all_exprds <- all_obs_exprds

# all_hicds=all_hicds[1]

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
      all_obs_fcc1 <- all_obs_fcc[all_obs_fcc == 1]
      stopifnot(names(all_obs_fcc1) %in% names(tad_size))
      obs_tadSize_fcc1 <- tad_size[names(tad_size) %in% names(all_obs_fcc1)]
      stopifnot(length(obs_tadSize_fcc1) == length(all_obs_fcc1))
      obs_nFCC1 <- length(all_obs_fcc1)
      obs_nFCC1size3 <- sum(obs_tadSize_fcc1 == 3)
      obs_meanSizeFCC1 <- mean(obs_tadSize_fcc1)
      obs_medianSizeFCC1 <- median(obs_tadSize_fcc1)
      
      fcc_file <- file.path(pipFolder, hicds, exprds, "8cOnlyFCC_runAllDown", "prodSignedRatio_permDT.Rdata")
      if(!file.exists(fcc_file)) return(NULL)
      stopifnot(file.exists(fcc_file))
      fcc_perm_dt <- get(load(fcc_file))
      stopifnot(ncol(fcc_perm_dt) >= keepPermut)
      
      keepCols <- sample(x=1:ncol(fcc_perm_dt), size = keepPermut)
      stopifnot(length(keepCols) == keepPermut)
      fcc_perm_dt <- fcc_perm_dt[,keepCols ]
      all_perm <- apply(fcc_perm_dt[,1:5], 2, function(x) {
        fcc1_idx <- which(x == 1)
        perm_tad_fcc1 <- rownames(fcc_perm_dt)[fcc1_idx]
        stopifnot(perm_tad_fcc1 %in% names(tad_size))
        perm_tadSize_fcc1 <- tad_size[names(tad_size) %in% perm_tad_fcc1]
        stopifnot(length(perm_tad_fcc1) == length(fcc1_idx))
        t(data.frame(
          perm_nFCC1 = length(fcc1_idx),
        perm_nFCC1size3 = sum(perm_tadSize_fcc1 == 3),
        perm_meanSizeFCC1 = mean(perm_tadSize_fcc1),
        perm_medianSizeFCC1 = median(perm_tadSize_fcc1)))
      })
      permMean_nFCC1 <- mean(all_perm[1,])
      permMean_nFCC1size3 <- mean(all_perm[2,])
      permMean_meanSizeFCC1 <- mean(all_perm[3,])
      permMean_medianSizeFCC1 <- mean(all_perm[4,])
      
      data.frame(
        hicds = hicds,
        exprds = exprds,
        obs_nFCC1 = obs_nFCC1,
        obs_nFCC1size3 = obs_nFCC1size3,
        obs_meanSizeFCC1 = obs_meanSizeFCC1,
        obs_medianSizeFCC1 = obs_medianSizeFCC1,
        permMean_nFCC1 = permMean_nFCC1,
        permMean_nFCC1size3 = permMean_nFCC1size3,
        permMean_meanSizeFCC1 = permMean_meanSizeFCC1,
        permMean_medianSizeFCC1 = permMean_medianSizeFCC1,
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
  # load("FCC_RATIONBRTADS_HEATMAP_PERMG2T//all_dt.Rdata")
}

# load("CHECK_FCC1/all_dt.Rdata")
all_dt$obs_ratioSize3 <- all_dt$obs_nFCC1size3/all_dt$obs_nFCC1
all_dt$permMean_ratioSize3 <- all_dt$permMean_nFCC1size3/all_dt$permMean_nFCC1

m_all_dt <- melt(all_dt, id=c("hicds", "exprds"))
m_all_dt$varType <- gsub("(.+?)_.+", "\\1", m_all_dt$variable)
m_all_dt$varLab <- gsub("(.+?)_(.+)", "\\2", m_all_dt$variable)

m_all_dt$cmpType <- all_cmps[paste0(m_all_dt$exprds)]

nDS <- length(unique(file.path(m_all_dt$hicds, m_all_dt$exprds)))

plotTit <- "# TADs with FCC=1"
subTit <- paste0("all datasets (n=", nDS, "); mean permut (# perm. = ", keepPermut, ")")

sub1_dt <- m_all_dt[grepl("nFCC", m_all_dt$varLab),]

plot_myBox <- function(p){
  p <- p +   
    geom_boxplot(notch = TRUE, outlier.shape=NA)+
    geom_point(position=position_jitterdodge(),  alpha=0.5) +
    eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")"))) + 
    eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")"))) + 
    my_box_theme+
    theme(
      legend.text=element_text(size=12),
      axis.line=element_line(),
      axis.text.x = element_text(size=14)
    )
  return(p)
}
# ggplot(sub1_dt, aes(x=varLab, color=varType, y=value)) +
#   ggtitle(plotTit, subtitle = subTit)+
#   labs(fill ="", color="", x="", y="# TADs" )+
#   geom_boxplot(notch = TRUE, outlier.shape=NA) +
#   geom_point(aes(color=cmpType),position=position_jitterdodge(),  alpha=0.5) +
#   eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")"))) +
#   eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")"))) +
#   my_box_theme+
#   theme(
#     legend.text=element_text(size=12),
#     axis.line=element_line(),
#     axis.text.x = element_text(size=14)


p1 <- plot_myBox(ggplot(sub1_dt, aes(x=varLab, color=varType, y=value)) + 
  ggtitle(plotTit, subtitle = subTit)+
  labs(fill ="", color="", x="", y="# TADs" ))
  
outFile <- file.path(outFolder, paste0("nFCC1_nFCC1size3_obsPerm_boxplot.", plotType))
ggsave(p1, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


sub2_dt <- m_all_dt[!grepl("nFCC", m_all_dt$varLab) & grepl("ratio", m_all_dt$varLab),]

plotTit <- "Ratio of FCC=1 TADs of size # genes=3"

p2 <- plot_myBox(ggplot(sub2_dt, aes(x=varLab, color=varType, y=value)) + 
  ggtitle(plotTit, subtitle = subTit)+
    labs(fill ="", color="", x="", y="Ratio of FCC=1 TADs" ))

outFile <- file.path(outFolder, paste0("FCC1ratioSize3_obsPerm_boxplot.", plotType))
ggsave(p2, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


sub3_dt <- m_all_dt[!grepl("nFCC", m_all_dt$varLab) & ! grepl("ratio", m_all_dt$varLab),]

plotTit <- "Mean and median # genes of FCC=1 TADs"

p3 <- plot_myBox(ggplot(sub3_dt, aes(x=varLab, color=varType, y=value)) + 
  ggtitle(plotTit, subtitle = subTit)+
    labs(fill ="", color="", x="", y="Ratio of FCC=1 TADs" ))

outFile <- file.path(outFolder, paste0("FCC1meanMedianSize_obsPerm_boxplot.", plotType))
ggsave(p3, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



