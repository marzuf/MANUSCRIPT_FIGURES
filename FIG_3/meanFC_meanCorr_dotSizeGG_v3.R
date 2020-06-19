

options(scipen=100)


# Rscript meanFC_meanCorr_dotSizeGG_v3.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc chr11_TAD390 chr10_TAD268
# Rscript meanFC_meanCorr_dotSizeGG_v3.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr11_TAD390 chr10_TAD268
# Rscript meanFC_meanCorr_dotSizeGG_v3.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr10_TAD16 chr17_TAD162
# Rscript meanFC_meanCorr_dotSizeGG_v3.R ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker chr10_TAD16 chr17_TAD162
# 

script_name <- "meanFC_meanCorr_dotSizeGG_v3.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(ggplot2)
require(reshape2)
require(ggpubr)
require(ggrepel)


registerDoMC(4)

source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")



outFolder <- "MEANFC_MEANCORR_DOTSIZEGG_V3"
dir.create(outFolder, recursive = TRUE)

all_cols[all_cols == "red"] <- "brown3"
all_cols[all_cols == "blue"] <- "darkblue"
all_cols[all_cols == "green"] <- "forestgreen"

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

source("../settings.R")


myWidthGG <- 9
myHeightGG <- 7

setDir <- "/media/electron"
setDir <- ""

mainFolder <- file.path(runFolder)
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(pipFolder)
stopifnot(dir.exists(pipFolder))


ex_hicds <- "ENCSR489OCU_NCI-H460_40kb"
ex_exprds <- "TCGAlusc_norm_lusc"
hicds_tit <- "ENCSR489OCU_NCI-H460_40kb"
exprds_tit <- "TCGAlusc_norm_lusc"
tads_to_annot <- c("chr11_TAD390", "chr10_TAD268")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2)
ex_hicds <- args[1]
ex_exprds <- args[2]
hicds_tit <- ex_hicds
exprds_tit <- ex_exprds
if(length(args) > 2) {
  tads_to_annot <- args[3:length(args)]  
}else{
  tads_to_annot <- NULL
}

stopifnot(ex_hicds %in% names(hicds_names))
hicds_tit <- hicds_names[paste0(ex_hicds)]

stopifnot(ex_exprds %in% names(exprds_names))
exprds_tit <- exprds_names[paste0(ex_exprds)]


plotCex <- 1.2

signifThresh <- 0.01

final_DT <- get(load(file.path(runFolder, "CREATE_FINAL_TABLE/all_result_dt.Rdata")))
final_DT_signif <- final_DT[final_DT$adjPvalComb <= signifThresh,]

cat(paste0(ex_hicds, "\n", ex_exprds, "\n"))

ex_DT <- final_DT[final_DT$hicds == ex_hicds &
                    final_DT$exprds == ex_exprds,
                  ]

stopifnot(nrow(ex_DT) > 0)

ex_DT <- ex_DT[order(ex_DT$adjPvalComb),]
ex_DT_signif <- ex_DT[ex_DT$adjPvalComb <= signifThresh,]
stopifnot(nrow(ex_DT_signif) > 0)

# if <= 0.01 => log10 <= 2 => cex = 0.7
ex_DT$adjPvalComb_log10 <- -log10(ex_DT$adjPvalComb)

outFile <- file.path(outFolder, paste0(ex_hicds, "_", ex_exprds, "_dt.Rdata"))
save(ex_DT, file=outFile, version=2)
save(ex_DT, file="ex_DT.Rdata", version=2)

x_lab <- "TAD mean LogFC"
y_lab <- "TAD mean intraCorr"

plotTit <- paste0(hicds_tit, " - ", exprds_tit)
subTit <- paste0("# TADs = ", nrow(ex_DT), "; # signif. TADs = ", nrow(ex_DT_signif))

ex_DT$adjPvalComb_log10_col <- ex_DT$adjPvalComb_log10
ex_DT$adjPvalComb_log10_col[ex_DT$adjPvalComb_log10_col <= -log10(0.05)] <- -log10(0.05)

ex_DT$adjPvalComb_crop <- ex_DT$adjPvalComb 
ex_DT$adjPvalComb_crop[ex_DT$adjPvalComb_crop >= 0.05] <- 0.05

ex_DT$adjPvalCombSignif <- ex_DT$adjPvalComb <= signifThresh


ex_DT$dotcol_labs1 <- ifelse(ex_DT$meanLogFC > 0, "logFC>0", "logFC<=0")
ex_DT$dotcol_labs2 <- ifelse(ex_DT$adjPvalComb <= 0.01, "adjP<=0.01",
                             ifelse(ex_DT$adjPvalComb <= 0.05, "adjP<=0.05", NA))

ex_DT$dotcol_lab <- interaction(ex_DT$dotcol_labs1, ex_DT$dotcol_labs2)
ex_DT$dotcol_lab <- as.character(ex_DT$dotcol_lab)
stopifnot(which(is.na(ex_DT$dotcol_labs2)) == which(ex_DT$adjPvalComb > 0.05))
ex_DT$dotcol_lab[is.na(ex_DT$dotcol_labs2)] <- "adjP>0.05"
stopifnot(!is.na(ex_DT$dotcol_lab))




# library(yarrr)
# basel_pal_transp <- piratepal(palette = "basel",  trans = .8)
# save(basel_pal_transp, file=file.path(outFolder, "basel_pal_transp.Rdata"), version=2)
basel_pal_transp <- get(load(file.path(gsub("V3", "V2", outFolder), "basel_pal_transp.Rdata")))
# basel_pal <- piratepal(palette = "basel")
# save(basel_pal, file=file.path(outFolder, "basel_pal.Rdata"), version=2)
basel_pal <- get(load(file.path(gsub("V3", "V2", outFolder), "basel_pal.Rdata")))
# xmen_pal_transp <- piratepal(palette = "xmen",  trans = .8)
# save(xmen_pal_transp, file=file.path(outFolder, "xmen_pal_transp.Rdata"), version=2)
xmen_pal_transp <- get(load(file.path(gsub("V3", "V2", outFolder), "xmen_pal_transp.Rdata")))

strongRed <-  as.character(basel_pal[2])
lightRed <- as.character(basel_pal_transp[2])

strongBlue <- as.character(basel_pal[1])
lightBlue <- as.character(basel_pal_transp[1])

mycols <- setNames(c(strongRed, lightRed, strongBlue, lightBlue, "darkgrey"), 
                   c("logFC>0.adjP<=0.01","logFC>0.adjP<=0.05",
                   "logFC<=0.adjP<=0.01","logFC<=0.adjP<=0.05", 
                   "adjP>0.05"))

mycols_labs <- setNames(c("mean LogFC > 0 &\nadj. p-val <=0.01",
                          "mean LogFC > 0 &\nadj. p-val <=0.05",
                          "mean LogFC <= 0 &\nadj. p-val <=0.01",
                          "mean LogFC <= 0 &\nadj. p-val <=0.05", 
                          "adj. p-val > 0.05"),
                   c("logFC>0.adjP<=0.01","logFC>0.adjP<=0.05",
                     "logFC<=0.adjP<=0.01","logFC<=0.adjP<=0.05", 
                     "adjP>0.05"))



# stopifnot(!is.na(ex_DT))

fc_corr_p <- ggplot(ex_DT, aes(x=meanLogFC, y = meanCorr, fill = dotcol_lab, color=dotcol_lab )) +
  ggtitle(plotTit, subtitle = subTit)+
  
  scale_fill_manual(values=mycols, labels = mycols_labs)+  
                    # labels=setNames(names(mycols), mycols)
  scale_color_manual(values=mycols,
                    labels=setNames(names(mycols), mycols), guide=F)+
  
  # scale_color_manual(values=ratioDown_cols,
  #                    labels=setNames(names(ratioDown_cols),ratioDown_cols)) +
  guides(color=F)+
  
  geom_point(data=ex_DT[ex_DT$dotcol_lab=="adjP>0.05",],
             shape=21,
             aes(size = adjPvalComb_log10_col, color = dotcol_lab), stroke=1) +
  geom_point(data=ex_DT[ex_DT$dotcol_lab!="adjP>0.05",],
             shape=21,
             aes(size = adjPvalComb_log10_col, color = dotcol_lab), stroke=1)+
  
  labs( size= "TAD adj. p-val", fill="")+#, color="TAD ratio\ndown-reg. genes")+
#  labs( fill = "TAD adj. p-val",size="TAD adj. p-val")+#, color="TAD ratio\ndown-reg. genes")+
  geom_label_repel(data = ex_DT[ex_DT$region %in% tads_to_annot,],
    aes(x= meanLogFC, 
        y=meanCorr, label=region),
    inherit.aes = F,
    min.segment.length = unit(0, 'lines'),
    force = 10,
    segment.size = 0.5,
    nudge_y = 0.1
  ) + 
  scale_x_continuous(name=paste0(x_lab),
                     breaks = scales::pretty_breaks(n = 10))+
  scale_y_continuous(name=paste0(y_lab),
                     breaks = scales::pretty_breaks(n = 10))+
  scale_size(range=c(0.5,8),#expand=c(2,0),
             breaks=c(-log10(0.05), -log10(0.01), -log10(0.005), -log10(0.001)),
             labels=c(">= 0.05","0.01"," 0.005","0.001"),
             guide="legend") +
  my_box_theme + 
  
  guides(size = guide_legend(order=1))+
  guides(fill = guide_legend(override.aes = list(color="white", size=4),order=2))+
  # guides(color = guide_legend(order=3))+
  
  theme(axis.line = element_line(),
        axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5, color="black"),
        axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5, color="black"),
        panel.grid = element_blank(),
        legend.key = element_rect(fill = NA),
        # legend.text =  element_text(size=10),
        legend.title =  element_text(size=12,face="bold"),
        panel.grid.major.y =  element_blank(),
        panel.grid.minor.y =  element_blank(),
        legend.text = element_text(margin = margin(b = 0.1, unit = 'in'), size=10)
        ) +
  geom_hline(yintercept = 0, linetype=2, color="darkgrey")+
  geom_vline(xintercept = 0, linetype=2, color="darkgrey")


outFile <- file.path(outFolder, paste0(ex_hicds, "_", ex_exprds, "_meanIntraCorr_meanLogFC_dotplot_with_signifGG.", plotType))
ggsave(fc_corr_p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))





#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
cat("***** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


