hicds <- "ENCSR489OCU_NCI-H460_40kb"
exprds <- "TCGAluad_norm_luad"
tad_to_show <- "chr11_TAD390"

plotType <- "svg"

source("../settings.R")
require(ggsci)
require(ggrepel)
source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")

# Rscript geneRank_tadRank_selectTAD.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc chr11_TAD390
# Rscript geneRank_tadRank_selectTAD.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc chr10_TAD268
# Rscript geneRank_tadRank_selectTAD.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr11_TAD390
# Rscript geneRank_tadRank_selectTAD.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr10_TAD268
# 
# Rscript geneRank_tadRank_selectTAD.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr10_TAD16
# Rscript geneRank_tadRank_selectTAD.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr17_TAD162
# Rscript geneRank_tadRank_selectTAD.R ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker chr10_TAD16
# Rscript geneRank_tadRank_selectTAD.R ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker chr17_TAD162
# 

hicds="ENCSR489OCU_NCI-H460_40kb" 
exprds="TCGAluad_nonsmoker_smoker" 
tad_to_show="chr17_TAD162"

outFolder <- file.path("GENERANK_TADRANK_SELECTTAD")
dir.create(outFolder, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
hicds <- args[1]
exprds <- args[2]
tad_to_show <- args[3]

stopifnot(hicds %in% names(hicds_names))
stopifnot(exprds %in% names(exprds_names))

hicds_lab <- hicds_names[paste0(hicds)]
exprds_lab <- exprds_names[paste0(exprds)]

inDT <- get(load(file.path(runFolder, "GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata")))

ds_dt <- inDT[inDT$hicds == hicds & inDT$exprds == exprds,]

tad_dt <- ds_dt[,c("hicds", "exprds", "region", "tad_adjCombPval")]
tad_dt <- unique(tad_dt)
tad_dt$tad_rank <- rank(tad_dt$tad_adjCombPval, ties="min")
tad_dt$rev_tad_rank <- rank(-tad_dt$tad_adjCombPval, ties="min")

ds_dt <- merge(ds_dt, tad_dt[, c("hicds", "exprds", "region", "rev_tad_rank")], all=TRUE, by=c("hicds", "exprds", "region"))

nTADs <- length(unique(ds_dt$region))

head(ds_dt[order(ds_dt$tad_adjCombPval),])


geneSignifThresh <- 0.01
tadSignifThresh <- 0.01

geneRankSignifLimit <- max(ds_dt$gene_rank[ds_dt$adj.P.Val <=geneSignifThresh])
tadRankSignifLimit <- max(ds_dt$tad_rank[ds_dt$tad_adjCombPval <=geneSignifThresh])

geneRankSignifLimit_resc <- geneRankSignifLimit/max(ds_dt$gene_rank)
tadRankSignifLimit_resc <- tadRankSignifLimit/max(ds_dt$tad_rank)


stopifnot(!is.na(geneRankSignifLimit))
stopifnot(!is.na(tadRankSignifLimit))

tad_to_show_col <-  "red"
hide_col <- "grey"
showCex <- 1.2

nTop <- 10

top_reg <- unique(ds_dt$region[ds_dt$tad_rank <= nTop])
last_reg <- unique(ds_dt$region[ds_dt$rev_tad_rank <= nTop])

ds_dt$showCol <- ifelse(ds_dt$region == tad_to_show, tad_to_show_col, hide_col)

stopifnot(!is.na(ds_dt$showCol))
plot_dt <- ds_dt[ds_dt$region %in% top_reg | ds_dt$region %in% last_reg ,]

plotTit <- paste0(hicds_lab, "\n" , exprds_lab)

#########################################################################################################################
#### BARPLOTS - all gene ranks
#########################################################################################################################
maxGeneRank <- max(ds_dt$gene_rank)
maxTADrank <- max(ds_dt$tad_rank)

ds_dt$gene_rank_rel <- ds_dt$gene_rank/maxGeneRank
ds_dt$tad_rank_rel <- ds_dt$tad_rank/maxTADrank


geneBar_pos <- 1
tadBar_pos <- 2
axisOffset <- 0.5
labOffset <- 0.02

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", tad_to_show, "_geneRank_tadRank_showBars_segments.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))

par(bty="n", family=fontFamily)
plot(NULL,
     xlab="",
     ylab="",
     axes=F,
     # main = paste0(hicds, " - " , exprds),
     main = plotTit,
     # sub=paste0(),
     cex.axis=plotCex,
     cex.lab = plotCex,
     cex.main = plotCex,
     xlim = c(geneBar_pos-axisOffset, tadBar_pos+axisOffset),
     # ylim = c(-axisOffset, 1+axisOffset)
     ylim = c(1+axisOffset, -axisOffset )
)
# mtext(side=3, paste0(hicds, " - ", exprds), cex=plotCex)
mtext(side=3, paste0(tad_to_show), cex=plotCex-0.2, line=-1)

segments(x0=geneBar_pos,
         x1=tadBar_pos,
         y0=ds_dt$gene_rank_rel[ds_dt$showCol == tad_to_show_col],
         y1=ds_dt$tad_rank_rel[ds_dt$showCol == tad_to_show_col],
         col = tad_to_show_col
)
# add bar for the genes
segments(x0=geneBar_pos, y0=0, x1=geneBar_pos, y1=1, lwd=5)
# add bar for the TADs
segments(x0=tadBar_pos, y0=0, x1=tadBar_pos, y1=1, lwd=5)

text(
  x = c(geneBar_pos, tadBar_pos),
  # y = 1 + 0.2,
  y = 0- 0.2,
  labels=c("Gene ranks", "TAD ranks"),
  cex = plotCex
)
show_dt <- ds_dt[ds_dt$region == tad_to_show,]

text(
  y = c(unique(show_dt$tad_rank_rel), unique(show_dt$gene_rank_rel)),
  x = c(rep(tadBar_pos, length(unique(show_dt$tad_rank_rel)))+labOffset, rep(geneBar_pos, length(unique(show_dt$gene_rank_rel)))-labOffset),
  labels=c(unique(show_dt$tad_rank), unique(show_dt$gene_rank)),
  cex = 1,
  pos=c(rep(4, length(unique(show_dt$tad_rank_rel))), rep(2, length(unique(show_dt$gene_rank_rel))))
)
segments(
  y0 = c(unique(show_dt$tad_rank_rel), unique(show_dt$gene_rank_rel)),
  y1 = c(unique(show_dt$tad_rank_rel), unique(show_dt$gene_rank_rel)),
  x0 = c(rep(tadBar_pos, length(unique(show_dt$tad_rank_rel))), rep(geneBar_pos, length(unique(show_dt$gene_rank_rel)))),
  x1 = c(rep(tadBar_pos, length(unique(show_dt$tad_rank_rel)))+labOffset, rep(geneBar_pos, length(unique(show_dt$gene_rank_rel)))-labOffset),
  lwd=1
)

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

##################################################################### GGREPEL VERSION
tad_y <- 0
gene_y <- 0.5

axisOffset <- 0.5

show_dt <- ds_dt[ds_dt$region == tad_to_show,]
show_dt$TAD_y <- tad_y
show_dt$gene_y <- gene_y

show_dt$tadSignifLab <- ifelse(show_dt$tad_adjCombPval <= 0.001, "***", 
                                      ifelse(show_dt$tad_adjCombPval <= 0.01, "**",
                                             ifelse(show_dt$tad_adjCombPval <= 0.05, "*", "")))
show_dt$geneSignifLab <- ifelse(show_dt$adj.P.Val <= 0.001, "***", 
                               ifelse(show_dt$adj.P.Val <= 0.01, "**",
                                      ifelse(show_dt$adj.P.Val <= 0.05, "*", "")))


show_dt_m1 <- show_dt[,c("entrezID", "tad_rank_rel", "tad_rank", "TAD_y", "tadSignifLab")]
colnames(show_dt_m1)[2:5] <- c("rank_x_rel", "rank_x", "y_pos", "signif_lab") 
show_dt_m1$signif_rank_lab <- paste0(show_dt_m1$rank_x, "\n", show_dt_m1$signif_lab)

show_dt_m2 <- show_dt[,c("entrezID", "gene_rank_rel", "gene_rank", "gene_y", "geneSignifLab")]
colnames(show_dt_m2)[2:5] <- c("rank_x_rel","rank_x", "y_pos", "signif_lab") 
show_dt_m2$signif_rank_lab <- paste0(show_dt_m2$signif_lab, "\n", show_dt_m2$rank_x)

show_dt_m <- rbind(show_dt_m1,show_dt_m2)

myTit <- paste0(tad_to_show)
myTit <- paste0(hicds_lab, " - ", exprds_lab, ": ", tad_to_show)
subTit <- paste0(hicds_lab, " - ", exprds_lab)
subTit <- ""

geneLab <- "Gene ranks"
tadLab <- "TAD ranks"

repel_offset <- 0.4
repel_offset_axis <- 0.2

rankLabSize <- 5
scaleLabSize <- 3

relativeMin <- 0
relativeMax <- 1

show_dt_m$y_nudge <- ifelse(show_dt_m$y_pos == tad_y, -repel_offset, repel_offset) 

save(show_dt_m, file="show_dt_m.Rdata", version=2)

signif_bar_limit_offset <- 0.1
signif_bar_col <- "forestgreen"
feature_dot_col <- "red"


  

ranks_p <- ggplot(show_dt_m, aes(x=rank_x_rel, y =y_pos, label=show_dt_m$rank_x))+
  ggtitle(myTit, subtitle=subTit)+
  scale_y_continuous(name="", breaks=c(tad_y, gene_y), labels=c(tadLab, geneLab), limits=c(tad_y-axisOffset,gene_y+axisOffset))+
  geom_point(color=feature_dot_col)+
  labs(y="", x="")+
  theme_void() +
  theme(
    text = element_text(family=fontFamily),
    axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
    axis.text.y = element_text(size=14, hjust=0.5, vjust=0.5, face="bold"),
    axis.text.x = element_blank(),
    plot.title = element_text(hjust=-0.1, vjust=0, size = 16, face="bold"),
    plot.subtitle = element_text(hjust=-0.1, vjust=0, size = 14, face="italic"),
    legend.title = element_text(face="bold")
  )+
  geom_segment(data=show_dt, aes(x=tad_rank_rel, y=TAD_y, xend=gene_rank_rel, yend=gene_y), inherit.aes=F) +
  geom_segment(aes(x=0, xend=1, y=tad_y, yend = tad_y), colour="grey", inherit.aes = F, linetype=3) +
  geom_segment(aes(x=0, xend=1, y=gene_y, yend = gene_y), colour="grey", inherit.aes = F, linetype=3) 

ranks_p_bar <- ranks_p + 
  geom_segment(x=geneRankSignifLimit_resc, 
               xend=geneRankSignifLimit_resc, y=gene_y-signif_bar_limit_offset, 
               yend=gene_y+signif_bar_limit_offset, lineend="square", color=signif_bar_col, alpha=0.2, size=1)+
  geom_segment(x=tadRankSignifLimit_resc, 
               xend=tadRankSignifLimit_resc, y=tad_y-signif_bar_limit_offset, 
               yend=tad_y+signif_bar_limit_offset, lineend="square", color=signif_bar_col, alpha=0.2, size=1)
  
  
  
lab_dt_1 <- data.frame(rank_x_rel = c(relativeMin, relativeMax, relativeMin, relativeMax), 
                     rank_x=c(1, length(unique(ds_dt$region)), 1, length(unique(ds_dt$entrezID))),
                     y_pos = c(tad_y, tad_y, gene_y, gene_y),
                     signif_lab = "",
                     signif_rank_lab = "")
lab_dt_1$y_nudge <- ifelse(lab_dt_1$y_pos == tad_y, -repel_offset_axis, repel_offset_axis) 
lab_dt_2 <- show_dt_m[,c("rank_x_rel","rank_x", "y_pos", "y_nudge", "signif_lab", "signif_rank_lab")]
lab_dt_1$labSize <- scaleLabSize
lab_dt_1$labCol <- "darkgrey"
lab_dt_2$labSize <- rankLabSize
lab_dt_2$labCol <- "black"

save(ds_dt, file="ds_dt.Rdata", version=2)
save(lab_dt_1, file="lab_dt_1.Rdata", version=2)
save(lab_dt_2, file="lab_dt_2.Rdata", version=2)


lab_dt_all <- rbind(lab_dt_1, lab_dt_2)
lab_dt_all <- unique(lab_dt_all)

lab_dt_all$labCol2 <- lab_dt_all$labCol
lab_dt_all$labCol2[lab_dt_all$signif_lab != "" & lab_dt_all$signif_rank_lab != ""] <- "forestgreen"

save(lab_dt_all, file="lab_dt_all.Rdata", version=2)

all_p <- ranks_p +   
                  geom_text_repel(data=lab_dt_all, aes(x=rank_x_rel, y=y_pos, label=signif_rank_lab),inherit.aes = F,
                  nudge_y      = lab_dt_all$y_nudge,
                  direction    = "x",
                  size = lab_dt_all$labSize,
                  angle        = 0,
                  vjust        = 0,
                  hjust=0.5,
                  segment.size = 0.3, 
                  colour=lab_dt_all$labCol
  ) 

  
outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", tad_to_show, "_geneRank_tadRank_showBars_segments_GG_withSignifLabs.", plotType))
ggsave(plot = all_p, filename = outFile, height=myHeightGG/2.5, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))

all_p <- ranks_p +   
  geom_text_repel(data=lab_dt_all, aes(x=rank_x_rel, y=y_pos, label=rank_x),inherit.aes = F,
                  nudge_y      = lab_dt_all$y_nudge,
                  direction    = "x",
                  size = lab_dt_all$labSize,
                  angle        = 0,
                  vjust        = 0,
                  hjust=0.5,
                  segment.size = 0.3, 
                  colour=lab_dt_all$labCol2
  ) 


outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", tad_to_show, "_geneRank_tadRank_showBars_segments_GG_withSignifCols.", plotType))
ggsave(plot = all_p, filename = outFile, height=myHeightGG/2.5, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))


all_p <- ranks_p_bar +   
  geom_text_repel(data=lab_dt_all, aes(x=rank_x_rel, y=y_pos, label=rank_x),inherit.aes = F,
                  nudge_y      = lab_dt_all$y_nudge,
                  direction    = "x",
                  size = lab_dt_all$labSize,
                  angle        = 0,
                  vjust        = 0,
                  hjust=0.5,
                  segment.size = 0.3, 
                  colour=lab_dt_all$labCol
  ) 

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", tad_to_show, "_geneRank_tadRank_showBars_segments_GG_withSignifBars.", plotType))
ggsave(plot = all_p, filename = outFile, height=myHeightGG/2.5, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))




stop("-ok\n")

# 
# lab_dt <- data.frame(x = c(tad_y,gene_y,tad_y,gene_y),
#                      y = c(tad_y, tad_y, gene_y, gene_y),
#                      value=c(1, length(unique(ds_dt$region)), 1, length(unique(ds_dt$entrezID))))

# all_p <- ranks_p2 +   
#   geom_text_repel(data=lab_dt, aes(x=x, y=y, label=value),inherit.aes = F,
#   nudge_y      = c(-repel_offset_axis ,-repel_offset_axis, repel_offset_axis,repel_offset_axis),
#   direction    = "x",
#   size = scaleLabSize,
#   angle        = 0,
#   vjust        = 0,
#   hjust=0.5,
#   segment.size = 0.3, 
#   colour="darkgrey"
# ) 

# geom_text_repel(
#   # nudge_y      = c(-repel_offset ,-repel_offset, -repel_offset, repel_offset,repel_offset,repel_offset),
#   nudge_y      =  show_dt_m$y_nudge,
#   size = rankLabSize,
#   direction    = "x",
#   angle        = 0,
#   vjust        = 0,
#   hjust=0.5,
#   segment.size = 0.3
# ) +


