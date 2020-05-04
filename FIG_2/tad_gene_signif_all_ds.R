

# Rscript tad_vs_gene_signif_all_ds.R

require(doMC)
require(foreach)
require(ggplot2)
require(ggsci)
require(reshape2)
require(ggpubr)

registerDoMC(40)

plotType <- "svg"
source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")

# myWidth <- myWidth * 1.2

outFolder <- "TAD_VS_GENE_SIGNIF_ALL_DS"
dir.create(outFolder, recursive = TRUE)

buildData <- TRUE

ggsci_pal <- "d3"
ggsci_subpal <- ""

plotMargin <- c(1,2,1,1)


myHeightGG <- 5
myWidthGG <- 7


rankTieMeth <- "min"

genes_nTop <- 100
topTADs <- 10
lastTADs <- 10

tads_nTop <- 10
topGenes <- 100
lastGenes <- 100

yOffset <- 0.2 # for the balloonplot

inDT <- get(load(file.path(runFolder, "GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata")))


##############################################################
############################### PREP DATA FOR ALL CATEGORIES
##############################################################


inDT2 <- inDT
inDT2$dataset <- file.path(inDT2$hicds, inDT2$exprds)
inDT_dt <- do.call(rbind, by(inDT2, inDT2$dataset, function(x) {
  x$gene_rank_std <- x$gene_rank/max(x$gene_rank)
  x$gene_rank2 <- rank(x$adj.P.Val, ties=rankTieMeth)
  x$gene_rank_rev <- rank(-x$adj.P.Val, ties=rankTieMeth)

  
  tad_dt <- x[,c("region", "tad_adjCombPval", "tad_rank")]
  tad_dt <- unique(tad_dt)
  tad_dt$tad_rank_std <- tad_dt$tad_rank/max(tad_dt$tad_rank)
  tad_dt$tad_rank2 <- rank(tad_dt$tad_adjCombPval, ties=rankTieMeth)
  tad_dt$tad_rank_rev <- rank(-tad_dt$tad_adjCombPval, ties=rankTieMeth)
  tad_dt$tad_rank <- NULL
  tad_dt$tad_adjCombPval <- NULL
  merge(x, tad_dt, by=c("region"))
  
  
  
}))
stopifnot(inDT_dt$gene_rank==inDT_dt$gene_rank2)
stopifnot(! (inDT_dt$gene_rank_rev <= lastGenes & inDT_dt$gene_rank <= topGenes))
stopifnot(inDT_dt$tad_rank==inDT_dt$tad_rank2)
stopifnot(! (inDT_dt$tad_rank_rev <= lastTADs & inDT_dt$tad_rank <= topTADs))

inDT_dt$gene_rank_cat <- ifelse(inDT_dt$gene_rank_rev <= lastGenes, paste0("last ", lastGenes, " genes"), 
                                ifelse(inDT_dt$gene_rank <= topGenes, paste0("top ", topGenes, " genes"), "other")) 

inDT_dt$tad_rank_cat <- ifelse(inDT_dt$tad_rank_rev <= lastTADs, paste0("last ", lastTADs, " TADs"), 
                               ifelse(inDT_dt$tad_rank <= topTADs, paste0("top ", topTADs, " TADs"), "other")) 

stopifnot(!is.na(inDT_dt))

gene_levels <- c(paste0("last ", lastGenes, " genes"), "other",  paste0("top ", topGenes, " genes"))

tad_levels <- c(paste0("last ", lastTADs, " TADs"), "other", paste0("top ", topTADs, " TADs"))


##############################################################
############################### PREP DATA FOR ALL CATEGORIES
##############################################################

annotateSize <- 6
annotateSizeL <- 8

nDS <- length(unique(inDT_dt$dataset))



myTit <- "Genes signif. gene-level vs. TAD-level"

agg_dt <- aggregate(entrezID~gene_rank_cat + tad_rank_cat, data = inDT_dt, FUN=length)
colnames(agg_dt)[colnames(agg_dt) == "entrezID"] <- "nGenes"

agg_dt$gene_rank_cat <- factor(agg_dt$gene_rank_cat, levels=gene_levels)
agg_dt$tad_rank_cat <- factor(agg_dt$tad_rank_cat, levels=tad_levels)
stopifnot(!is.na(agg_dt))

agg_dt <- agg_dt[order(as.numeric(agg_dt$gene_rank_cat), as.numeric(agg_dt$tad_rank_cat)),]

balloon_p1 <- ggballoonplot(agg_dt, x="gene_rank_cat", y ="tad_rank_cat", fill = "nGenes")+
  ggtitle(myTit, subtitle = paste0("# datasets = ", nDS))+
  scale_fill_viridis_c(option = "C", trans="log10") +
  scale_size( trans="log10") +
  labs(size = "# genes", fill = "# genes")+
  theme(
    plot.title = element_text(size = 18, face ="bold", hjust=0.5),
    plot.subtitle = element_text(size = 16, face ="italic", hjust=0.5),
    axis.text =element_text(size=14)
  )


outFile <- file.path(outFolder, 
                     paste0("nSignifGenes_topGenes",genes_nTop, "_vs_topTADs", tads_nTop, "_balloonplot.", plotType))
ggsave(plot = balloon_p1, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

balloon_p1b <- balloon_p1 + annotate("text",
                      x=as.numeric(agg_dt$gene_rank_cat),
                      y=as.numeric(agg_dt$tad_rank_cat)+yOffset,
                      label = agg_dt$nGenes,
                      size = annotateSize
                      )

outFile <- file.path(outFolder, 
                     paste0("nSignifGenes_topGenes",genes_nTop, "_vs_topTADs", tads_nTop, "_balloonplotWithText.", plotType))
ggsave(plot = balloon_p1b, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, 
                     paste0("nSignifGenes_topGenes",genes_nTop, "_vs_topTADs", tads_nTop, "_table.txt"))
write.table(agg_dt, file=outFile, col.names = T, row.names = F, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))


agg_dt2 <- agg_dt[!(agg_dt$gene_rank_cat=="other" & agg_dt$tad_rank_cat=="other"),] 
balloon_p2 <- ggballoonplot(agg_dt2, x="gene_rank_cat", y ="tad_rank_cat", 
                            fill = "nGenes", size="nGenes", 
                            font.label=c(14, "bold", "red"), 
                            show.label=FALSE)+
  ggtitle(myTit, subtitle = paste0("# datasets = ", nDS))+
  scale_size( trans="log10")+
  scale_fill_viridis_c(option = "C", trans="log10")+
  # scale_fill_viridis_c(option = 'C', trans = "log10",  breaks = trans_breaks("log10", function(x) 10^x)) +
  labs(size = "# genes", fill = "# genes")+
  theme(
    plot.title = element_text(size = 18, face ="bold", hjust=0.5),
    plot.subtitle = element_text(size = 16, face ="italic", hjust=0.5),
    axis.text =element_text(size=14)
  )

outFile <- file.path(outFolder, 
                     paste0("nSignifGenes_topGenes",genes_nTop, "_vs_topTADs", tads_nTop, "_balloonplot_noOther.", plotType))
ggsave(plot = balloon_p2, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

balloon_p2b <- balloon_p2 + annotate("text",
                                     x=as.numeric(agg_dt2$gene_rank_cat),
                                     y=as.numeric(agg_dt2$tad_rank_cat)+yOffset,
                                     label = agg_dt2$nGenes,
                                     size = annotateSize)+
              annotate("text",
                       x=as.numeric(agg_dt$gene_rank_cat[(agg_dt$gene_rank_cat=="other" & agg_dt$tad_rank_cat=="other")] ),
                       y=as.numeric(agg_dt$tad_rank_cat[(agg_dt$gene_rank_cat=="other" & agg_dt$tad_rank_cat=="other")] ),
                       label = agg_dt$nGenes[(agg_dt$gene_rank_cat=="other" & agg_dt$tad_rank_cat=="other")],
                       size = annotateSizeL)
              

outFile <- file.path(outFolder, 
                     paste0("nSignifGenes_topGenes",genes_nTop, "_vs_topTADs", tads_nTop, "_balloonplot_noOtherWithText.", plotType))
ggsave(plot = balloon_p2b, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


rank_hex_p <- ggplot(inDT_dt, aes(x = gene_rank, y = tad_rank)) +
  ggtitle("TAD rank vs. gene rank") +
  scale_fill_viridis_c(option = "C")+
  stat_binhex(show.legend = T, bins = 20)+
  labs(fill = "", x ="Gene rank", y= "TAD rank") +
    my_box_theme +
  theme(  panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
          panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1))

outFile <- file.path(outFolder, 
                     paste0("TADrank_vs_geneRank_density_hexplot.", plotType))
ggsave(plot = rank_hex_p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



##############################################################
############################### stacked barplot top10 TADs
##############################################################


plot_dt <- inDT_dt[inDT_dt$tad_rank <= tads_nTop,]

  
stopifnot(!duplicated(plot_dt))
agg_dt <- aggregate(entrezID ~ tad_rank+gene_rank_cat, data = plot_dt, FUN=length)
colnames(agg_dt)[colnames(agg_dt) == "entrezID"] <- "nGenes"

geneCat_tadRank_p <- ggplot(agg_dt, aes(x=tad_rank, y = nGenes, color = gene_rank_cat, fill = gene_rank_cat)) +
  ggtitle("Gene-level vs. TAD-level ranks", subtitle="by TAD rank")+
  geom_bar(stat="identity")+
labs(fill ="Gene rank", x=paste0("TAD rank (top ", tads_nTop, ")"), y="# of genes")+
  scale_x_continuous(breaks=c(1:max(agg_dt$tad_rank)), labels = c(1:max(agg_dt$tad_rank)))+
  guides(color=FALSE)+
    eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")")))+
  eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")")))+
  my_box_theme +
  theme(
    # axis.line = element_line()
    axis.line.x = element_line()
        )
outFile <- file.path(outFolder, 
                     paste0("geneCat_by_TADrank_barplot.", plotType))
ggsave(plot = geneCat_tadRank_p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



##############################################################
############################### geom area  top10 TADs
##############################################################


plot_dt <- inDT_dt[inDT_dt$gene_rank <= genes_nTop,]



stopifnot(!duplicated(plot_dt))
agg_dt <- aggregate(entrezID ~ gene_rank+tad_rank_cat, data = plot_dt, FUN=length)
colnames(agg_dt)[colnames(agg_dt) == "entrezID"] <- "nGenes"

ggplot(agg_dt, aes(x=gene_rank, y = nGenes, color = tad_rank_cat, fill = tad_rank_cat)) +
  geom_bar(stat="identity")


tadCat_geneRank_p <- ggplot(agg_dt, aes(x=gene_rank, y = nGenes, color = tad_rank_cat, fill = tad_rank_cat)) +
  ggtitle("TAD-level vs. gene-level ranks", subtitle="by gene rank")+
  geom_bar(stat="identity")+
  labs(fill ="TAD rank", x=paste0("Gene rank (top ", genes_nTop, ")"), y="# of genes")+
  scale_x_continuous(breaks=c(1:max(agg_dt$gene_rank)), labels = c(1:max(agg_dt$gene_rank)))+
  guides(color=FALSE)+
  eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")")))+
  eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")")))+
  my_box_theme +
  theme(
    # axis.line = element_line()
    axis.line.x = element_line(),
    axis.text.x = element_blank(),
    axis.ticks =  element_blank()
  )

outFile <- file.path(outFolder, 
                     paste0("TADcat_by_geneRank_barplot.", plotType))
ggsave(plot = tadCat_geneRank_p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


##############################################################
############################### densplot top10 TADs
##############################################################



my_tit <- paste0("Gene rank vs. TAD rank")
my_sub <- paste0("(TAD rank <= ", tads_nTop, ")")

outFile <- file.path(outFolder, paste0("geneRank_tadRank_topTADs_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="l", family=fontFamily)
boxplot(
  plot_dt$gene_rank ~plot_dt$tad_rank,
  xlab = "TAD rank",
  ylab = "Gene rank",
  main = my_tit,
  cex.axis = axisCex,
  cex.main = mainCex,
  cex.lab = labCex 
)
mtext(side = 3, text = paste0(my_sub), cex = subCex)
foo <- dev.off()  
cat(paste0("... written: ", outFile, "\n"))



##############################################################
############################### densplot top10 TADs
##############################################################

tads_nTop <- 10

plot_dt <- inDT
plot_dt <- plot_dt[plot_dt$tad_rank <= tads_nTop,]

my_tit <- paste0("Gene rank vs. TAD rank")
my_sub <- paste0("(TAD rank <= ", tads_nTop, ")")

outFile <- file.path(outFolder, paste0("geneRank_tadRank_topTADs_boxplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="l", family=fontFamily)
boxplot(
  plot_dt$gene_rank ~plot_dt$tad_rank,
  xlab = "TAD rank",
  ylab = "Gene rank",
  main = my_tit,
  cex.axis = axisCex,
  cex.main = mainCex,
  cex.lab = labCex 
)
mtext(side = 3, text = paste0(my_sub), cex = subCex)
foo <- dev.off()  
cat(paste0("... written: ", outFile, "\n"))


##############################################################
############################### densplot top100 genes
##############################################################

genes_nTop <- 100

plot_dt <- inDT
plot_dt <- plot_dt[plot_dt$gene_rank <= genes_nTop,]

my_tit <- paste0("TAD rank vs. gene rank")
my_sub <- paste0("(gene rank <= ", genes_nTop, ")")


# outFile <- file.path(outFolder, paste0("tadRank_geneRank_topGenes_densplot.", plotType))
# do.call(plotType, list(outFile, height=myHeight, width=myWidth))
outFile <- file.path(outFolder, paste0("tadRank_geneRank_topGenes_densplot.", "png"))
do.call("png", list(outFile, height=400, width=400))
par(bty="l", family=fontFamily)
densplot(
  x = plot_dt$gene_rank,
  y = plot_dt$tad_rank,
  cex = 0.7,
  xlab = "Gene rank",
  ylab = "TAD rank",
  main = my_tit,
  cex.axis = axisCex,
  cex.main = mainCex,
  cex.lab = labCex
)
mtext(side = 3, text = paste0(my_sub), cex = subCex)
foo <- dev.off()  
cat(paste0("... written: ", outFile, "\n"))

# stop("-ok")


# inDT$dataset <- paste0(inDT$hicds, "\n", inDT$exprds)
inDT$dataset <- file.path(inDT$hicds,  inDT$exprds)

nSignif_dt <- do.call(rbind, by(inDT, inDT$dataset, function(sub_dt) {
    
  tad_genes <- sub_dt$entrezID[sub_dt$tad_adjCombPval <= tadSignifThresh]
  gene_genes <- sub_dt$entrezID[sub_dt$adj.P.Val <= geneSignifThresh]
  
  tadOnly_genes <- setdiff(tad_genes, gene_genes)
  tadAndGenes_genes <- intersect(tad_genes, gene_genes)
  geneOnly_genes <- setdiff(gene_genes, tad_genes)
  
  stopifnot(length(tadOnly_genes) + length(tadAndGenes_genes) + length(geneOnly_genes) == length(unique(c(tad_genes, gene_genes))))
  
  data.frame(
    `TAD-level`=length(tadOnly_genes),
    `gene-level`=length(geneOnly_genes),
    `TAD+gene-level`=length(tadAndGenes_genes), check.names = FALSE)
  
}))
nSignif_dt_s <- nSignif_dt

signif_order <- colnames(nSignif_dt_s)

tmp_dt <- nSignif_dt_s
tmp_dt$nTot <- rowSums(tmp_dt)
tmp_dt <- tmp_dt[order(tmp_dt$nTot, decreasing = FALSE),]
tmp_dt$dataset <- paste0(dirname(rownames(tmp_dt)), "\n", basename(rownames(tmp_dt))) 
ntot_ds_order <- tmp_dt$dataset


mycols_ntot <- all_cols[all_cmps[basename(rownames(tmp_dt))]]

nDS <- length(ntot_ds_order)

countCmp <- setNames(as.numeric(table(all_cmps[basename(rownames(tmp_dt)) ])), names(table(all_cmps[basename(rownames(tmp_dt)) ])))
legTitle <- "Signif."
legDT <- data.frame(cmpType = names(all_cols), cmpColor = as.character(all_cols))
legDT <- legDT[rev(order(legDT$cmpType)),]
legDT$count <- countCmp[legDT$cmpType]
legDT$legLabel <- paste0(legDT$cmpType, " (", legDT$count, ")")

my_sub <- paste0("gene signif.: p-val <= ", geneSignifThresh, "; TAD signif.: p-val <= ", tadSignifThresh, "\n")






# 
# 
# ##############################################################
# ############################### barplot fract.
# ##############################################################
# 
# nSignif_dt <- data.frame(t(apply(nSignif_dt_s, 1, FUN=function(x) x/sum(x))), check.names = FALSE)
# stopifnot(abs(rowSums(nSignif_dt) - 1) < 1e-10)
# 
# nSignif_dt$dataset <- paste0(dirname(rownames(nSignif_dt)), "\n", basename(rownames(nSignif_dt)))
# 
# plot_dt <- melt(nSignif_dt, id = "dataset")
# plot_dt$dataset <- factor(plot_dt$dataset, levels=ntot_ds_order)
# plot_dt$labSymb <- labsymbol
# 
# my_tit <- "Fract. signif. genes TAD-level/gene-level"
# 
# plot_dt$variable <- factor(as.character(plot_dt$variable), levels = signif_order)
# stopifnot(!is.na(plot_dt$variable))
# 
# signif_fract_plot_tmp <- ggplot(plot_dt, aes(x=dataset, y=value, fill=variable, color=variable)) + 
#   geom_bar(position="stack", stat="identity") +
#   ggtitle(paste0(my_tit), 
#           subtitle = paste0(my_sub) )+
#   labs(fill=paste0(legTitle)) + 
#   guides(color=FALSE)+
#   # coord_cartesian(expand=FALSE)+
#   coord_cartesian(clip = 'off', expand=FALSE) +
#   scale_y_continuous(name=paste0("Fraction of signif. genes"),
#                      limits = c(0,1), 
#                      breaks = seq(from=0, to=1, by=0.1),
#                      labels = seq(from=0, to=1, by=0.1))+
#   eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")")))+
#   eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")")))+
#   theme( # Increase size of axis lines
#     plot.margin = unit(plotMargin, "lines"), # top, right, bottom, and left 
#     plot.title = element_text(hjust = 0.5, face = "bold", size=16, family=fontFamily),
#     plot.subtitle = element_text(hjust = 0, vjust=2,face = "italic", size =8, family=fontFamily, lineheight = 1.75),
#     panel.grid = element_blank(),
#     # panel.grid.major.y = element_line(colour = "grey"),
#     # panel.grid.minor.y = element_line(colour = "grey"),
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     axis.line.x = element_line(size = .2, color = "black"),
#     axis.line.y = element_line(size = .2, color = "black"),
#     axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12, family=fontFamily),
#     axis.title.y = element_text(color="black", size=14, family=fontFamily),
#     axis.title.x = element_text(color="black", size=14, family=fontFamily),
#     panel.border = element_blank(),
#     panel.background = element_rect(fill = "transparent"),
#     legend.background =  element_rect(),
#     legend.key = element_blank(),
#     legend.title = element_text(face="bold")
#   ) +
#   geom_text(data=tmp_dt, aes(x = tmp_dt$dataset, y=1, 
#                         label=sprintf("%.0f", tmp_dt$nTot)),
#             inherit.aes=FALSE, angle=90, size=3, 
#             vjust=0.5, hjust=0)+
#   theme(
#     # legend.position = c(.95, .95),
#     # legend.box.just = "right",
#     # legend.margin = margin(6, 6, 6, 6),
#     legend.justification = c("right", "top")
#   )+ 
#   geom_text(data=legDT, aes(label = legDT$legLabel, x = 59, y =c(0, 0.05, 0.1)),
#             vjust = 0, hjust=0,
#             inherit.aes = FALSE, color = legDT$cmpColor)
# 
# 
# 
# 
# signif_fract_plot_lab <- signif_fract_plot_tmp +
#   scale_x_discrete(name=paste0("(all datasets - n=", nDS, ")")) +
#   theme(    axis.text.x = element_text(color=mycols_ntot, hjust=1,vjust = 0.5, size=7, angle=90, family=fontFamily) )
#   
# 
# outFile <- file.path(outFolder, paste0("all_ds_fract_signif_genes_withLabs_geneSignif",geneSignifThresh, "_tadSignif", tadSignifThresh, "_barplot.", plotType))
# ggsave(plot = signif_fract_plot_lab, filename = outFile, height=myHeightGG*1.5, width = myWidthGG*2)
# cat(paste0("... written: ", outFile, "\n"))
# 
# 
# signif_fract_plot_symb <- signif_fract_plot_tmp +
#   scale_x_discrete(labels=plot_dt$labSymb, name=paste0("(all datasets - n=", nDS, ")")) +
#   theme(    axis.text.x = element_text(color=mycols_ntot, hjust=0.5,vjust = 0.5, size=7, angle=90, family=fontFamily) )
#   
#   
# outFile <- file.path(outFolder, paste0("all_ds_fract_signif_genes_withSymb_geneSignif",geneSignifThresh, "_tadSignif", tadSignifThresh, "_barplot.", plotType))
# ggsave(plot = signif_fract_plot_symb, filename = outFile, height=myHeightGG, width = myWidthGG*2)
# cat(paste0("... written: ", outFile, "\n"))
# 
# ##############################################################
# ############################### barplot nbr
# ##############################################################
# nSignif_dt <- nSignif_dt_s
# nSignif_dt$dataset <- paste0(dirname(rownames(nSignif_dt)), "\n", basename(rownames(nSignif_dt)))
# 
# plot_dt <- melt(nSignif_dt, id = "dataset")
# 
# plot_dt$dataset <- factor(plot_dt$dataset, levels=ntot_ds_order)
# plot_dt$labSymb <- labsymbol
# 
# my_tit <- "# signif. genes TAD-level/gene-level"
# plot_dt$variable <- factor(as.character(plot_dt$variable), levels = signif_order)
# stopifnot(!is.na(plot_dt$variable))
# 
# 
# signif_nbr_plot_tmp <- ggplot(plot_dt, aes(x=dataset, y=value, fill=variable, color=variable)) + 
#   geom_bar(position="stack", stat="identity") +
#   ggtitle(paste0(my_tit), 
#           subtitle = paste0(my_sub) ) +
#   labs(fill=paste0(legTitle)) + 
#   guides(color=FALSE)+
#   # coord_cartesian(expand=FALSE)+
#   coord_cartesian(clip = 'off', expand=FALSE) +
#   # scale_y_continuous(name=paste0("# of signif. genes"),
#   #                    limits = c(0,1), 
#   #                    breaks = seq(from=0, to=1, by=0.1),
#   #                    labels = seq(from=0, to=1, by=0.1))+
#   scale_y_continuous(name=paste0("# of signif. genes"), breaks = scales::pretty_breaks(n = 10))+
#   eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")")))+
#   eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")")))+
#   theme( # Increase size of axis lines
#     plot.margin = unit(plotMargin, "lines"), # top, right, bottom, and left 
#     plot.title = element_text(hjust = 0.5, face = "bold", size=16, family=fontFamily),
#     plot.subtitle = element_text(hjust = 0, vjust=2,face = "italic", size =8, family=fontFamily, lineheight = 1.75),
#     panel.grid = element_blank(),
#     panel.grid.major.y = element_line(colour = "grey"),
#     panel.grid.minor.y = element_line(colour = "grey"),
#     # panel.grid.major.y = element_blank(),
#     # panel.grid.minor.y = element_blank(),
#     axis.line.x = element_line(size = .2, color = "black"),
#     axis.line.y = element_line(size = .2, color = "black"),
#     axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12, family=fontFamily),
#     axis.title.y = element_text(color="black", size=14, family=fontFamily),
#     axis.title.x = element_text(color="black", size=14, family=fontFamily),
#     panel.border = element_blank(),
#     panel.background = element_rect(fill = "transparent"),
#     legend.background =  element_rect(),
#     legend.key = element_blank(),
#     legend.title = element_text(face="bold")
#   ) +
#   geom_text(data=tmp_dt, aes(x = tmp_dt$dataset, y=max(tmp_dt$nTot), 
#                         label=sprintf("%.0f", tmp_dt$nTot)),
#             inherit.aes=FALSE, angle=90, size=3, 
#             vjust=0.5, hjust=0)+
#   theme(
#     # legend.position = c(.95, .95),
#     # legend.box.just = "right",
#     # legend.margin = margin(6, 6, 6, 6),
#     legend.justification = c("right", "top")
#   )+ 
#   geom_text(data=legDT, aes(label = legDT$legLabel, x = 59, y =c(0, 500, 1000)),
#             vjust = 0, hjust=0,
#             inherit.aes = FALSE, color = legDT$cmpColor)# +
#   # geom_hline(yintercept=7462)+
#   # geom_hline(yintercept=487)
# 
# signif_nbr_plot_withLabs <- signif_nbr_plot_tmp +
#   scale_x_discrete(name=paste0("(all datasets - n=", nDS, ")")) +
#   theme(    axis.text.x = element_text(color=mycols_ntot, hjust=1,vjust = 0.5, size=7, angle=90, family=fontFamily) )
# 
# outFile <- file.path(outFolder, paste0("all_ds_nbr_signif_genes_withLabs_geneSignif",geneSignifThresh, "_tadSignif", tadSignifThresh, "_barplot.", plotType))
# ggsave(plot = signif_nbr_plot_withLabs, filename = outFile, height=myHeightGG*1.5, width = myWidthGG*2)
# cat(paste0("... written: ", outFile, "\n"))
# 
# 
# signif_nbr_plot_withSymb <- signif_nbr_plot_tmp+
#   scale_x_discrete(labels=plot_dt$labSymb, name=paste0("(all datasets - n=", nDS, ")")) +
#   theme(    axis.text.x = element_text(color=mycols_ntot, hjust=0.5,vjust = 0.5, size=7, angle=90, family=fontFamily) )
# 
# outFile <- file.path(outFolder, paste0("all_ds_nbr_signif_genes_withSymb_geneSignif",geneSignifThresh, "_tadSignif", tadSignifThresh, "_barplot.", plotType))
# ggsave(plot = signif_nbr_plot_withSymb, filename = outFile, height=myHeightGG, width = myWidthGG*2)
# cat(paste0("... written: ", outFile, "\n"))
# 
