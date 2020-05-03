

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

myWidth <- myWidth * 1.2

outFolder <- "TAD_VS_GENE_SIGNIF_ALL_DS"
dir.create(outFolder, recursive = TRUE)

buildData <- TRUE

ggsci_pal <- "d3"
ggsci_subpal <- ""

plotMargin <- c(1,2,1,1)


inDT <- get(load("../../v2_Yuanlong_Cancer_HiC_data_TAD_DA/GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"))

tmp <- inDT[inDT$hicds==hicds & inDT$exprds== exprds,  ]

##############################################################
############################### PREP DATA FOR ALL CATEGORIES
##############################################################
tads_nTop <- 10

topGenes <- 100
lastGenes <- 100

genes_nTop <- 100

topTADs <- 10
lastTADs <- 10

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

##############################################################
############################### PREP DATA FOR ALL CATEGORIES
##############################################################
agg_dt <- aggregate(entrezID~gene_rank_cat + tad_rank_cat, data = inDT_dt, FUN=length)
colnames(agg_dt)[colnames(agg_dt) == "entrezID"] <- "nGenes"

ggballoonplot(agg_dt, x="gene_rank_cat", y ="tad_rank_cat", fill = "nGenes")+
  scale_fill_viridis_c(option = "C", trans="log")

agg_dt2 <- agg_dt[!(agg_dt$gene_rank_cat=="other" & agg_dt$tad_rank_cat=="other"),] 
ggballoonplot(agg_dt2, x="gene_rank_cat", y ="tad_rank_cat", fill = "nGenes")+
  scale_fill_viridis_c(option = "C", trans="log")

inDT_dt$foo_val <- TRUE

agg_dt3 <- aggregate(entrezID~gene_rank+tad_rank +dataset, FUN=length,data=inDT_dt)

ggplot(agg_dt3, aes(x = gene_rank, y = tad_rank, z =entrezID)) +
  # stat_summary_hex( function(z) log(sum(z))) +
  stat_summary_hex() +
  scale_fill_viridis_c(option = "C", trans="log")
  labs(fill = "Log counts")


##############################################################
############################### stacked barplot top10 TADs
##############################################################

tads_nTop <- 10

topGenes <- 100
lastGenes <- 100

plot_dt <- inDT_dt[inDT_dt$tad_rank <= tads_nTop,]

  
stopifnot(!duplicated(plot_dt))
agg_dt <- aggregate(entrezID ~ tad_rank+gene_rank_cat, data = plot_dt, FUN=length)
colnames(agg_dt)[colnames(agg_dt) == "entrezID"] <- "nGenes"

ggplot(agg_dt, aes(x=tad_rank, y = nGenes, color = gene_rank_cat, fill = gene_rank_cat)) +
  geom_bar(stat="identity")



##############################################################
############################### geom area  top10 TADs
##############################################################

rankTieMeth <- "min"

genes_nTop <- 100

topTADs <- 10
lastTADs <- 10



plot_dt <- inDT_dt[inDT_dt$gene_rank <= genes_nTop,]



stopifnot(!duplicated(plot_dt))
agg_dt <- aggregate(entrezID ~ gene_rank+tad_rank_cat, data = plot_dt, FUN=length)
colnames(agg_dt)[colnames(agg_dt) == "entrezID"] <- "nTADs"

ggplot(agg_dt, aes(x=gene_rank, y = nTADs, color = tad_rank_cat, fill = tad_rank_cat)) +
  geom_bar(stat="identity")




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


outFile <- file.path(outFolder, paste0("tadRank_geneRank_topGenes_densplot.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
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
