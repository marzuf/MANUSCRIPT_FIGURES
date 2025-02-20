

# Rscript tad_gene_signif_missedTADs.R

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

outFolder <- "TAD_GENE_SIGNIF_MISSEDTADS"
dir.create(outFolder, recursive = TRUE)

buildData <- TRUE

ggsci_pal <- "d3"
ggsci_subpal <- ""

plotMargin <- c(1,2,1,1)

myHeightGG <- 5
myWidthGG <- 7

rankTieMeth <- "min"

nTopTADs <- 10
nTopGenes <- 100

yOffset <- 0.4 # for the balloonplot

inDT <- get(load(file.path(runFolder, "GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata")))
inDT$tad_id <- file.path(inDT$hicds, inDT$exprds, inDT$region)
tadRanks <- setNames(inDT$tad_rank, as.character(inDT$tad_id))
nDS <- length(unique(dirname(inDT$tad_id)))

baloonMin <- 3
baloonMax <- 15

geneSignifThresh <- 0.01
tadSignifThresh <- 0.01

##############################################################
############################### barplot missed top TADs - top genes
##############################################################

tmp_dt <- inDT[inDT$tad_rank <= nTopTADs,]
ntoprank_dt <- aggregate(tad_id~hicds+exprds, data = tmp_dt, FUN=function(x)length(unique(x)))
ntoprank_dt <- ntoprank_dt[order(ntoprank_dt$tad_id, decreasing = TRUE),]
ntoprank_dt$dataset <- file.path(ntoprank_dt$hicds, ntoprank_dt$exprds)
ntoptads_byDS <- setNames(ntoprank_dt$tad_id, ntoprank_dt$dataset)

# found the top TADs that have no top 100 genes:

# top 10 tads
top10_tads <- unique(inDT$tad_id[inDT$tad_rank <= nTopTADs])

# tads of the top 100 genes
top100genes_tads <- unique(inDT$tad_id[inDT$gene_rank <= nTopGenes])

top_missed_dt <- data.frame(missed_id=setdiff(top10_tads, top100genes_tads), stringsAsFactors=FALSE)
top_missed_dt$dataset <- dirname(top_missed_dt$missed_id)
stopifnot(!duplicated(top_missed_dt$missed_id))

top_agg_missed_dt <- aggregate(missed_id ~ dataset, data=top_missed_dt, FUN=function(x)length(unique(x)))

stopifnot(inDT$tad_rank[inDT$tad_id %in% top_missed_dt$missed_id] <= nTopTADs)
stopifnot(inDT$gene_rank[inDT$tad_id %in% top_missed_dt$missed_id] > nTopGenes)

plotTit <- paste0("Top ", nTopTADs, " TADs without any top ", nTopGenes, " genes")
subTit <- paste0("all datasets - n=", nDS)

top_agg_missed_dt <- top_agg_missed_dt[order(top_agg_missed_dt$missed_id, decreasing = T),]
top_agg_missed_dt$dataset <- factor(top_agg_missed_dt$dataset, levels=top_agg_missed_dt$dataset)
top_agg_missed_dt$cmpType <- all_cmps[paste0(basename(as.character(top_agg_missed_dt$dataset)))]
top_agg_missed_dt$dsCols <- all_cols[all_cmps[paste0(basename(as.character(top_agg_missed_dt$dataset)))]]

top_top_bar <- ggbarplot(top_agg_missed_dt, x="dataset", y="missed_id", fill = "cmpType" ) + 
  ggtitle(plotTit, subtitle = subTit)+
  scale_y_continuous(name="# \"missed\" TADs", breaks = scales::pretty_breaks(n = 8), expand=c(0,0))+
  scale_fill_manual(values = all_cols, labels = cmp_names)+
  labs(fill="", x="")+
  my_box_theme + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

outFile <- file.path(outFolder, 
                     paste0("topTADsWithoutTopGenes_nTopTADs", nTopTADs, "_nTopGenes", nTopGenes, "_barplot.", plotType))
ggsave(plot = top_top_bar, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

save(top_agg_missed_dt, file="top_agg_missed_dt.Rdata", version=2)
save(ntoptads_byDS, file="ntoptads_byDS", version=2)

top_agg_missed_dt$nTopTADs <- ntoptads_byDS[paste0(top_agg_missed_dt$dataset)]
stopifnot(!is.na(top_agg_missed_dt$nTopTADs))
top_agg_missed_dt$missed_id_ratio <- top_agg_missed_dt$missed_id/top_agg_missed_dt$nTopTADs 

save(top_agg_missed_dt, file="top_agg_missed_dt.Rdata", version=2)

top_agg_missed_dt <- top_agg_missed_dt[order(top_agg_missed_dt$missed_id_ratio, decreasing = T),]
top_agg_missed_dt$dataset <- factor(as.character(top_agg_missed_dt$dataset), levels=as.character(top_agg_missed_dt$dataset))

top_top_bar <- ggbarplot(top_agg_missed_dt, x="dataset", y="missed_id_ratio", fill = "cmpType" ) + 
  ggtitle(plotTit, subtitle = subTit)+
  scale_y_continuous(name="Ratio \"missed\" TADs", breaks = scales::pretty_breaks(n = 8), expand=c(0,0))+
  scale_fill_manual(values = all_cols, labels = cmp_names)+
  labs(fill="", x="")+
  my_box_theme + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

outFile <- file.path(outFolder, 
                     paste0("topTADsWithoutTopGenes_ratio_nTopTADs", nTopTADs, "_nTopGenes", nTopGenes, "_barplot.", plotType))
ggsave(plot = top_top_bar, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


# 
xxx=tmp_dt[tmp_dt$hicds=="ENCSR489OCU_NCI-H460_40kb" & tmp_dt$exprds=="TCGAluad_norm_luad",]
unique(xxx[,c("tad_id", "tad_rank")])
xxx= inDT[ inDT$hicds=="ENCSR489OCU_NCI-H460_40kb" &  inDT$exprds=="TCGAluad_norm_luad" &
             inDT$tad_rank <= nTopTADs,]
stopifnot(sum(aggregate(gene_rank ~ hicds+exprds+region, data=xxx, FUN=min)[,"gene_rank"] > nTopGenes) == 
  top_agg_missed_dt$missed_id[top_agg_missed_dt$dataset==file.path("ENCSR489OCU_NCI-H460_40kb","TCGAluad_norm_luad")])

##############################################################
############################### min gene rank of signif. tads
##############################################################
top_inDT <- inDT[inDT$tad_rank <= nTopTADs,]
min_rank_tad_dt <- aggregate(gene_rank~hicds+exprds+region+tad_id+tad_rank, data=top_inDT, FUN=min)
min_rank_tad_dt$dataset <- file.path(min_rank_tad_dt$hicds, min_rank_tad_dt$exprds)

tmp_mean <- aggregate(gene_rank~dataset, data=min_rank_tad_dt, FUN=mean)
tmp_mean <- tmp_mean[order(tmp_mean$gene_rank, decreasing = TRUE),]

min_rank_tad_dt$dataset <- factor(as.character(min_rank_tad_dt$dataset), levels=as.character(tmp_mean$dataset))
min_rank_tad_dt$cmpType <- all_cmps[paste0(basename(as.character(min_rank_tad_dt$dataset)))]

plotTit <- paste0("Gene best rank of top ", nTopTADs, " TADs")

min_geneRank_box_byDS <- ggplot(min_rank_tad_dt, aes(x=dataset, y=gene_rank, fill = cmpType )) + 
  geom_boxplot()+
  ggtitle(plotTit, subtitle = subTit)+
  scale_y_continuous(name="Min gene rank of top TADs", breaks = scales::pretty_breaks(n = 8))+
  scale_fill_manual(values = all_cols, labels = cmp_names)+
  labs(fill="", x="")+
  my_box_theme + 
  theme(
    axis.line=element_line(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

outFile <- file.path(outFolder, 
                     paste0("minGeneRankTopTADs_nTopTADs", nTopTADs, "_nTopGenes", nTopGenes, "_byDS_barplot.", plotType))
ggsave(plot = min_geneRank_box_byDS, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

min_rank_tad_dt$tad_rank <- factor(as.character(min_rank_tad_dt$tad_rank), levels=as.character(1:nTopTADs))

min_rank_tad_dt$gene_rank_log10 <- log10(min_rank_tad_dt$gene_rank)

min_geneRank_box_byRank <- ggplot(min_rank_tad_dt, aes(x=tad_rank, y=gene_rank, fill = cmpType )) + 
  geom_boxplot()+
  ggtitle(plotTit, subtitle = subTit)+
  scale_y_continuous(name="Min gene rank of top TADs", breaks = scales::pretty_breaks(n = 8))+
  scale_fill_manual(values = all_cols, labels = cmp_names)+
  labs(fill="", x="")+
  my_box_theme + 
  theme(
    axis.line=element_line(),
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

outFile <- file.path(outFolder, 
                     paste0("minGeneRankTopTADs_nTopTADs", nTopTADs, "_nTopGenes", nTopGenes, "_byRank_barplot.", plotType))
ggsave(plot = min_geneRank_box_byRank, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


min_geneRank_box_byRankLog10 <- ggplot(min_rank_tad_dt, aes(x=tad_rank, y=gene_rank_log10, fill = cmpType )) + 
  geom_boxplot()+
  ggtitle(plotTit, subtitle = subTit)+
  scale_y_continuous(name="Min gene rank of top TADs", breaks = scales::pretty_breaks(n = 8))+
  scale_fill_manual(values = all_cols, labels = cmp_names)+
  labs(fill="", x="")+
  my_box_theme + 
  theme(
    axis.line=element_line(),
    # axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

outFile <- file.path(outFolder, 
                     paste0("minGeneRankTopTADs_nTopTADs", nTopTADs, "_nTopGenes", nTopGenes, "_byRank_log10_barplot.", plotType))
ggsave(plot = min_geneRank_box_byRankLog10, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

plotTit <- paste0("Gene best rank of top ", nTopTADs, " TADs")

min_geneRank_box_byDSLog10 <- ggplot(min_rank_tad_dt, aes(x=dataset, y=gene_rank_log10, fill = cmpType )) + 
  geom_boxplot()+
  ggtitle(plotTit, subtitle = subTit)+
  scale_y_continuous(name="Min gene rank of top TADs", breaks = scales::pretty_breaks(n = 8))+
  scale_fill_manual(values = all_cols, labels = cmp_names)+
  labs(fill="", x="")+
  my_box_theme + 
  theme(
    axis.line=element_line(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

outFile <- file.path(outFolder, 
                     paste0("minGeneRankTopTADs_nTopTADs", nTopTADs, "_nTopGenes", nTopGenes, "_byDS_log10_barplot.", plotType))
ggsave(plot = min_geneRank_box_byDSLog10, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



##############################################################
############################### barplot missed top TADs - top genes by rank
##############################################################

# found the top TADs that have no top 100 genes:

# top 10 tads
top10_tads <- inDT$tad_id[inDT$tad_rank <= nTopTADs]

# tads of the top 100 genes
top100genes_tads <- inDT$tad_id[inDT$gene_rank <= nTopGenes]

rank_top_missed_dt <- data.frame(missed_id=setdiff(top10_tads, top100genes_tads), stringsAsFactors=FALSE)
rank_top_missed_dt$dataset <- dirname(rank_top_missed_dt$missed_id)
rank_top_missed_dt$tad_rank <- tadRanks[paste0(rank_top_missed_dt$missed_id)]
stopifnot(!is.na(rank_top_missed_dt$tad_rank))

rank_top_agg_missed_dt <- aggregate(missed_id ~ tad_rank+dataset, data=rank_top_missed_dt, FUN=length)
rank_top_agg_missed_dt$cmpType <- all_cmps[paste0(basename(as.character(rank_top_agg_missed_dt$dataset)))]
rank_top_agg_missed_dt$dsCols <- all_cols[all_cmps[paste0(basename(as.character(rank_top_agg_missed_dt$dataset)))]]

rank_top_agg_missed_dt$tad_rank <- factor(rank_top_agg_missed_dt$tad_rank, levels=1:10)

plotTit <- paste0("Top ", nTopTADs, " TADs without any top ", nTopGenes, " genes")
subTit <- paste0("by rank - all datasets - n=", nDS)


top_top_box <- ggplot(rank_top_agg_missed_dt, aes(x=tad_rank, y=missed_id)) + 
  # geom_boxplot() +
  geom_boxplot(notch = FALSE, outlier.shape=NA, fill=NA)+
  geom_point(aes(color=cmpType),   alpha=0.5)    +
  ggtitle(plotTit, subtitle = subTit)+
  labs(color="")+
  theme(legend.key = element_rect(fill = NA))+
  scale_color_manual(values = all_cols, labels = cmp_names)+
  scale_y_continuous(name="# \"missed\" TADs", breaks = scales::pretty_breaks(n = 8))+
  scale_x_discrete(name="TAD rank", labels=1:nTopTADs)+
  my_box_theme# + 
  # theme(
  #   axis.text.x = element_blank(),
  #   axis.ticks.x = element_blank()
  # )


outFile <- file.path(outFolder, 
                     paste0("topTADsWithoutTopGenes_nTopTADs", nTopTADs, "_nTopGenes", nTopGenes, "_byRank_boxplot.", plotType))
ggsave(plot = top_top_box, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

##############################################################
############################### barplot missed signif. tads - top genes
##############################################################

# found the top TADs that have no top 100 genes:

# signif tads
signif_tads <- inDT$tad_id[inDT$tad_adjCombPval <= tadSignifThresh]
nSignif <- table(dirname(signif_tads))

# tads of the top 100 genes
top100genes_tads <- inDT$tad_id[inDT$gene_rank <= nTopGenes]

signif_missed_dt <- data.frame(missed_id=setdiff(signif_tads, top100genes_tads), stringsAsFactors=FALSE)
signif_missed_dt$dataset <- dirname(signif_missed_dt$missed_id)

signif_agg_missed_dt <- aggregate(missed_id ~ dataset, data=signif_missed_dt, FUN=length)

signif_agg_missed_dt$totSignif <- as.numeric(nSignif[paste0(signif_agg_missed_dt$dataset)])

signif_agg_missed_dt$ratioMissed <- signif_agg_missed_dt$missed_id/signif_agg_missed_dt$totSignif

signif_agg_missed_dt <- signif_agg_missed_dt[order(signif_agg_missed_dt$ratioMissed, decreasing = T),]
signif_agg_missed_dt$dataset <- factor(signif_agg_missed_dt$dataset, levels=signif_agg_missed_dt$dataset)
signif_agg_missed_dt$cmpType <- all_cmps[paste0(basename(as.character(signif_agg_missed_dt$dataset)))]
signif_agg_missed_dt$dsCols <- all_cols[all_cmps[paste0(basename(as.character(signif_agg_missed_dt$dataset)))]]


plotTit <- paste0("Signif. TADs without any top ", nTopGenes, " genes")
subTit <- paste0("TAD adj. p-val <= ", tadSignifThresh, "; all datasets - n=", nDS)


signif_top_barp <- ggbarplot(signif_agg_missed_dt, x="dataset", y="ratioMissed", fill = "cmpType" ) + 
  ggtitle(plotTit, subtitle = subTit)+
  scale_y_continuous(name="Ratio of signif. TADs", breaks = scales::pretty_breaks(n = 8), expand=c(0,0))+
  scale_fill_manual(values = all_cols, labels = cmp_names)+
  labs(fill="", x="")+
  my_box_theme + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


outFile <- file.path(outFolder, 
                     paste0("ratioSignifTADsWithoutTopGenes_tadSignif", tadSignifThresh, "_nTopGenes", nTopGenes, "_barplot.", plotType))
ggsave(plot = signif_top_barp, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



