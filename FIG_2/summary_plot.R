
# Rscript summary_plot.R 

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

# geneSignifThresh and tadSignifThresh loaded from settings.R

nTopGenes <- 100

all_dt2 <- get(load(file.path(runFolder, "CREATE_FINAL_TABLE", "all_result_dt.Rdata")))
signif_dt2 <- all_dt2[all_dt2$adjPvalComb <= tadSignifThresh,]
signif_dt2$tad_id <- file.path(signif_dt2$hicds, signif_dt2$exprds, signif_dt2$region)
signif_dt2$dataset <- file.path(signif_dt2$hicds, signif_dt2$exprds)
signif_dt2 <- unique(signif_dt2)
nSignif_dt <- aggregate(tad_id~dataset, data=signif_dt2, FUN=length)
colnames(nSignif_dt)[colnames(nSignif_dt) == "tad_id"] <- "nSignifTADs"

all_dt <- get(load(file.path(runFolder, "GENE_RANK_TAD_RANK", "all_gene_tad_signif_dt.Rdata")))
all_dt$tad_id <- file.path(all_dt$hicds, all_dt$exprds, all_dt$region)
all_dt$dataset <- file.path(all_dt$hicds, all_dt$exprds)
tadSignif_dt <- all_dt[all_dt$tad_adjCombPval <= tadSignifThresh,]

tadSignif_geneSignif_dt <- all_dt[all_dt$tad_adjCombPval <= tadSignifThresh & all_dt$adj.P.Val <= geneSignifThresh,]
geneSignif_nSignif_dt <- aggregate(tad_id~dataset, data=tadSignif_geneSignif_dt, FUN=function(x)length(unique(x)))
colnames(geneSignif_nSignif_dt)[colnames(geneSignif_nSignif_dt) == "tad_id"] <- "nSignifTADs_withSignifGenes"

tadSignif_geneTop_dt <- all_dt[all_dt$tad_adjCombPval <= tadSignifThresh & all_dt$gene_rank <= nTopGenes,]
geneTop_nSignif_dt <- aggregate(tad_id~dataset, data=tadSignif_geneTop_dt, FUN=function(x)length(unique(x)))
colnames(geneTop_nSignif_dt)[colnames(geneTop_nSignif_dt) == "tad_id"] <- "nSignifTADs_withTopGenes"

count_dt <- merge(nSignif_dt,  merge(geneSignif_nSignif_dt, geneTop_nSignif_dt, by="dataset", all=T),by="dataset", all=T)
count_dt[is.na(count_dt)] <- 0
count_dt$nSignifTADs_noSignifGenes <- count_dt$nSignifTADs-count_dt$nSignifTADs_withSignifGenes
count_dt$nSignifTADs_noTopGenes <- count_dt$nSignifTADs-count_dt$nSignifTADs_withTopGenes

stopifnot(setequal(signif_dt2$tad_id, tadSignif_dt$tad_id))

outFolder <- "SUMMARY_PLOT"
dir.create(outFolder, recursive = TRUE)

buildData <- FALSE


count_dt <- count_dt[order(count_dt$nSignifTADs),]
ds_order <- as.character(count_dt$dataset)

m_count_dt <- melt(count_dt, id="dataset")
m_count_dt$dataset <- factor(m_count_dt$dataset, levels = ds_order)
stopifnot(!is.na(m_count_dt$dataset))


ggplot(m_count_dt[m_count_dt$variable %in% c("nSignifTADs_withSignifGenes", "nSignifTADs_noSignifGenes"),], aes(x=dataset, y = value, 
  geom_bar(stat="identity", position="stack") +
  scale_x_discrete(labels=function(x) gsub("/", "\n", x)) +
  coord_flip()+
  my_box_theme + 
  theme(legend.position = "top",
        axis.title.x = element_blank())





















ggsci_pal <- "lancet"
ggsci_subpal <- ""

all_hicds <- all_obs_hicds
all_exprds <- all_obs_exprds

# all_hicds=all_hicds[1:2]

fcc_thresh <- 1