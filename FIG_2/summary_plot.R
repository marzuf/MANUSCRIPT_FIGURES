
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

count_dt$dataset_name <- file.path(hicds_names[dirname(count_dt$dataset)], 
                                exprds_names[basename(count_dt$dataset)])
dsName_order <-  as.character(count_dt$dataset_name)

# m_count_dt <- melt(count_dt, id="dataset")
# m_count_dt$dataset <- factor(m_count_dt$dataset, levels = ds_order)
count_dt$dataset <- NULL
m_count_dt <- melt(count_dt, id="dataset_name")
m_count_dt$dataset_name <- factor(m_count_dt$dataset_name, levels = dsName_order)

# stopifnot(!is.na(m_count_dt$dataset))
stopifnot(!is.na(m_count_dt$dataset_name))

my_pal <- "nord::victory_bonds"
# my_pal <- "wesanderson::FantasticFox1"
withsignif <- paletteer::paletteer_d(my_pal)[3]
nosignif <- paletteer::paletteer_d(my_pal)[4] 

plot_names <- c("nSignifTADs_withSignifGenes", "nSignifTADs_noSignifGenes")

my_cols <- setNames(c(withsignif, nosignif), plot_names)
my_cols_names <- setNames(c("with signif. genes", "without signif. genes"), plot_names)

plotTit <- ""
subTit <- ""

# m_count_dt$ds_nbr <- as.numeric(m_count_dt$dataset)
# m_count_dt$ds_nbr2 <- m_count_dt$ds_nbr *2

ggplot(m_count_dt[m_count_dt$variable %in% c(plot_names),], 
       # aes(x=dataset, y = value, fill=variable))+ 
    aes(x=dataset_name, y = value, fill=variable))+ 
  ggtitle(plotTit, subtitle = paste0(subTit))+
  geom_bar(stat="identity", position="stack") +
  scale_fill_manual(values=my_cols, labels=my_cols_names)+
  labs(fill = "Signif. TADs:", x = "", y="# of TADs") +
  scale_x_discrete(labels=function(x) gsub("/", "\n", x)) +
  scale_y_continuous(breaks=seq(from=10, 70, by=10))+
  coord_flip()+
  my_box_theme + 
  theme(
    plot.subtitle = element_text(hjust=0, face="italic"),
    panel.grid.major.y =element_blank(),
    panel.grid.major.x =element_line(color="darkgrey"),
    panel.grid.minor.x =element_line(color="darkgrey"),
      legend.position = "top",
        axis.title.x = element_text(),
    axis.text.y = element_text(hjust=1, size=6),
    axis.line.x = element_line(colour="darkgrey"),
        axis.title.y = element_blank()
        )


#####################3

check_dt <- do.call(rbind, by(all_dt, all_dt$dataset, function(x) {
  
  nSignifTADs <- length(unique(x$region[x$tad_adjCombPval<=tadSignifThresh]  ))
  
  nSignifTADs_withSignifGenes <- length(unique(x$region[x$tad_adjCombPval<=tadSignifThresh & x$adj.P.Val <= geneSignifThresh]  ))
  
  nSignifTADs_withTopGenes <- length(unique(x$region[x$tad_adjCombPval<=tadSignifThresh & x$gene_rank <= nTopGenes]  ))
  
  nSignifTADs_noTopGenes <- length(setdiff(x$region[x$tad_adjCombPval<=tadSignifThresh],
                                                x$region[x$tad_adjCombPval<=tadSignifThresh & x$gene_rank <= nTopGenes]))
  
  nSignifTADs_noSignifGenes<- length(setdiff(x$region[x$tad_adjCombPval<=tadSignifThresh],
                                          x$region[x$tad_adjCombPval<=tadSignifThresh & x$adj.P.Val <= geneSignifThresh]))
  
  data.frame(
   dataset=unique(x$dataset),
    nSignifTADs=nSignifTADs  ,
     nSignifTADs_withSignifGenes=nSignifTADs_withSignifGenes,
      nSignifTADs_withTopGenes=nSignifTADs_withTopGenes,
      nSignifTADs_noSignifGenes=nSignifTADs_noSignifGenes,
      nSignifTADs_noTopGenes=nSignifTADs_noTopGenes,
      stringsAsFactors = FALSE
)
}))

rownames(check_dt) <- NULL
check_dt <- check_dt[order(as.character(check_dt$dataset)),]
count_dt <- count_dt[order(as.character(count_dt$dataset)),]
stopifnot(all.equal(check_dt, count_dt))



















ggsci_pal <- "lancet"
ggsci_subpal <- ""

all_hicds <- all_obs_hicds
all_exprds <- all_obs_exprds

# all_hicds=all_hicds[1:2]

fcc_thresh <- 1