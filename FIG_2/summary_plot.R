
# Rscript summary_plot.R 

require(doMC)
require(foreach)
require(ggplot2)
require(ggsci)
require(colorRamps)
require(reshape2)

require(ggpubr)

setDir <- ""

registerDoMC(50)

plotType <- "svg"

source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")

settingFolder <- file.path(runFolder, "PIPELINE", "INPUT_FILES")

myWidthGG <- 9
myHeightGG <- 9

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
count_dt_s <-count_dt
########3


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

########
count_dt <-count_dt_s 



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

max_y <- max(count_dt$nSignifTADs)

my_pal <- "nord::victory_bonds"
# my_pal <- "wesanderson::FantasticFox1"
# withsignif <- paletteer::paletteer_d(my_pal)[3]
withsignif <- "#E19600FF" 
# nosignif <- paletteer::paletteer_d(my_pal)[4] 
nosignif <- "#193264FF "

withsignif <- "#E19600FF" 
nosignif <- "#193264FF"


plot_names <- c("nSignifTADs_withSignifGenes", "nSignifTADs_noSignifGenes")

my_cols <- setNames(c(withsignif, nosignif), plot_names)
my_cols_names <- setNames(c("with signif. genes", "without any signif. genes"), plot_names)

plotTit <- ""
subTit <- ""
plotTit <- paste0("Differentially activated TADs and signif. genes")
subTit <- paste0("gene adj. p-val <= ", geneSignifThresh, "; TAD adj. p-val <= ", tadSignifThresh)

# m_count_dt$ds_nbr <- as.numeric(m_count_dt$dataset)
# m_count_dt$ds_nbr2 <- m_count_dt$ds_nbr *2

break_step <- 5
samp_y_start_offset <- 4 + break_step
samp_y_start <- 70
samp_axis_offset <- 30
axis_lim <- 65

text_label_dt <- do.call(rbind, lapply(ds_order, function(x) {
  settingFile <- file.path(settingFolder, dirname(x), paste0("run_settings_", basename(x), ".R"))
  stopifnot(file.exists(settingFile))
  source(settingFile)
  samp1 <- get(load(file.path(setDir, sample1_file)))
  samp2 <- get(load(file.path(setDir, sample2_file)))
  data.frame(
    dataset=x,
    n_samp1=length(samp1),
    n_samp2=length(samp2),
    cond1=cond1,
    cond2=cond2,
    stringsAsFactors = FALSE
  )
}))
text_label_dt$dataset <- factor(text_label_dt$dataset, levels=ds_order)
stopifnot(!is.na(text_label_dt$dataset))
text_label_dt$y_pos <- max_y + samp_y_start_offset
text_label_dt$y_pos <- samp_y_start
text_label_dt$x_pos <- as.numeric(text_label_dt$dataset)
# text_label_dt$samp_lab <- paste0("(", text_label_dt$n_samp1, " - ", text_label_dt$n_samp2, ")")
text_label_dt$samp_lab <- paste0("(", text_label_dt$n_samp1, " vs. ", text_label_dt$n_samp2, ")")


add_axes <- function(p) {
  return(
    p +   geom_text(data=text_label_dt, aes(x=x_pos, y= y_pos, label=samp_lab), inherit.aes = FALSE, hjust=0) +
      scale_y_continuous(
        breaks=seq(from=break_step, to=samp_y_start+break_step, by=break_step), 
        labels = c(seq(from=break_step, to=axis_lim, by=break_step),"",  "\t\t(# samp.)"), 
        # breaks=seq(from=5, to=axis_lim, by=5), 
        expand=c(0,0), 
                         limits=c(0, max_y+ samp_axis_offset)) +
      with_samp_theme +
      geom_hline(yintercept = seq(from=break_step, axis_lim, by=break_step), color="darkgrey") +
      geom_segment(x= 0.5, xend=0.5, yend=axis_lim, y=0, color="darkgrey")
    
  )
}

no_samp_theme <- theme(
  plot.subtitle = element_text(hjust=0.5, face="italic"),
  panel.grid.major.y =element_blank(),
  panel.grid.major.x =element_line(color="darkgrey"),
  panel.grid.minor.x =element_line(color="darkgrey"),
  legend.position = "top",
  legend.text = element_text(size=14),
  axis.title.x = element_text(),
  axis.text.y = element_text(hjust=1, size=8),
  axis.line.x = element_line(colour="darkgrey"),
  axis.title.y = element_blank()
) 

p_signif <- ggplot(m_count_dt[m_count_dt$variable %in% c(plot_names),], 
       # aes(x=dataset, y = value, fill=variable))+ 
    aes(x=dataset_name, y = value, fill=variable))+ 
  ggtitle(plotTit, subtitle = paste0(subTit))+
  geom_bar(stat="identity", position="stack") +
  scale_fill_manual(values=my_cols, labels=my_cols_names)+
  # labs(fill = "Signif. TADs:", x = "", y="# of TADs") +
  labs(fill = "", x = "", y="# of TADs") +
  # scale_x_discrete(labels=function(x) gsub("/", "\n", x)) +
  scale_x_discrete(labels=function(x) gsub("/", " - ", x)) +
  # scale_y_continuous(breaks=seq(from=10, 70, by=5), expand=c(0,0))+
  coord_flip()+
  my_box_theme + 
  no_samp_theme


with_samp_theme <-  theme(
  plot.subtitle = element_text(hjust=0.5, face="italic"),
  axis.line.x =element_blank(),
  panel.grid.major.y =element_blank(),
  panel.grid.major.x =element_blank(),
  panel.grid.minor.x =element_blank()
) 

p_signif2 <- add_axes(p_signif)
# + geom_text(data=text_label_dt, aes(x=x_pos, y= y_pos, label=samp_lab), inherit.aes = FALSE, hjust=0) +
#   scale_y_continuous(breaks=seq(from=5, 65, by=5), expand=c(0,0),
#                      limits=c(0, max_y + samp_axis_offset)) +
#   with_samp_theme +
#   geom_hline(yintercept = seq(from=5, 65, by=5), color="darkgrey")
# p_signif2 <- p_signif2+geom_segment(x= 0.5, xend=0.5, yend=65, y=0, color="darkgrey")

outFile <- file.path(outFolder, paste0("summary_plot_signifTADs_signifGenes.", plotType))
ggsave(p_signif2, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

##############################################################
##############################################################

plot_names <- c("nSignifTADs_withTopGenes", "nSignifTADs_noTopGenes")

my_cols <- setNames(c(withsignif, nosignif), plot_names)
my_cols_names <- setNames(c("with top genes", "without any top genes"), plot_names)

plotTit <- ""
subTit <- ""
plotTit <- paste0("Differentially activated TADs and top DE genes")
subTit <- paste0("gene rank <= ", nTopGenes, "; TAD adj. p-val <= ", tadSignifThresh)

p_signif <- ggplot(m_count_dt[m_count_dt$variable %in% c(plot_names),], 
                   # aes(x=dataset, y = value, fill=variable))+ 
                   aes(x=dataset_name, y = value, fill=variable))+ 
  ggtitle(plotTit, subtitle = paste0(subTit))+
  geom_bar(stat="identity", position="stack") +
  scale_fill_manual(values=my_cols, labels=my_cols_names)+
  # labs(fill = "Signif. TADs:", x = "", y="# of TADs") +
  labs(fill = "", x = "", y="# of TADs") +
  # scale_x_discrete(labels=function(x) gsub("/", "\n", x)) +
  scale_x_discrete(labels=function(x) gsub("/", " - ", x)) +
  # scale_y_continuous(breaks=seq(from=10, 70, by=5), expand=c(0,0))+
  coord_flip()+
  my_box_theme + 
  no_samp_theme



p_signif2 <- add_axes(p_signif)


outFile <- file.path(outFolder, paste0("summary_plot_signifTADs_topGenes.", plotType))
ggsave(p_signif2, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

#####################








