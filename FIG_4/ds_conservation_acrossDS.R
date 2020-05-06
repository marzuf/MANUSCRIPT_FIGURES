options(scipen=100)

SSHFS=F

setDir <- "/media/electron"
setDir <- ""

# Rscript ds_conservation_acrossDS.R

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

script_name <- "ds_conservation_acrossDS.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(RColorBrewer)
require(reshape2)
require(ggplot2)
require(ggpubr)
# require(gplots)
registerDoMC(4)
source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")

outFolder <- "DS_CONSERVATION_ACROSSDS"
dir.create(outFolder, recursive=TRUE)

entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)
entrez2symb <- setNames(entrez2symb_dt$symbol, entrez2symb_dt$entrezID)


buildTable <- TRUE

plotType <- "svg"

source("../settings.R")

myHeightGG <- 5
myWidthGG <- 7

tad_pval <- 0.01
minMatchBp_ratio <- 0.8
minMatch_genes <- 3

inFile <- file.path(runFolder, 
                    "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2",
                    paste0("plot_matching_dt_tadsadjPvalComb", tad_pval, "_minBpRatio", minMatchBp_ratio, "_minInterGenes", minMatch_genes, ".Rdata"))

matching_dt <- get(load(inFile))

dt <- t(matching_dt)
m_dt <- melt(dt)
colnames(m_dt) <- c("region", "dataset", "conserved")
m_dt$dataset <- as.character(m_dt$dataset)
m_dt$hicds <- dirname(m_dt$dataset)
m_dt$exprds <- basename(m_dt$dataset)
m_dt$cmpType <- all_cmps[paste0(m_dt$exprds)]
m_dt$cmpCol <- all_cols[m_dt$cmpType]

agg_dt <- aggregate(conserved ~ region, data=m_dt, FUN=sum)
agg_dt <- agg_dt[order(agg_dt$conserved),]
reg_levels <- as.character(agg_dt$region)
stopifnot(agg_dt$conserved >= 2)
# agg_dt$region_rank <- 

# plot(
#   x = 1:nrow(agg_dt), 
#   y = agg_dt$conserved,
#   type="b"
# )

agg_cmp_dt <- aggregate(conserved ~ region + cmpType, data=m_dt, FUN=sum)
agg_cmp_dt$cmpCol <- all_cols[agg_cmp_dt$cmpType]

atLeast <- 10

conserv_atLeast <- agg_dt$region[agg_dt$conserved >= atLeast] 

plot_dt <- agg_cmp_dt[agg_cmp_dt$region %in% conserv_atLeast,]

plot_dt$region <- factor(plot_dt$region, levels = reg_levels)
stopifnot(!is.na(plot_dt$region))

plot_dt$cmpType_lab <- cmp_names[paste0(plot_dt$cmpType)]
stopifnot(!is.na(plot_dt$cmpType_lab))

all_cols_lab <- all_cols
names(all_cols_lab) <- cmp_names[names(all_cols_lab)]

cmp_levels <- sort(cmp_names)

stopifnot(cmp_levels %in% names(all_cols_lab))

plot_dt$region_rank <- rank(as.numeric(plot_dt$region))

my_ylab <- "# datasets conserved"
my_xlab <- "Ranked conserved regions"
myTit <- "Signif. conserved regions"
subTit <- paste0("TAD adj. p-val <= ", tad_pval, "\nconserv. match ratio bp >= ", minMatchBp_ratio, " and match genes >= ", minMatch_genes)
   

plot_dt$cmpType_lab <- factor(plot_dt$cmpType_lab, levels=cmp_levels)
stopifnot(!is.na(plot_dt$cmpType_lab))

save(plot_dt, file=file.path(outFolder, "plot_dt.Rdata"), version=2)

bar_p <- ggplot(plot_dt, aes(x=region_rank, y=conserved, fill = cmpType_lab, color=cmpType_lab)) +
  ggtitle(myTit, subtitle = subTit)+
  geom_bar(stat="identity") + 
  scale_y_continuous(name=my_ylab, breaks= scales::pretty_breaks(n = 5))+
  scale_x_continuous(name=my_xlab, breaks= noZero_breaks(n=5))+
  scale_fill_manual(values=all_cols_lab[cmp_levels], labels=cmp_levels)+
  scale_color_manual(values=all_cols_lab[cmp_levels], labels=cmp_levels)+
  labs(x=my_xlab, y=my_ylab, color="", fill="")+
  my_box_theme + 
  theme(
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x  =  element_blank(),
    # axis.line = element_line()
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    axis.line.x = element_line()
  )

outFile <- file.path(outFolder, paste0("nConserved_byRegion_allCond.", plotType))
ggsave(bar_p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



#####

allCmpTypes <- unique(all_cmps)
cmpT=allCmpTypes[1]
byCmp_dt <- foreach(cmpT = allCmpTypes, .combine='rbind') %dopar% {
  inFile <- file.path(runFolder, 
                      "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2",
                      cmpT,
                      paste0("plot_matching_dt_tadsadjPvalComb", tad_pval, "_minBpRatio", minMatchBp_ratio, "_minInterGenes", minMatch_genes, ".Rdata"))
  
  matching_dt <- get(load(inFile))
  
  dt <- t(matching_dt)
  m_dt <- melt(dt)
  colnames(m_dt) <- c("region", "dataset", "conserved")
  m_dt$dataset <- as.character(m_dt$dataset)
  m_dt$hicds <- dirname(m_dt$dataset)
  m_dt$exprds <- basename(m_dt$dataset)
  m_dt$cmpType <- all_cmps[paste0(m_dt$exprds)]
  stopifnot(as.character(m_dt$cmpType) == cmpT)
  agg_dt <- aggregate(conserved ~ region, data=m_dt, FUN=sum)
  agg_dt <- agg_dt[order(agg_dt$conserved),]
  reg_levels <- as.character(agg_dt$region)
  stopifnot(agg_dt$conserved >= 2)
  agg_dt$cmpType <- cmpT
  agg_dt$reg_rank <- 1:nrow(agg_dt)
  agg_dt
}
# byCmp_dt$reg_rank <- 1:nrow(byCmp_dt)
# byCmp_dt$reg_rank <- factor(as.character(byCmp_dt$reg_rank))

my_ylab <- "# datasets conserved"
my_xlab <- "Ranked conserved regions"
myTit <- "Signif. conserved regions"
subTit <- paste0("TAD adj. p-val <= ", tad_pval, "\nconserv. match ratio bp >= ", minMatchBp_ratio, " and match genes >= ", minMatch_genes)

byCmp_dt$cmpType_lab <- cmp_names[paste0(byCmp_dt$cmpType)]
stopifnot(!is.na(byCmp_dt$cmpType_lab))

byCmp_dt$cmpType_lab <- factor(byCmp_dt$cmpType_lab, levels=cmp_levels)
stopifnot(!is.na(byCmp_dt$cmpType_lab))

save(byCmp_dt, file=file.path(outFolder, "byCmp_dt.Rdata"), version=2)


p_wrap <- ggplot(byCmp_dt, aes(x=reg_rank, y=conserved, group=1))+
  ggtitle(myTit, subtitle = subTit)+
  geom_line()+ 
  facet_wrap(. ~ cmpType_lab, scales="free") +
  # scale_y_continuous(name=my_ylab)+
  scale_y_continuous(name=my_ylab, breaks= scales::pretty_breaks(n = 5))+
  scale_x_continuous(name=my_xlab, breaks= noZero_breaks(n=5))+
  my_box_theme + 
  theme(
    strip.text = element_text(size=10, face="bold"),# color=all_cols_lab),
    # strip.background = element_rect(color=all_cols_lab),
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x  =  element_blank(),
    # axis.line = element_line()
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    axis.line.x = element_line()
  )


outFile <- file.path(outFolder, paste0("nConserved_byRegion_wrapCond1.", plotType))
ggsave(p_wrap, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



g <- ggplot_gtable(ggplot_build(p_wrap))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- all_cols_lab[cmp_levels]
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  # curr_sub <- g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill 
  # g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- all_cols_lab[names(curr_sub)]
  k <- k+1
}

outFile <- file.path(outFolder, paste0("nConserved_byRegion_wrapCond.", plotType))
do.call(plotType, list(outFile, height=myHeightGG, width=myWidthGG))
grid::grid.draw(g)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



stop("--ok\n")


# * TRASH *#




myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.2

myHeightHeat <- myHeight * 1.8
myWidthHeat <- myWidth * 1.8

myWidthGG <- 12
myHeightGG <- 7

source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")


script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script11same_name <- "11sameNbr_runEmpPvalCombined"

tieMeth <- "min"


mainFolder <- file.path(".")
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)
file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

all_datasets <- unlist(lapply(1:length(all_exprds), function(x) file.path(names(all_exprds)[x], all_exprds[[x]])))

cat(paste0("n allDS = ", length(all_datasets), "\n"))

outFolder <- file.path("DS_CONSERVATION_ACROSSDS")
dir.create(outFolder, recursive = TRUE)


final_table_file <- file.path("CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_table_DT <- get(load(final_table_file))

final_table_DT$cmpType <- all_cmps[paste0(final_table_DT$exprds)]
stopifnot(!is.na(final_table_DT$cmpType))



tad_signif_col <- "tad_adjCombPval"
gene_signif_col <- "adj.P.Val"

signif_column <- "adjPvalComb"
signifThresh <- 0.01




tad_pval_thresh <- 0.01
gene_pval_thresh <- 0.05
tad_pval_thresh_log10 <- -log10(tad_pval_thresh)
gene_pval_thresh_log10 <- -log10(gene_pval_thresh)
file_suffix <- paste0("tad_pval_thresh", tad_pval_thresh, "_gene_pval_thresh", gene_pval_thresh)



minOverlapBpRatio <- 0.8
minIntersectGenes <- 3

nRegionLolli <- 10

cat(paste0("> signif_column\t=\t", signif_column, "\n"))
cat(paste0("> signifThresh\t=\t", signifThresh, "\n"))
cat(paste0("> minOverlapBpRatio\t=\t", minOverlapBpRatio, "\n"))
cat(paste0("> minIntersectGenes\t=\t", minIntersectGenes, "\n"))
cat(paste0("> nRegionLolli\t=\t", nRegionLolli, "\n"))


tad_plot_list <- list()

data_cmpType = "norm_vs_tumor"

for(data_cmpType in c("norm_vs_tumor", "subtypes", "wt_vs_mut", "")) {
  
  if(data_cmpType == "") {
    nDS <- length(unique(paste0(final_table_DT$hicds,final_table_DT$exprds)))  
  } else {
    nDS <- length(unique(paste0(final_table_DT$hicds[final_table_DT$cmpType==data_cmpType],final_table_DT$exprds[final_table_DT$cmpType==data_cmpType])))
  }
  stopifnot(nDS > 0)
  
  
  inFolder <- file.path("TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2", data_cmpType)
  stopifnot(dir.exists(inFolder))
  
  inFile <- file.path(inFolder, paste0(ifelse(data_cmpType=="", data_cmpType, paste0(data_cmpType, "_")), 
                                       "conserved_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio,
                                       "_minInterGenes", minIntersectGenes, ".Rdata"))
  stopifnot(file.exists(inFile))
  conserved_signif_tads <- get(load(inFile))
  
  countConserv <- abs(sort(-lengths(conserved_signif_tads)))
  stopifnot(countConserv > 1)
  
  countConserv_dt <- data.frame(countConserv)
  countConserv_dt$region <- factor(rownames(countConserv_dt), levels=rownames(countConserv_dt))
  countConserv_dt$conservRatio <- countConserv_dt$countConserv/nDS
  stopifnot(countConserv_dt$conservRatio <= 1)
  
  
  
  
  tad_plot_list[[paste0(ifelse(data_cmpType=="", "all", data_cmpType))]] <- ggbarplot(countConserv_dt, 
                                                                                      title=paste0("Conserv. signif. regions"),
                                                                                      subtitle=paste0("# DS = ", nDS, "; n = ", nrow(countConserv_dt)),
                                                                                      x = "region", 
                                                                                      y = "conservRatio",
                                                                                      xlab="TADs (sorted by conserv. ratio)",
                                                                                      ylab="conserv. ratio",
                                                                                      col="darkgrey", fill="darkgrey") +
    theme(
      axis.text.x=element_text(size=14),
      axis.text.y=element_text(size=14),
      plot.title=element_text(size=16, face = "bold"),
      axis.ticks.x = element_blank()) + 
    coord_cartesian(expand = F)
  
  
}




####################################################################################################################################### >>> prepare the data
if(buildTable) {
  cat("... start preparing data before matching \n")
  
  hicds = all_hicds[1]
  all_gene_tad_signif_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    
    hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(hicds_file))
    g2t_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    stopifnot(nrow(g2t_dt) > 0 )
    
    
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      
      regionList_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_regionList.Rdata")
      stopifnot(file.exists(regionList_file))
      all_regs <- get(load(regionList_file))
      
      
      geneList_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
      stopifnot(file.exists(geneList_file))
      pipeline_geneList <- get(load(geneList_file))
      
      stopifnot(pipeline_geneList %in% g2t_dt$entrezID)
      
      exprds_g2t_dt <- g2t_dt[g2t_dt$entrezID %in% pipeline_geneList,]
      
      
      comb_empPval_file <- file.path(pipFolder, hicds, exprds, script11same_name, "emp_pval_combined.Rdata" )
      stopifnot(file.exists(comb_empPval_file))
      comb_empPval <- eval(parse(text = load(paste0(comb_empPval_file))))
      stopifnot(setequal(all_regs, names(comb_empPval)))
      comb_empPval <- comb_empPval[all_regs]
      # ADJUST THE PVAL
      tad_adjCombPval <- p.adjust(comb_empPval, method="BH")
      stopifnot(names(tad_adjCombPval) == all_regs)
      
      
      ### retrieve limma signif genes
      topTable_DT_file <- file.path(pipFolder,  hicds, exprds, script1_name, "DE_topTable.Rdata")
      stopifnot(file.exists(topTable_DT_file))
      topTable_DT <- get(load(topTable_DT_file))
      topTable_DT$genes <- as.character(topTable_DT$genes)
      stopifnot(names(pipeline_geneList) %in% topTable_DT$genes)
      topTable_DT <- topTable_DT[topTable_DT$genes %in% names(pipeline_geneList),]
      topTable_DT$entrezID <- pipeline_geneList[paste0(topTable_DT$genes)]
      stopifnot(!is.na(topTable_DT$entrezID))
      stopifnot(!duplicated(topTable_DT$entrezID))
      stopifnot(topTable_DT$entrezID %in% pipeline_geneList)
      topTable_DT <- topTable_DT[topTable_DT$genes %in% names(pipeline_geneList),]
      topTable_DT$entrezID <- pipeline_geneList[paste0(topTable_DT$genes)]
      stopifnot(!is.na(topTable_DT$entrezID))
      stopifnot(!duplicated(topTable_DT$entrezID))
      stopifnot(exprds_g2t_dt$entrezID %in% topTable_DT$entrezID)
      
      
      tad_dt <- data.frame(region=names(tad_adjCombPval), tad_adjCombPval = as.numeric(tad_adjCombPval), stringsAsFactors=FALSE)
      
      tad_gene_dt <- merge(exprds_g2t_dt[,c("entrezID", "region")], tad_dt, by="region", all.x=TRUE, all.y=FALSE)
      
      out_dt <- merge(tad_gene_dt, topTable_DT[,c("entrezID", "logFC", "adj.P.Val")], by="entrezID", all.x=TRUE, all.y=FALSE )
      out_dt <- unique(out_dt)
      stopifnot(!duplicated(out_dt$entrezID))
      stopifnot(exprds_g2t_dt$entrezID %in% out_dt$entrezID)
      
      out_dt$gene_rank <- rank(out_dt$adj.P.Val, ties=tieMeth)
      out_dt$tad_rank <- rank(out_dt$tad_adjCombPval, ties=tieMeth)
      
      out_dt_cols <- colnames(out_dt)
      out_dt$hicds <- hicds
      out_dt$exprds <- exprds
      
      out_dt[, c("hicds", "exprds", out_dt_cols)]
      
    } # end-foreach iterating over exprds
    exprds_dt
  } # end-foreach iterating over hicds
  
  outFile <- file.path(outFolder, paste0("all_gene_tad_signif_dt.Rdata"))
  save(all_gene_tad_signif_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  # outFile <- "GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata"
  outFile <- file.path(outFolder, paste0( "all_gene_tad_signif_dt.Rdata"))
  cat("... load data\n")
  all_gene_tad_signif_dt <- get(load(outFile))
}


all_gene_tad_signif_dt$gene_signif <- all_gene_tad_signif_dt[,paste0(gene_signif_col)] <= gene_pval_thresh
all_gene_tad_signif_dt$cmpType <- all_cmps[all_gene_tad_signif_dt$exprds]  
stopifnot(!is.na(all_gene_tad_signif_dt$cmpType))

gene_plot_list <- list()

data_cmpType = "norm_vs_tumor"

all_signif_dt <- all_gene_tad_signif_dt[all_gene_tad_signif_dt$gene_signif,]

for(data_cmpType in c("norm_vs_tumor", "subtypes", "wt_vs_mut", "")) {
  
  
  if(data_cmpType != "") {
    signif_dt <- all_signif_dt[all_signif_dt$cmpType == data_cmpType,]
    nDS <- length(unique(paste0(all_gene_tad_signif_dt$hicds[all_gene_tad_signif_dt$cmpType==data_cmpType],
                                all_gene_tad_signif_dt$exprds[all_gene_tad_signif_dt$cmpType==data_cmpType])))
    
    
    
  } else {
    signif_dt <- all_signif_dt  
    nDS <- length(unique(paste0(all_gene_tad_signif_dt$hicds,all_gene_tad_signif_dt$exprds)))  
  }
  stopifnot(nDS > 0)
  
  conserved_signif_genes <- setNames(as.numeric(table(signif_dt$entrezID)),
                                     as.character(names(table(signif_dt$entrezID))))
  
  conserved_signif_genes <- conserved_signif_genes[conserved_signif_genes > 1]
  
  countConserv <- abs(sort(-(conserved_signif_genes)))
  
  countConserv_dt <- data.frame(countConserv)
  countConserv_dt$gene <- factor(rownames(countConserv_dt), levels=rownames(countConserv_dt))
  countConserv_dt$conservRatio <- countConserv_dt$countConserv/nDS
  stopifnot(countConserv_dt$conservRatio <= 1)
  
  
  
  gene_plot_list[[paste0(ifelse(data_cmpType=="", "all", data_cmpType))]] <-  ggbarplot(countConserv_dt, 
                                                                                        title=paste0("Conserv. signif. genes - ", data_cmpType),
                                                                                        subtitle=paste0("# DS = ", nDS, "; n = ", nrow(countConserv_dt)),
                                                                                        x = "gene", 
                                                                                        y = "conservRatio",
                                                                                        xlab="genes (sorted by conserv. ratio)",
                                                                                        ylab="conserv. ratio",
                                                                                        col="darkgrey", fill="darkgrey") +
    theme(axis.text.x=element_text(size=14),
          axis.text.y=element_text(size=14),
          plot.title=element_text(size=16, face = "bold"),
          axis.ticks.x = element_blank()) + 
    coord_cartesian(expand = F)
  
}

outFile <- file.path(outFolder, paste0( "gene_plot_list.Rdata"))
save(gene_plot_list, file=outFile, version=2)
cat("... load data\n")
gene_plot_list <- get(load(outFile))

outFile <- file.path(outFolder, paste0( "tad_plot_list.Rdata"))
cat("... load data\n")
save(tad_plot_list, file=outFile, version=2)
tad_plot_list <- get(load(outFile))



arrAll <- ggarrange(
  ggarrange(
    gene_plot_list[["all"]],
    tad_plot_list[["all"]],
    nrow=2
  ),
  ggarrange(
    gene_plot_list[["norm_vs_tumor"]],
    tad_plot_list[["norm_vs_tumor"]],
    nrow=2
  ),
  ggarrange(
    gene_plot_list[["subtypes"]],
    tad_plot_list[["subtypes"]],
    nrow=2
  ),
  ggarrange(
    gene_plot_list[["wt_vs_mut"]],
    tad_plot_list[["wt_vs_mut"]],
    nrow=2
  ),
  ncol = 4 
)


outFile <- file.path(outFolder, paste0("conserved_genes_regions_all_barplot.", plotType))
ggsave(plot = arrAll, filename = outFile, height=myHeightGG*2, width = myWidthGG*2)
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))

