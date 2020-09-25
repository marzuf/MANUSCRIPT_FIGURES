

########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript look_TAD_expression_withRank_withMutStatus_v4.R\n"))

script_name <- "look_TAD_expression_withRank_withMutStatus_v4.R"

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
require(ggsci)
ggsci_pal <- "d3"
ggsci_subpal <- ""
require(reshape2)
log10_offset <- 0.01


require(ggsci)
require(ggpubr)
statTest <- "wilcox.test"


# Rscript look_TAD_expression_withRank_withMutStatus_v4.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr10_TAD16

hicds="ENCSR489OCU_NCI-H460_40kb"
exprds="TCGAluad_mutKRAS_mutEGFR"
tad_to_plot="chr10_TAD16"

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
hicds <- args[1]
exprds <- args[2]
tad_to_plot <- args[3]

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4


source("../settings.R")

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))


if(tad_to_plot == "chr17_TAD162") myWidthGG <- myWidthGG * 1.1

boxSpacing <- 0.2
mutSpacing <- 0.075

my_pal <- setNames(c(
  pal_lancet()(1)[1],
  pal_lancet()(4)[4],
  pal_lancet()(3)[3]),c("mutKRAS/withMut", "mutKRAS/noMut", "mutEGFR/noMut"))


myHeightGG <- 5.5
myWidthGG <- 7


######################################################################################################

cat("load inDT \n")
inDT <- get(load(file.path(runFolder, "GENE_RANK_TAD_RANK", "all_gene_tad_signif_dt.Rdata")))
inDT <- inDT[inDT$hicds == hicds & inDT$exprds == exprds,]

tad_plot_rank <- unique(inDT$tad_rank[inDT$hicds == hicds & inDT$exprds == exprds & inDT$region == tad_to_plot])
stopifnot(!is.na(tad_plot_rank))
stopifnot(length(tad_plot_rank) == 1)

stopifnot(hicds %in% names(hicds_names))
hicds_lab <- hicds_names[paste0(hicds)]

stopifnot(exprds %in% names(exprds_names))
exprds_lab <- exprds_names[paste0(exprds)]

plotSub <- paste0(tad_to_plot, " - rank: ", tad_plot_rank)

my_xlab <- "TAD genes (ordered by start positions)"
#my_ylab <- "RNA-seq expression count [log10]"
my_ylab <- "RNA-seq TPM [log10]" # corrected 15.03.2020


mainFolder <- runFolder
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path(mainFolder, "PIPELINE", "INPUT_FILES")

outFolder <- file.path("LOOK_TAD_EXPRESSION_WITHRANK_WITHMUTSTATUS_V4")
dir.create(outFolder, recursive = TRUE)

mutSamples <- get(load(file.path(runFolder, "NFE2L2_KEAP1_MUTSAMPLES", "mut_samples.Rdata")))

entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
gff_dt <- gff_dt[order(gff_dt$chromo, gff_dt$start, gff_dt$end, gff_dt$entrezID),]

g2t_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
g2t_dt <- read.delim(g2t_file, col.names=c("entrezID", "chromo", "start", "end", "region"), header=FALSE, stringsAsFactors = FALSE)
g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
stopifnot(!duplicated(g2t_dt$entrezID))

stopifnot(tad_to_plot %in% g2t_dt$region)

settingFile <- file.path(settingFolder, hicds, paste0("run_settings_", exprds, ".R"))
stopifnot(file.exists(settingFile))
source(settingFile)

samp1 <- get(load(file.path(setDir, sample1_file)))
samp2 <- get(load(file.path(setDir, sample2_file)))

geneList <- get(load(file.path(pipFolder, hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")))
stopifnot(geneList %in% g2t_dt$entrezID)

toplot_dt <- g2t_dt[g2t_dt$region == tad_to_plot & g2t_dt$entrezID %in% geneList,]
toplot_dt <- toplot_dt[order(toplot_dt$start, toplot_dt$end),]

stopifnot(toplot_dt$entrezID %in% geneList)
plot_geneList <- geneList[geneList %in% toplot_dt$entrezID]

count_file <- file.path(pipFolder, hicds, exprds, "0_prepGeneData", "rna_fpkmDT.Rdata")  # corrected here 11.03.2020 discussion to marco which data to plot
stopifnot(file.exists(count_file))
fpkm_dt <- get(load(count_file))
# changed here 11.03.2020 -> should be normalized sample-wise
fpkm_dt2 <- apply(fpkm_dt, 2, function(x)x/sum(x))
# stopifnot(colSums(fpkm_dt2) == 1)
stopifnot(abs(colSums(fpkm_dt2) - 1) <= 10^-4)
# and then multiply by 10^6 to have FPKM
fpkm_dt2 <- fpkm_dt2*10^6
fpkm_dt2 <- data.frame(fpkm_dt2, check.names = FALSE)
stopifnot(dim(fpkm_dt) == dim(fpkm_dt2))
stopifnot(rownames(fpkm_dt) == rownames(fpkm_dt2))
stopifnot(colnames(fpkm_dt) == colnames(fpkm_dt2))

fpkm_dt <- fpkm_dt2
stopifnot(names(geneList) %in% rownames(fpkm_dt))

fpkm_plot_dt <- fpkm_dt[rownames(fpkm_dt) %in% names(plot_geneList),]
fpkm_plot_dt$entrezID <- plot_geneList[paste0(rownames(fpkm_plot_dt))]
stopifnot(!is.na(fpkm_plot_dt$entrez))

stopifnot(samp1 %in% colnames(fpkm_plot_dt))
stopifnot(samp2 %in% colnames(fpkm_plot_dt))

fpkm_plot_dt <- fpkm_plot_dt[,c("entrezID", samp1, samp2)]
m_fpkm_dt <- melt(fpkm_plot_dt, id="entrezID")

stopifnot(m_fpkm_dt$entrezID %in% toplot_dt$entrezID)

toplot_dt <- merge(m_fpkm_dt, toplot_dt, merge="entrezID", all.x=TRUE, all.y=FALSE)
stopifnot(!is.na(toplot_dt))
toplot_dt$cond <- ifelse(toplot_dt$variable %in% samp1, cond1, 
                         ifelse(toplot_dt$variable %in% samp2, cond2, NA ))
stopifnot(!is.na(toplot_dt$cond))
toplot_dt$symbol <- entrez2symb[paste0(toplot_dt$entrezID)]
stopifnot(!is.na(toplot_dt$symbol))

toplot_dt <- toplot_dt[order(toplot_dt$chromo, toplot_dt$start, toplot_dt$end, toplot_dt$value),]
toplot_dt$value_log10 <- log10(toplot_dt$value + log10_offset)

withRank_toplot_dt2 <- do.call(rbind, by(toplot_dt, list(toplot_dt$symbol), function(x) {
  x$cond <- factor(x$cond, levels=c(cond1,cond2))
  dt <- x[order(-as.numeric(x$cond), x$value, decreasing = TRUE),]
  dt$samp_rank <- 1:nrow(dt)
  dt
}))
withRank_toplot_dt2$hicds <- hicds
withRank_toplot_dt2$exprds <- exprds

cat("merge withRank and inDT \n")

withRank_toplot_dt2 <- merge(withRank_toplot_dt2, inDT, all.x=TRUE, all.y=FALSE, by=intersect(colnames(inDT), colnames(withRank_toplot_dt2)))


tmp <- withRank_toplot_dt2[,c("symbol", "start", "end", "gene_rank")]
tmp <- unique(tmp)
tmp <- tmp[order(tmp$start, tmp$end),]
tmp$lab <- paste0(tmp$symbol, "\n(rank: ", tmp$gene_rank, ")")

withRank_toplot_dt2$symbol <- factor(withRank_toplot_dt2$symbol, levels=tmp$symbol)

withRank_toplot_dt2$cond_dots <- withRank_toplot_dt2$cond
withRank_toplot_dt2$cond_dots <- ifelse(as.character(withRank_toplot_dt2$variable) %in% mutSamples, "withMut", as.character(withRank_toplot_dt2$cond_dots))
withRank_toplot_dt2$cond_dots <- factor(withRank_toplot_dt2$cond_dots, levels = c(cond1,cond2, "withMut"))
stopifnot(!is.na(withRank_toplot_dt2$cond_dots))

withRank_toplot_dt2$cond_sh <- withRank_toplot_dt2$cond
withRank_toplot_dt2$cond_sh <- ifelse(as.character(withRank_toplot_dt2$variable) %in% mutSamples, "withMut", "noMut")
withRank_toplot_dt2$cond_sh <- factor(withRank_toplot_dt2$cond_sh, levels = c("noMut", "withMut"))
stopifnot(!is.na(withRank_toplot_dt2$cond_sh))

withRank_toplot_dt2$symbol_lab <- paste0(withRank_toplot_dt2$symbol, "\n(rank: ", withRank_toplot_dt2$gene_rank, ")")
withRank_toplot_dt2$symbol_lab <- factor(withRank_toplot_dt2$symbol_lab, levels=tmp$lab)
stopifnot(!is.na(withRank_toplot_dt2$symbol_lab))

mutcols <- setNames(c("#FC7715FF", "#AFBF41FF"), c("noMut", "withMut"))

withRank_toplot_dt2 <- withRank_toplot_dt2[order(withRank_toplot_dt2$cond_sh),]

withRank_toplot_dt2$symbol_pos <- as.numeric(withRank_toplot_dt2$symbol_lab)
withRank_toplot_dt2$cond_pos <- ifelse(as.numeric(withRank_toplot_dt2$cond) == 1, -boxSpacing, boxSpacing)
withRank_toplot_dt2$box_pos <- withRank_toplot_dt2$symbol_pos + withRank_toplot_dt2$cond_pos

withRank_toplot_dt2b <- do.call(rbind, by(withRank_toplot_dt2, withRank_toplot_dt2$box_pos, function(x) {
  if(length(unique(as.numeric(x$cond_sh))) == 1) {
    x$mut_pos <-0  # if there is only wt or mut -> put at the middle of the box
  } else {
    x$mut_pos <- ifelse(as.numeric(x$cond_sh) == 1, -mutSpacing, mutSpacing)
  }
  x
}))
withRank_toplot_dt2b$dot_pos <- withRank_toplot_dt2b$box_pos + withRank_toplot_dt2b$mut_pos

tad_geneExpr_withMut_dt <- withRank_toplot_dt2b
saveFile <- file.path(outFolder, paste0("fig3BwithMut_", hicds, "_", exprds, "_", tad_to_plot, "_tad_geneExpr_withMut_dt.Rdata"))
save(tad_geneExpr_withMut_dt, file=saveFile, version=2)
cat(paste0("... written:" , saveFile, "\n"))

plot_levels <- c("mutKRAS/withMut", "mutKRAS/noMut", "mutEGFR/noMut")  

tad_geneExpr_withMut_dt$plot_cond <- file.path(as.character(tad_geneExpr_withMut_dt$cond),
                                               as.character(tad_geneExpr_withMut_dt$cond_sh) )

tad_geneExpr_withMut_dt$plot_cond <- factor(tad_geneExpr_withMut_dt$plot_cond, levels=plot_levels)

# three boxplots per gene (egfr, kras-only, kras+nrf2) and compute a p-value for each comparison 
# (in particular egfr vs. kras-only)

subTit <- paste0(tad_to_plot, " (rank: ", tad_plot_rank, ")")
plotTit <- paste0(hicds_lab, " - ", exprds_lab)

my_cmps <- combn(levels(tad_geneExpr_withMut_dt$plot_cond), 2, simplify = FALSE)

myG_theme <-   theme( 
  text = element_text(family=fontFamily),
  plot.title = element_text(hjust = 0.5, face = "bold", size=16),
  plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
  # panel.grid = element_blank(),
  # panel.grid.major.y = element_line(colour = "grey"),
  # panel.grid.minor.y = element_line(colour = "grey"),
  axis.line.x= element_line(size = .2, color = "black"),
  axis.line.y = element_line(size = .2, color = "black"),
  axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12),
  # axis.text.x =element_text(color="black", hjust=0.5,vjust = 0.5, size=12, face="bold"),
  axis.text.x = element_blank(),
  # axis.ticks.x = element_blank(),
  axis.title.y = element_text(color="black", size=14),
  # axis.title.x = element_text(color="black", size=14),
  axis.title.x = element_blank(),
  # panel.border = element_blank(),
  # panel.background = element_rect(fill = "transparent"),
  legend.background =  element_rect(),
  legend.text = element_text(size=11),
  legend.key = element_blank(),
  legend.title = element_text(face="bold", size=12)
)

# Box plot facetted by "symbol"
p <- ggboxplot(tad_geneExpr_withMut_dt, 
               x = "plot_cond",
               y = "value_log10",
               color = "plot_cond",
               add = "jitter",
               facet.by = "symbol", 
               short.panel.labs = TRUE) +
  ggtitle(plotTit, subtitle = subTit)+
  guides(shape=F, fill=F)+
  scale_y_continuous(name=paste0(my_ylab),
                     breaks = scales::pretty_breaks(n = 10))+
  scale_color_manual(values=my_pal)+
  # labs(color="foo")+
  labs(color=paste0("(", statTest,")")) +
  myG_theme + 
  stat_compare_means(comparisons=my_cmps, 
                     method=statTest,
                     p.adjust.method = "none") #, p.adjust.method)

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", tad_to_plot, 
                                       "_allSamples_exprValues_boxplot_facet_withTest.", plotType))
ggsave(plot = p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


stopifnot(sum(grepl("mutKRAS", as.character(tad_geneExpr_withMut_dt$plot_cond))) ==
  sum(as.character(tad_geneExpr_withMut_dt$cond) == "mutKRAS"))
stopifnot(sum(grepl("mutEGFR", as.character(tad_geneExpr_withMut_dt$plot_cond))) ==
            sum(as.character(tad_geneExpr_withMut_dt$cond) == "mutEGFR"))



##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))

