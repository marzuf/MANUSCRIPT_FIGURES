########################################################################################################################################################################################

# => v2: it is more correct to take the fpkmDT.Rdata !!!
# discussion with Marco 10.03.2020:
# should be normalized sample-wise to sum up to 1 !!!

startTime <- Sys.time()
cat(paste0("> Rscript look_TAD_expression_withRank_v2_withSignif.R\n"))

script_name <- "look_TAD_expression_withRank_v2_withSignif.R"

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

# Rscript look_TAD_expression_withRank_v2_withSignif.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr10_TAD16
# Rscript look_TAD_expression_withRank_v2_withSignif.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr17_TAD162
# Rscript look_TAD_expression_withRank_v2_withSignif.R ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker chr10_TAD16
# Rscript look_TAD_expression_withRank_v2_withSignif.R ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker chr17_TAD162
# Rscript look_TAD_expression_withRank_v2_withSignif.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr11_TAD390
# Rscript look_TAD_expression_withRank_v2_withSignif.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr10_TAD268
# Rscript look_TAD_expression_withRank_v2_withSignif.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc chr11_TAD390
# Rscript look_TAD_expression_withRank_v2_withSignif.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc chr10_TAD268
# 
# FOR CONSERVED REGION 130 - top1 all (GIMAP)
# Rscript look_TAD_expression_withRank_v2_withSignif.R ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf chr7_TAD568
# FOR CONSERVED REGION 42 - top 2 all (AKR1C)
# Rscript look_TAD_expression_withRank_v2_withSignif.R LG2_40kb TCGAluad_mutKRAS_mutEGFR chr10_TAD15
# FOR CONSERVED REGION 11 - top3 all (CD8)
# Rscript look_TAD_expression_withRank_v2_withSignif.R Panc1_rep12_40kb TCGApaad_wt_mutKRAS chr1_TAD544


# FOR CONSERVED REGION 22 - top1 norm_tumor GIMAP
# Rscript look_TAD_expression_withRank_v2_withSignif.R ENCSR444WCZ_A549_40kb TCGAlusc_norm_lusc chr7_TAD553
# FOR CONSERVED REGION 9 - top2 norm_tumor M1T
# Rscript look_TAD_expression_withRank_v2_withSignif.R GSE118514_RWPE1_40kb TCGAprad_norm_prad chr16_TAD163
# FOR CONSERVED REGION 10 - top3 norm_tumor KRT
# Rscript look_TAD_expression_withRank_v2_withSignif.R ENCSR346DCU_LNCaP_40kb TCGAprad_norm_prad chr17_TAD146

# TAD with highest number of missed low fc same dir
# Rscript look_TAD_expression_withRank_v2_withSignif.R LG1_40kb TCGAluad_nonsmoker_smoker chr19_TAD43
# example of TAD with recurrently missed
# Rscript look_TAD_expression_withRank_v2_withSignif.R LG1_40kb TCGAluad_nonsmoker_smoker chr10_TAD313

# example for 1 dataset
#missedLowFCdiffDir: "chr4_TAD635"
#missedLowFCsameDir:  "chr10_TAD313"
#highFCsameDir  & ex_dt$tadSignifOnly ][1]  "chr10_TAD16"
## ex_dt$region[ex_dt$highFCdiffDir & ex_dt$tadSignifOnly][1]
# chr1_TAD566  ENCSR079VIJ_G401_40kb/TCGAkich_norm_kich
#ex_dt[ex_dt$region == "chr4_TAD635", c("ds_id", "geneSignif", "tadSignif", "symbol")]
#ex_dt[ex_dt$region == "chr10_TAD313", c("ds_id", "geneSignif", "tadSignif", "symbol")]
#ex_dt[ex_dt$region == "chr10_TAD16", c("ds_id", "geneSignif", "tadSignif", "symbol")]


#Rscript look_TAD_expression_withRank_v2_withSignif.R LG1_40kb TCGAluad_nonsmoker_smoker chr19_TAD43
#Rscript look_TAD_expression_withRank_v2_withSignif.R LG1_40kb TCGAluad_nonsmoker_smoker chr10_TAD313
#Rscript look_TAD_expression_withRank_v2_withSignif.R LG1_40kb TCGAluad_nonsmoker_smoker chr4_TAD635
#Rscript look_TAD_expression_withRank_v2_withSignif.R LG1_40kb TCGAluad_nonsmoker_smoker chr10_TAD16
#Rscript look_TAD_expression_withRank_v2_withSignif.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich chr1_TAD566



hicds="ENCSR489OCU_NCI-H460_40kb"
exprds="TCGAlusc_norm_lusc"
tad_to_plot="chr11_TAD390"

col1 <- pal_futurama()(5)[1]
col2 <- pal_futurama()(5)[5]
col1 <- pal_aaas()(5)[4]
col2 <- pal_npg()(5)[5]

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
hicds <- args[1]
exprds <- args[2]
tad_to_plot <- args[3]


my_xlab <- "TAD genes (ordered by start positions)"
my_ylab <- "RNA-seq TPM [log10]"

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4
myHeightGG <- 6
myWidthGG <- 7.5

source("../settings.R")

cat("load inDT \n")
inDT <- get(load(file.path(runFolder, "/GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata")))
inDT <- inDT[inDT$hicds == hicds & inDT$exprds == exprds,]

tad_plot_rank <- unique(inDT$tad_rank[inDT$hicds == hicds & inDT$exprds == exprds & inDT$region == tad_to_plot])
stopifnot(!is.na(tad_plot_rank))
stopifnot(length(tad_plot_rank) == 1)

plotSub <- paste0(tad_to_plot, " - rank: ", tad_plot_rank)


stopifnot(hicds %in% names(hicds_names))
hicds_lab <- hicds_names[paste0(hicds)]

stopifnot(exprds %in% names(exprds_names))
exprds_lab <- exprds_names[paste0(exprds)]

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

mainFolder <- file.path(runFolder)
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
settingFolder <- file.path(mainFolder, "PIPELINE", "INPUT_FILES")

outFolder <- file.path("LOOK_TAD_EXPRESSION_WITH_RANK_V2_WITHSIGNIF")
dir.create(outFolder, recursive = TRUE)

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

# count_file <- file.path(pipFolder, hicds, exprds, "1_runGeneDE", "DE_rnaseqDT.Rdata")
count_file <- file.path(pipFolder, hicds, exprds, "0_prepGeneData", "rna_fpkmDT.Rdata")  # corrected here 11.03.2020
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

stopifnot(cond1 %in% names(cond1_names))
cond1_lab <- cond1_names[paste0(cond1)]
stopifnot(cond2 %in% names(cond2_names))
cond2_lab <- cond2_names[paste0(cond2)]
stopifnot(!is.na(cond1_lab))
stopifnot(!is.na(cond2_lab))

stopifnot(!is.na(toplot_dt))
toplot_dt$cond <- ifelse(toplot_dt$variable %in% samp1, cond1, 
                         ifelse(toplot_dt$variable %in% samp2, cond2, NA ))
#toplot_dt$cond <- ifelse(toplot_dt$variable %in% samp1, cond1_lab, 
#                         ifelse(toplot_dt$variable %in% samp2, cond2_lab, NA ))

check_dt <- toplot_dt[,c("variable", "cond")]
cat("table(check_dt$cond)\n")
table(check_dt$cond)

cat(paste0(length(samp1), "\n"))
cat(paste0(length(samp2), "\n"))


stopifnot(!is.na(toplot_dt$cond))
toplot_dt$symbol <- entrez2symb[paste0(toplot_dt$entrezID)]
stopifnot(!is.na(toplot_dt$symbol))

toplot_dt <- toplot_dt[order(toplot_dt$chromo, toplot_dt$start, toplot_dt$end, toplot_dt$value),]
toplot_dt$value_log10 <- log10(toplot_dt$value + log10_offset)

# withRank_toplot_dt <- do.call(rbind, by(toplot_dt, list(toplot_dt$symbol, toplot_dt$cond), function(x) {
#   dt <- x[order(x$value, decreasing = TRUE),]
#   dt$samp_rank <- 1:nrow(dt)
#   dt
# }))

all_conds <- c(cond1, cond2)

withRank_toplot_dt2 <- do.call(rbind, by(toplot_dt, list(toplot_dt$symbol), function(x) {
  x$cond <- factor(x$cond, levels=all_conds)
  dt <- x[order(-as.numeric(x$cond), x$value, decreasing = TRUE),]
  dt$samp_rank <- 1:nrow(dt)
  dt
}))
withRank_toplot_dt2$hicds <- hicds
withRank_toplot_dt2$exprds <- exprds

cat("merge withRank and inDT \n")

withRank_toplot_dt2 <- merge(withRank_toplot_dt2, inDT, all.x=TRUE, all.y=FALSE, by=intersect(colnames(inDT), colnames(withRank_toplot_dt2)))
withRank_toplot_dt2$pval_lab <- ifelse(withRank_toplot_dt2$adj.P.Val <= 0.001, "***", 
                       ifelse(withRank_toplot_dt2$adj.P.Val <= 0.01, "**",
                              ifelse(withRank_toplot_dt2$adj.P.Val <= 0.05, "*", "")))


tmp <- withRank_toplot_dt2[,c("symbol", "start", "end", "gene_rank", "pval_lab")]
stopifnot(!is.na(tmp$pval_lab))
tmp <- unique(tmp)
tmp <- tmp[order(tmp$start, tmp$end),]
tmp$lab <- paste0(tmp$symbol, "\n(rank: ", tmp$gene_rank, ")\n", tmp$pval_lab)

withRank_toplot_dt2$symbol <- factor(withRank_toplot_dt2$symbol, levels=tmp$symbol)
withRank_toplot_dt2$cond <- factor(withRank_toplot_dt2$cond, levels = c(cond1,cond2))

withRank_toplot_dt2$symbol_lab <- paste0(withRank_toplot_dt2$symbol, "\n(rank: ", withRank_toplot_dt2$gene_rank, ")\n", withRank_toplot_dt2$pval_lab )
withRank_toplot_dt2$symbol_lab <- factor(withRank_toplot_dt2$symbol_lab, levels=tmp$lab)
stopifnot(!is.na(withRank_toplot_dt2$symbol_lab))

save(withRank_toplot_dt2, file ="withRank_toplot_dt2.Rdata")

subTit <- paste0(tad_to_plot, " (rank: ", tad_plot_rank, ")")


cat("table(withRank_toplot_dt2$cond)\n")
table(withRank_toplot_dt2$cond)


check_dt <- withRank_toplot_dt2[,c("cond", "variable")]
check_dt <- unique(check_dt)
stopifnot(sum(table(check_dt$cond) [grepl(paste0("^", cond1, "$"), names(table(check_dt$cond)))]) == length(samp1))

stopifnot(sum(table(check_dt$cond) [grepl(paste0("^", cond1, "$"), names(table(check_dt$cond)))]) == length(samp1))
stopifnot(sum(table(check_dt$cond) [grepl(paste0("^", cond2, "$"), names(table(check_dt$cond)))]) == length(samp2))


check_dt$cond_labels <- ifelse(as.character(check_dt$cond) == cond1, cond1_lab, 
                                    ifelse(as.character(check_dt$cond) == cond2, cond2_lab, NA))
stopifnot(!is.na(check_dt$cond_labels))

all_conds <- c(cond1_lab, cond2_lab)

cond_labels <- paste0(all_conds, " (n=" , table(check_dt$cond_labels) [all_conds], ")")


p_var_boxplot <-  ggplot(withRank_toplot_dt2, aes(x = symbol_lab, y = value_log10, color = cond)) + 
  # geom_boxplot(notch = TRUE, outlier.shape=NA)+
  # geom_jitter(aes(colour = cond), position=position_jitterdodge())+
  geom_boxplot(notch = TRUE, outlier.shape=NA)+
  geom_point(aes(fill=cond), position=position_jitterdodge(),  alpha=0.5) +
  
  
 
  ggtitle(paste0(hicds_lab, " - ", exprds_lab), subtitle = paste0(subTit))+
  scale_x_discrete(name=my_xlab)+
  scale_y_continuous(name=paste0(my_ylab),
                     breaks = scales::pretty_breaks(n = 20))+
  
  scale_color_manual(values=c(col1, col2), labels=cond_labels)+
  scale_fill_manual(values=c(col1, col2), labels=cond_labels)+
  
  labs(fill  = paste0("Cond."), color=paste0("Cond.")) +
  theme( 
	text = element_text(family=fontFamily),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "grey"),
    panel.grid.minor.y = element_line(colour = "grey"),
    axis.line.x= element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .2, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12),
    axis.text.x =element_text(color="black", hjust=0.5,vjust = 0.5, size=12, face="bold"),
    # axis.ticks.x = element_blank(),
    axis.title.y = element_text(color="black", size=13),
    axis.title.x = element_text(color="black", size=13),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.text = element_text(size=12),
    legend.key = element_blank(),
    legend.title = element_text(face="bold")
  )

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", tad_to_plot, "_allSamples_exprValues_boxplot.", plotType))
ggsave(plot = p_var_boxplot, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))





##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))

