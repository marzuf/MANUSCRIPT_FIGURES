# output format for juicebox should be:

#chr1  x1	 x2	    chr2   y1	      y2	  name	score  strand1  strand2	color [optional fields]
chr5   85000000  89000000   chr5   85000000   89000000    .     .      .        .       255,0,0


########################################################################################################################################################################################

# => v2: it is more correct to take the fpkmDT.Rdata !!!
# discussion with Marco 10.03.2020:
# should be normalized sample-wise to sum up to 1 !!!

startTime <- Sys.time()
cat(paste0("> Rscript look_geneTAD_expression_withRank_v2_noFilter.R\n"))

script_name <- "look_geneTAD_expression_withRank_v2_noFilter.R"

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

# Rscript look_geneTAD_expression_withRank_v2.R ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutCTNNB1 IFIT3
# Rscript look_geneTAD_expression_withRank_v2.R ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_wt_mutBRAF HLA-DQA1


setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
gff_dt <- gff_dt[order(gff_dt$chromo, gff_dt$start, gff_dt$end, gff_dt$entrezID),]

hicds="ENCSR489OCU_NCI-H460_40kb"
exprds="TCGAlusc_norm_lusc"
gene_to_plot="IFIT3"

col1 <- pal_futurama()(5)[1]
col2 <- pal_futurama()(5)[5]
col1 <- pal_aaas()(5)[4]
col2 <- pal_npg()(5)[5]

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
hicds <- args[1]
exprds <- args[2]
gene_to_plot <- as.character(args[3])

entrez_to_plot <- NA

if(!gene_to_plot %in% gff_dt$entrezID) stopifnot(gene_to_plot %in% gff_dt$symbol)

if(gene_to_plot %in% gff_dt$symbol) {
  entrez_to_plot <- gff_dt$entrezID[gff_dt$symbol == gene_to_plot]
}
if(gene_to_plot %in% gff_dt$entrezID) {
  entrez_to_plot <- gene_to_plot
  gene_to_plot <- gff_dt$symbol[gff_dt$entrezID == entrez_to_plot]
}
stopifnot(!is.na(entrez_to_plot))
stopifnot( length(entrez_to_plot) == 1)

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
inDT <- get(load(file.path(runFolder, "GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata")))
inDT <- inDT[inDT$hicds == hicds & inDT$exprds == exprds,]

stopifnot(entrez_to_plot %in% inDT$entrezID)
tad_to_plot <- inDT$region[inDT$entrezID == entrez_to_plot]
stopifnot(!is.na(tad_to_plot))
stopifnot( length(tad_to_plot) == 1)


tad_plot_rank <- unique(inDT$tad_rank[inDT$hicds == hicds & inDT$exprds == exprds & inDT$region == tad_to_plot])
stopifnot(!is.na(tad_plot_rank))
stopifnot(length(tad_plot_rank) == 1)

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

outFolder <- file.path("LOOK_GENETAD_EXPRESSION_WITH_RANK_V2")
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
