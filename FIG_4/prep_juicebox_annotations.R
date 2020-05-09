# output format for juicebox should be:

# #chr1  x1	 x2	    chr2   y1	      y2	  name	score  strand1  strand2	color [optional fields]
# chr5   85000000  89000000   chr5   85000000   89000000    .     .      .        .       255,0,0


########################################################################################################################################################################################

# => v2: it is more correct to take the fpkmDT.Rdata !!!
# discussion with Marco 10.03.2020:
# should be normalized sample-wise to sum up to 1 !!!

startTime <- Sys.time()
cat(paste0("> Rscript prep_juicebox_annotations.R\n"))

script_name <- "prep_juicebox_annotations.R"

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))


# Rscript prep_juicebox_annotations.R <hicds> <features_to_plot>

# Rscript prep_juicebox_annotations.R K562_40kb GIMAP1 GIMAP2 GIMAP4 GIMAP5 GIMAP6 GIMAP7

plotType <- "foo"

source("../settings.R")

outFolder <- "PREP_JUICEBOX_ANNOTATIONS"
dir.create(outFolder, recursive = TRUE)

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)
gff_dt <- gff_dt[order(gff_dt$chromo, gff_dt$start, gff_dt$end, gff_dt$entrezID),]


args <- commandArgs(trailingOnly = TRUE)

tads_to_plot <- args[grepl("_TAD", args)]
hicds <- args[args %in% all_obs_hicds]
symbol_to_plot <- args[args %in% gff_dt$symbol]
entrez_to_plot <- args[args %in% gff_dt$entrezID]

if(length(symbol_to_plot) > 0) {
  stopifnot(symbol_to_plot %in% gff_dt$symbol)
}

if(length(entrez_to_plot) > 0) {
  stopifnot(entrez_to_plot %in% gff_dt$entrezID)
  newsymb <- gff_dt$symbol[gff_dt$entrezID == entrez_to_plot]
  symbol_to_plot <- c(symbol_to_plot, entrez_to_plot)
}

symbol_to_plot <- unique(symbol_to_plot)
entrez_to_plot <- unique(entrez_to_plot)
tads_to_plot <- unique(tads_to_plot)

tad_color <- "255,0,0"
gene_color <- "0,255,0"


plot_dt_2d <- data.frame(
  chr1=character(),
  x1=numeric(),
  x2=numeric(),
  chr2=character(),
  y1=numeric(),
  y2=numeric(),
  name=character(),
  score=character(),
  strand1=character(),
  strand2=character(),
  color=character(),
  stringsAsFactors = FALSE
)

plot_dt_1d <- data.frame(
  chr=character(),
  start=numeric(),
  end=numeric(),
  color=character(),
  stringsAsFactors = FALSE
)

if(length(tads_to_plot) > 0) {
  for(p_tad in tads_to_plot) {
    
    stopifnot(!is.na(hicds))
    
    tad_file <- file.path(runFolder, hicds, "genes2tad", "all_assigned_regions.txt")
    stopifnot(file.exists(tad_file))
    
    tad_dt <- read.delim(tad_file, header=FALSE, stringsAsFactors = FALSE, col.names = c("chromo", "region", "start", "end"))
    stopifnot(p_tad %in% tad_dt$region)
    
    tad_start <- tad_dt$start[tad_dt$region == p_tad]
    stopifnot(!is.na(tad_start))
    stopifnot(length(tad_start) == 1)
    tad_start <- tad_start - 1
    
    tad_end <- tad_dt$end[tad_dt$region == p_tad]
    stopifnot(!is.na(tad_end))
    stopifnot(length(tad_end) == 1)
    
    tad_chr <- tad_dt$chromo[tad_dt$region == p_tad]
    stopifnot(!is.na(tad_chr))
    stopifnot(length(tad_chr) == 1)
    
    chromo <- gsub("(chr.+?)_.+", "\\1", p_tad)
    
    stopifnot(chromo == tad_chr)
    
    # #chr1  x1	 x2	    chr2   y1	      y2	  name	score  strand1  strand2	color [optional fields]
    # chr5   85000000  89000000   chr5   85000000   89000000    .     .      .        .       255,0,0
    new_plot_dt_2d <- data.frame(
      chr1=chromo,
      x1=tad_start,
      x2=tad_end,
      chr2=chromo,
      y1=tad_start,
      y2=tad_end,
      name=p_tad,
      score= ".",
      strand1=".",
      strand2=".",
      color=tad_color,
      stringsAsFactors = FALSE
    )
    
    plot_dt_2d <- rbind(plot_dt_2d, new_plot_dt_2d)
    
  }
  
  outFile <- file.path(outFolder, paste0(hicds, "_", paste0(tads_to_plot, collapse="_"), "_juicebox_2d_annot.txt"))
  write.table(plot_dt_2d, col.names = TRUE, row.names = FALSE, sep="\t", quote=FALSE, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
}



if(length(symbol_to_plot) > 0) {
  for(p_symb in symbol_to_plot) {
    
    
    symb_start <- gff_dt$start[gff_dt$symbol == p_symb]
    stopifnot(!is.na(symb_start))
    stopifnot(length(symb_start) == 1)
    symb_start <- symb_start - 1
    
    symb_end <- gff_dt$end[gff_dt$symbol == p_symb]
    stopifnot(!is.na(symb_end))
    stopifnot(length(symb_end) == 1)
    
    symb_chr <- gff_dt$chromo[gff_dt$symbol == p_symb]
    stopifnot(!is.na(symb_chr))
    stopifnot(length(symb_chr) == 1)
    
    # BED format
    # #chr start end color
    # chr5   85000000  89000000 0,255,0
    new_plot_dt_1d <- data.frame(
      chr=symb_chr,
      start=symb_start,
      end=symb_end,
      color=gene_color,
      stringsAsFactors = FALSE
    )
    
    plot_dt_1d <- rbind(plot_dt_1d, new_plot_dt_1d)
    
  }
  
  outFile <- file.path(outFolder, paste0(hicds, "_", paste0(symbol_to_plot, collapse="_"), "_juicebox_1d_annot.txt"))
  write.table(plot_dt_1d, col.names = FALSE, row.names = FALSE, sep="\t", quote=FALSE, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
}



