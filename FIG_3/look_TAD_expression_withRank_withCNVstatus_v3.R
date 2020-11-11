#inkscape ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_chr11_TAD390_allSamples_exprValues_boxplot.svg --export-pdf=ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_chr11_TAD390_allSamples_exprValues_boxplot.pdf
#inkscape ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_chr10_TAD268_allSamples_exprValues_boxplot.svg --export-pdf=ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_chr10_TAD268_allSamples_exprValues_boxplot.pdf
#inkscape ENCSR489OCU_NCI-H460_40kb_TCGAluad_mutKRAS_mutEGFR_chr10_TAD16_allSamples_exprValues_boxplot.svg --export-pdf=ENCSR489OCU_NCI-H460_40kb_TCGAluad_mutKRAS_mutEGFR_chr10_TAD16_allSamples_exprValues_boxplot.pdf
#inkscape ENCSR489OCU_NCI-H460_40kb_TCGAluad_mutKRAS_mutEGFR_chr17_TAD162_allSamples_exprValues_boxplot.svg --export-pdf=ENCSR489OCU_NCI-H460_40kb_TCGAluad_mutKRAS_mutEGFR_chr17_TAD162_allSamples_exprValues_boxplot.pdf


########################################################################################################################################################################################
startTime <- Sys.time()
cat(paste0("> Rscript look_TAD_expression_withRank_withCNVstatus_v3.R\n"))

script_name <- "look_TAD_expression_withRank_withCNVstatus_v3.R"

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


# Rscript look_TAD_expression_withRank_withMutStatus_v3.R <hicds> <exprds> <gene_symbol>
# Rscript look_TAD_expression_withRank_withMutStatus_v3.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc chr11_TAD390 # MMPP
# Rscript look_TAD_expression_withRank_withMutStatus_v3.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc chr10_TAD268 # SFTPA
# Rscript look_TAD_expression_withRank_withMutStatus_v3.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr11_TAD390
# Rscript look_TAD_expression_withRank_withMutStatus_v3.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr10_TAD268

# ENCSR489OCU_NCI-H460_40kb	TCGAluad_mutKRAS_mutEGFR	chr10_TAD16	AKR1C1,AKR1C2,AKR1C3
# ENCSR489OCU_NCI-H460_40kb	TCGAluad_mutKRAS_mutEGFR	chr17_TAD162	HOXB2,HOXB3,HOXB4,HOXB5,HOXB6,HOXB7
# ENCSR489OCU_NCI-H460_40kb	TCGAluad_nonsmoker_smoker	chr10_TAD16	AKR1C1,AKR1C2,AKR1C3
# ENCSR489OCU_NCI-H460_40kb	TCGAluad_nonsmoker_smoker	chr17_TAD162	HOXB2,HOXB3,HOXB4,HOXB5,HOXB6,HOXB7
# 
# ENCSR489OCU_NCI-H460_40kb	TCGAluad_norm_luad	chr11_TAD390	MMP1,MMP12,MMP13
# ENCSR489OCU_NCI-H460_40kb	TCGAluad_norm_luad	chr10_TAD268	BEND3P3,SFTPA1,SFTPA2
# ENCSR489OCU_NCI-H460_40kb	TCGAlusc_norm_lusc	chr11_TAD390	MMP1,MMP10,MMP12,MMP13,MMP3
# ENCSR489OCU_NCI-H460_40kb	TCGAlusc_norm_lusc	chr10_TAD268	BEND3P3,SFTPA1,SFTPA2


# Rscript look_TAD_expression_withRank_withCNVstatus_v3.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr10_TAD16
# Rscript look_TAD_expression_withRank_withMutStatus_v3.R ENCSR489OCU_NCI-H460_40kb TCGAluad_mutKRAS_mutEGFR chr17_TAD162
# Rscript look_TAD_expression_withRank_withMutStatus_v3.R ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker chr10_TAD16
# Rscript look_TAD_expression_withRank_withMutStatus_v3.R ENCSR489OCU_NCI-H460_40kb TCGAluad_nonsmoker_smoker chr17_TAD162
# Rscript look_TAD_expression_withRank_withMutStatus_v3.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr11_TAD390
# Rscript look_TAD_expression_withRank_withMutStatus_v3.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad chr10_TAD268
# Rscript look_TAD_expression_withRank_withMutStatus_v3.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc chr11_TAD390
# Rscript look_TAD_expression_withRank_withMutStatus_v3.R ENCSR489OCU_NCI-H460_40kb TCGAlusc_norm_lusc chr10_TAD268
# 

hicds="ENCSR489OCU_NCI-H460_40kb"
exprds="TCGAluad_mutKRAS_mutEGFR"
tad_to_plot="chr10_TAD16"

hicds="ENCSR489OCU_NCI-H460_40kb"
exprds="TCGAlusc_norm_lusc"
tad_to_plot="chr11_TAD390"

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 3)
hicds <- args[1]
exprds <- args[2]
tad_to_plot <- args[3]

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 500, 7)
myWidth <- myHeight
plotCex <- 1.4
myHeightGG <- 6
myWidthGG <- 7.5


source("../settings.R")

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")
registerDoMC(ifelse(SSHFS, 2, 80))

col1 <- pal_futurama()(5)[1]
col2 <- pal_futurama()(5)[5]
col1 <- pal_aaas()(5)[4]
col2 <- pal_npg()(5)[5]
col1 <- exprBox_cond1Col
col2 <- exprBox_cond2Col

if(tad_to_plot == "chr17_TAD162") myWidthGG <- myWidthGG * 1.1


# cnv_dt <- get(load("../../v2_Yuanlong_Cancer_HiC_data_TAD_DA/CHECK_AKR1C_DEL/all_cnv_dt.Rdata"))
cnv_dt <- get(load("../../v2_Yuanlong_Cancer_HiC_data_TAD_DA/ENCSR489OCU_NCI-H460_40kb_TCGAluad_mutKRAS_mutEGFR_all_cnv_dt.Rdata"))
colnames(cnv_dt)[colnames(cnv_dt) == "value"] <- "cnv_status"

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

outFolder <- file.path("LOOK_TAD_EXPRESSION_WITHRANK_WITHCNVSTATUS_V3")
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

#count_file <- file.path(pipFolder, hicds, exprds, "1_runGeneDE", "DE_rnaseqDT.Rdata")
#count_file <- file.path(pipFolder, hicds, exprds, "0_prepGeneData", "rna_fpkmDT.Rdata")  # corrected here 15.03.2020
#stopifnot(file.exists(count_file))
#fpkm_dt <- get(load(count_file))
#stopifnot(names(geneList) %in% rownames(fpkm_dt))

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

save(withRank_toplot_dt2, file =file.path(outFolder, "withRank_toplot_dt2.Rdata"))

subTit <- paste0(tad_to_plot, " (rank: ", tad_plot_rank, ")")

# mutcols <- as.character(paletteer::paletteer_d("ggthemes::excel_Atlas"))[2:3]
mutcols <- setNames(c("#FC7715FF", "#AFBF41FF"), c("noMut", "withMut"))

withRank_toplot_dt2 <- withRank_toplot_dt2[order(withRank_toplot_dt2$cond_sh),]

withRank_toplot_dt2$cond <- as.character(withRank_toplot_dt2$cond)
stopifnot(withRank_toplot_dt2$cond %in% cnv_dt$cond)
withCNVwithRank_toplot_dt2 <- merge(withRank_toplot_dt2, cnv_dt, all.x=TRUE, all.y=FALSE, by=c("hicds","exprds","variable","cond","symbol"))

save(withCNVwithRank_toplot_dt2, file =file.path(outFolder, "withCNVwithRank_toplot_dt2.Rdata"))

withRank_toplot_dt2_s <- withRank_toplot_dt2

withRank_toplot_dt2 <- withCNVwithRank_toplot_dt2

withRank_toplot_dt2$cnv_status <- as.character(withRank_toplot_dt2$cnv_status)
withRank_toplot_dt2$cnv_status <- factor(withRank_toplot_dt2$cnv_status, levels=as.character(sort(as.numeric(as.character(unique(withRank_toplot_dt2$cnv_status))))))

withRank_toplot_dt2$cond <- factor(as.character(withRank_toplot_dt2$cond), levels=c(cond1,cond2))
stopifnot(!is.na(withRank_toplot_dt2$cond))

withRank_toplot_dt2$all_conds <- paste0(as.character(withRank_toplot_dt2$cond), "_", as.character(withRank_toplot_dt2$cond_sh))

all_genes <- unique(as.character(withRank_toplot_dt2$symbol))
mygene=all_genes[1]
# all_conds <- unique(as.character(withRank_toplot_dt2$cond))
# mycond=all_conds[1]
all_conds <- unique(as.character(withRank_toplot_dt2$all_conds))
mycond=all_conds[1]

cat("by cond\n")
table(withRank_toplot_dt2$symbol, withRank_toplot_dt2$cnv_status, withRank_toplot_dt2$cond)

cat("by all_conds\n")
table(withRank_toplot_dt2$symbol, withRank_toplot_dt2$cnv_status, withRank_toplot_dt2$all_conds)

for(i in unique(withRank_toplot_dt2$all_conds)) {
  cat(paste0("> ", i,"\n"))
  print(round(prop.table(table(withRank_toplot_dt2[withRank_toplot_dt2$all_conds==i,]$symbol, withRank_toplot_dt2[withRank_toplot_dt2$all_conds == i,]$cnv_status),1),4))
  
}
stop("-ok\n")
my_cols <- c("-2" = "darkgreen",
             "-1" = "orange",
             "0" = "grey",
             "1"= "darkorchid3",
             "2" ="dodgerblue")

for(mygene in all_genes) {
  
  sub_dt <- withRank_toplot_dt2[as.character(withRank_toplot_dt2$symbol) == mygene,]
  
  # sub_dt <- sub_dt[as.character(sub_dt$cnv_status) == "2" | as.character(sub_dt$cnv_status) == "-2" ,]
  
  for(mycond in all_conds){
    
    mysub <- paste0(mygene, " - ", mycond)
  
    p3 <- ggdensity(
          sub_dt[sub_dt$all_conds == mycond,],
          # sub_dt[sub_dt$cond == mycond,],
                    x = "value_log10",
                    y = "..density..",
                    # combine = TRUE,                  # Combine the 3 plots
                    xlab = "gene expr.",
                    # add = "median",                  # Add median line.
                    rug = FALSE,                      # Add marginal rug
                    color = "cnv_status",
                    fill = "cnv_status",
                    palette = "jco"
    ) +
      ggtitle(paste0(hicds_lab, " - ", exprds_lab), subtitle = paste0(mysub))+
      scale_color_manual(values=my_cols)+
      scale_fill_manual(values=my_cols)  +
      # labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
      guides(color=FALSE)+
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
      scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      theme(
        # text = element_text(family=fontFamily),
        panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
        panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major.x =  element_blank(),
        panel.grid.minor.x =  element_blank(),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
        axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
        axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
        plot.title = element_text(hjust=0.5, size = 16, face="bold"),
        plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
        legend.title = element_text(face="bold")
      ) 
    
    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", tad_to_plot, "_", mygene, "_", mycond, "_density.", plotType))
    ggsave(p3, file=outFile, height=5.5, width=7)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    p3 <- ggboxplot(
      sub_dt[sub_dt$all_conds == mycond,],
      # sub_dt[sub_dt$cond == mycond,],
      y = "value_log10",
      x = "cnv_status",
      # combine = TRUE,                  # Combine the 3 plots
      ylab = "gene expr.",
      # add = "median",                  # Add median line.
      xlab="CNV status",
      rug = FALSE,                      # Add marginal rug
      color = "cnv_status",
      # fill = "cnv_status",
      add = "jitter",
      palette = "jco"
    ) +
      ggtitle(paste0(hicds_lab, " - ", exprds_lab), subtitle = paste0(mysub))+
      scale_color_manual(values=my_cols)+
      # scale_fill_manual(values=my_cols)  +
      # labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
      guides(color=FALSE, shape=FALSE, FILL=FALSE)+
      scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
      # scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
      theme(
        # text = element_text(family=fontFamily),
        panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
        panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major.x =  element_blank(),
        panel.grid.minor.x =  element_blank(),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
        axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
        axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
        plot.title = element_text(hjust=0.5, size = 16, face="bold"),
        plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
        legend.title = element_text(face="bold")
      )  + rremove("legend")
    
    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", tad_to_plot, "_", mygene, "_", mycond, "_boxplot.", plotType))
    ggsave(p3, file=outFile, height=5.5, width=7)
    cat(paste0("... written: ", outFile, "\n"))
    
    
    
    
  }
}




myG_theme <-   theme( 
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
    axis.title.y = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.text = element_text(size=12),
    legend.key = element_blank(),
    legend.title = element_text(face="bold", size=12)
  )


p_var_boxplot <- ggplot(withRank_toplot_dt2, aes(x = symbol_lab, y = value_log10, col = cond)) + 
  geom_boxplot(notch = TRUE, outlier.shape=NA)+
  # geom_jitter(aes(fill =cond_sh, shape=cond_sh, group=cond), 
  #             alpha=0.5,
  #             position=position_jitterdodge())+
  # 
  geom_jitter(aes(fill =cnv_status, shape=cnv_status, group=cond), 
              alpha=0.5,
              position=position_jitterdodge())+
  
  
  ggtitle(paste0(hicds_lab, " - ", exprds_lab), subtitle = paste0(subTit))+
  scale_x_discrete(name=my_xlab)+
  scale_y_continuous(name=paste0(my_ylab),
                     breaks = scales::pretty_breaks(n = 20))+
  
  
  # scale_shape_manual(
  #   values = c(21,24),
  #   breaks = c("noMut", "withMut"),
  #   labels = c("not mut.", "mut.")
  # )+
  # 
  # scale_fill_manual(values=mutcols, 
  #                   breaks=c("noMut", "withMut"),
  #                   labels=c("not mut.", "mut."))+
  
  
  scale_color_manual(values=c(col1, col2))+
  # labs(fill  = paste0("KEAP1|NFE2L2"), color=paste0("Cond."), shape=paste0("KEAP1|NFE2L2")) +
  labs(fill  = paste0("CNV status"), color=paste0("Cond."), shape=paste0("CNV status")) +
	myG_theme


outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", tad_to_plot, "_allSamples_exprValues_boxplot_vShape.", plotType))
ggsave(plot = p_var_boxplot, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))

boxSpacing <- 0.2
mutSpacing <- 0.075

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


p_var_boxplot_v2 <- ggplot(withRank_toplot_dt2b, aes(x = box_pos, y = value_log10, group=box_pos, color=cond)) +
  ggtitle(paste0(hicds_lab, " - ", exprds_lab), subtitle = paste0(subTit))+
  geom_boxplot(notch = TRUE, outlier.shape=NA) + 
  # geom_jitter(aes(x=dot_pos, y = value_log10, group=dot_pos, shape=cond_sh, fill=cond_sh),
  #             alpha=0.5,
  #             width=0.05)+
  # 
  geom_jitter(aes(x=dot_pos, y = value_log10, group=dot_pos, shape=cnv_status, fill=cnv_status),
              alpha=0.5,
              width=0.05)+
  
  
  scale_color_manual(values=c(col1, col2))+
  
  # scale_fill_manual(values=mutcols, 
  #                   breaks=c("noMut", "withMut"),
  #                   labels=c("not mut.", "mut."))+
  # scale_shape_manual(
  #   values = c(21,24),
  #   breaks = c("noMut", "withMut"),
  #   labels = c("not mut.", "mut.")
  # )+
  
  # labs(fill  = paste0("KEAP1|NFE2L2"), color=paste0("Cond."), shape=paste0("KEAP1|NFE2L2")) +
  labs(fill  = paste0("CNV status"), color=paste0("Cond."), shape=paste0("CNV status")) +
  
  scale_x_continuous(name=my_xlab, breaks=unique(withRank_toplot_dt2b$symbol_pos), 
                     labels = as.character(levels(withRank_toplot_dt2b$symbol_lab)))+
  scale_y_continuous(name=paste0(my_ylab),
                     breaks = scales::pretty_breaks(n = 20))+
	myG_theme

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", tad_to_plot, "_allSamples_exprValues_boxplot_vShape_v2.", plotType))
ggsave(plot = p_var_boxplot_v2, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))

p_var_boxplot_v2 <- ggplot(withRank_toplot_dt2b, aes(x = box_pos, y = value_log10, group=box_pos, color=cond)) +
  ggtitle(paste0(hicds_lab, " - ", exprds_lab), subtitle = paste0(subTit))+
  geom_boxplot(notch = TRUE, outlier.shape=NA) + 
  
  # geom_jitter(aes(x=dot_pos, y = value_log10, group=dot_pos, shape=cond_sh, fill=cond), 
  #             alpha=0.5,
  #             width=0.05)+
  geom_jitter(aes(x=dot_pos, y = value_log10, group=dot_pos, shape=cnv_status, fill=cond), 
              alpha=0.5,
              width=0.05)+
  
  
  scale_color_manual(values=c(col1, col2))+
  scale_fill_manual(values=c(col1, col2))+
  # scale_fill_manual(values=mutcols, 
  #                   breaks=c("noMut", "withMut"),
  #                   labels=c("not mut.", "mut."))+
  # scale_shape_manual(
  #   values = c(21,24),
  #   breaks = c("noMut", "withMut"),
  #   labels = c("not mut.", "mut.")
  # )+
  
  # labs(fill  = paste0("KEAP1|NFE2L2"), color=paste0("Cond."), shape=paste0("KEAP1|NFE2L2")) +
  labs(fill  = paste0("CNV status"), color=paste0("Cond."), shape=paste0("CNV status")) +
  
  scale_x_continuous(name=my_xlab, breaks=unique(withRank_toplot_dt2b$symbol_pos), 
                     labels = as.character(levels(withRank_toplot_dt2b$symbol_lab)))+
  scale_y_continuous(name=paste0(my_ylab),
                     breaks = scales::pretty_breaks(n = 20))+
myG_theme

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", tad_to_plot, "_allSamples_exprValues_boxplot_vShape_v3.", plotType))
ggsave(plot = p_var_boxplot_v2, filename = outFile, height=myHeightGG, width = myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))


# paletteer::paletteer_d("ggthemes::Classic_Green_Orange_6") 1:2
# paletteer::paletteer_d("lisa::KarlZerbe") 4:5
# paletteer::paletteer_d("lisa::FrancescoXanto")4:5
# paletteer::paletteer_d("lisa::PabloPicasso_1")3:4
# paletteer::paletteer_d("lisa::ClaesOldenburg")3:4
# paletteer::paletteer_d("lisa::ClaudeMonet")2:3
# paletteer::paletteer_d("lisa::KazimirMalevich")3:4
# paletteer::paletteer_d("lisa::GustavKlimt")2:3
# paletteer::paletteer_d("lisa::HelenFrankenthaler")1:2
# paletteer::paletteer_d("lisa::RobertDelaunay")2:3
# paletteer::paletteer_d("lisa::JackBush")3:4






##############################
cat("***** DONE: ", script_name, "\n")

cat(paste0(startTime, "\n", Sys.time(), "\n"))

