# Rscript aucFCC_nSignif.R 

require(doMC)
require(foreach)
require(ggplot2)
require(ggsci)

require(ggpubr)
require(colorRamps)

registerDoMC(40)

plotType <- "svg"

source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")


outFolder <- "AUCFCC_NSIGNIF"
dir.create(outFolder, recursive = TRUE)

buildData <- TRUE

tadSignifThresh <- 0.01
geneSignifThresh <- 0.01


ggsci_pal <- "lancet"
ggsci_subpal <- ""


myWidthGG <- 10
myHeightGG <- 6

all_cols <- sort(all_cols)

inDT <- get(load(file.path(runFolder, "GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata")))

geneSignif_dt <- inDT[,c("hicds", "exprds", "entrezID", "adj.P.Val")]
geneSignif_dt <- unique(geneSignif_dt)
agg_geneSignif_dt <- aggregate(adj.P.Val ~ hicds+exprds, FUN=function(x) sum(x <= geneSignifThresh), data=geneSignif_dt)
colnames(agg_geneSignif_dt)[3] <- "nSignifGenes"

tadSignif_dt <- inDT[,c("hicds", "exprds", "region", "tad_adjCombPval")]
tadSignif_dt <- unique(tadSignif_dt)
agg_tadSignif_dt <- aggregate(tad_adjCombPval ~ hicds+exprds, FUN=function(x) sum(x <= tadSignifThresh), data=tadSignif_dt)
colnames(agg_tadSignif_dt)[3] <- "nSignifTADs"

auc_ratio_file <- file.path("../FIG_1", "FCC_WAVE_PLOT_NOABS/all_fcc_dt.Rdata")
stopifnot(file.exists(auc_ratio_file))
x <- get(load(auc_ratio_file))
x$dataset <- file.path(x$hicds,  x$exprds)
x <- x[order(x$fcc_auc, decreasing=TRUE),]
fcc_ds_order <- x$dataset

dsCols <- all_cols[all_cmps[basename(fcc_ds_order)]]
stopifnot(!is.na(dsCols))

agg_tadSignif_dt$cmpType <- all_cmps[agg_tadSignif_dt$exprds]
agg_tadSignif_dt$dsCols <- all_cols[all_cmps[basename(agg_tadSignif_dt$exprds)]]

agg_geneSignif_dt$cmpType <- all_cmps[agg_geneSignif_dt$exprds]
agg_geneSignif_dt$dsCols <- all_cols[all_cmps[basename(agg_geneSignif_dt$exprds)]]


### FCC AUC ~ nSIGNIF TADs
m_tads_dt <- merge(agg_tadSignif_dt, x[,c("hicds", "exprds", "fcc_auc")], all.x=TRUE, all.y=TRUE)
stopifnot(!is.na(m_tads_dt))

my_x <- m_tads_dt$nSignifTADs
my_y <- m_tads_dt$fcc_auc
my_xlab <- "# signif. TADs"
my_ylab <- "FCC AUC ratio"

nDS <- length(unique(file.path(m_tads_dt$hicds, m_tads_dt$exprds)))

plotTit <- "# signif. TADs and FCC AUC ratio"
subTit <- paste0("adj. p-val <= ", tadSignifThresh, "; all datasets - n", nDS)

outFile <- file.path(outFolder, paste0("aucFCC_nSignifTADs.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="L")

plot(
  x= my_x,
  y = my_y,
  main = paste0(plotTit), 
  col = m_tads_dt$dsCols,
  xlab=my_xlab,
  ylab=my_ylab,
  pch=16,
  cex=0.7,
  cex.main=plotCex,
  cex.lab=plotCex,
  cex.axis=plotCex
)
mtext(side=3, text = subTit, font=3)
abline(lm(my_y~my_x), col="darkgrey", lty=2)
addCorr(x =my_x, y=my_y, legPos = "bottomright",bty="n")
legend("topleft",
       col=all_cols,
       pch=16,
       legend=cmp_names[names(all_cols)], bty="n")
       
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


### FCC AUC ~ nSIGNIF GENES
m_genes_dt <- merge(agg_geneSignif_dt, x[,c("hicds", "exprds", "fcc_auc")], all.x=TRUE, all.y=TRUE)
stopifnot(!is.na(m_genes_dt))


my_x <- m_genes_dt$nSignifGenes
my_y <- m_genes_dt$fcc_auc
my_xlab <- "# signif. genes"
my_ylab <- "FCC AUC ratio"

nDS <- length(unique(file.path(m_genes_dt$hicds, m_genes_dt$exprds)))

plotTit <- "# signif. genes and FCC AUC ratio"
subTit <- paste0("adj. p-val <= ", geneSignifThresh, "; all datasets - n", nDS)

outFile <- file.path(outFolder, paste0("aucFCC_nSignifGenes.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="L")

plot(
  x= my_x,
  y = my_y,
  main = paste0(plotTit), 
  col = m_genes_dt$dsCols,
  xlab=my_xlab,
  ylab=my_ylab,
  pch=16,
  cex=0.7,
  cex.main=plotCex,
  cex.lab=plotCex,
  cex.axis=plotCex
)
mtext(side=3, text = subTit, font=3)
abline(lm(my_y~my_x), col="darkgrey", lty=2)
addCorr(x =my_x, y=my_y, legPos = "bottomleft",bty="n")
legend("topright",
       col=all_cols,
       pch=16,
       legend=cmp_names[names(all_cols)], bty="n")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

m_genes_dt$nSignifGenes_log10 <- log10(m_genes_dt$nSignifGenes)
my_x <- m_genes_dt$nSignifGenes_log10
my_y <- m_genes_dt$fcc_auc
my_xlab <- "# signif. genes [log10]"
my_ylab <- "FCC AUC ratio"

nDS <- length(unique(file.path(m_genes_dt$hicds, m_genes_dt$exprds)))

plotTit <- "# signif. genes and FCC AUC ratio"
subTit <- paste0("adj. p-val <= ", geneSignifThresh, "; all datasets - n", nDS)

outFile <- file.path(outFolder, paste0("aucFCC_nSignifGenesLog10.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="L")

plot(
  x= my_x,
  y = my_y,
  main = paste0(plotTit), 
  col = m_genes_dt$dsCols,
  xlab=my_xlab,
  ylab=my_ylab,
  pch=16,
  cex=0.7,
  cex.main=plotCex,
  cex.lab=plotCex,
  cex.axis=plotCex
)
mtext(side=3, text = subTit, font=3)
abline(lm(my_y~my_x), col="darkgrey", lty=2)
addCorr(x =my_x, y=my_y, legPos = "bottomleft",bty="n")
legend("topright",
       col=all_cols,
       pch=16,
       legend=cmp_names[names(all_cols)], bty="n")

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




