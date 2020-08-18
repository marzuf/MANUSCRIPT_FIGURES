
# Rscript example_fcc_dist_withPurity_final.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad
# Rscript example_fcc_dist.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich
# Rscript example_fcc_dist.R LG1_40kb TCGAluad_mutKRAS_mutEGFR # highest FCC AUC
# Rscript example_fcc_dist.R GSE105381_HepG2_40kb TCGAlihc_wt_mutCTNNB1
# Rscript example_fcc_dist.R ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf # lowest FCC AUC


require(doMC)
require(foreach)
require(ggplot2)
require(ggsci)
require(mapplots)
require(RColorBrewer)
require(ggpubr)

registerDoMC(40)

plotType <- "svg"

source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")


outFolder <- "EXAMPLE_FCC_DIST_WITHPURITY_FINAL"
dir.create(outFolder, recursive = TRUE)

buildData <- TRUE

hicds <- "ENCSR489OCU_NCI-H460_40kb"
exprds <- "TCGAluad_norm_luad"


args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
hicds <- args[1]
exprds <- args[2]
stopifnot(!is.na(hicds))
stopifnot(!is.na(exprds))

stopifnot(hicds %in% names(hicds_names))
stopifnot(exprds %in% names(exprds_names))

hicds_lab <- hicds_names[paste0(hicds)]
exprds_lab <- exprds_names[paste0(exprds)]


purity_dt <- get(load(file.path(runFolder,"ALLTADS_AND_PURITY_FINAL/aran/CPE/log10/all_ds_corrPurity_dt.Rdata" )))
purity_dt <- purity_dt[purity_dt$dataset == file.path(hicds, exprds),]
if(nrow(purity_dt) == 0){cat("unavailable purity data for this datastet\n");stop("-ok\n")}
agg_purity_dt <- aggregate(purityCorr~region+dataset, data=purity_dt, FUN=mean)
agg_purity_dt$hicds <- dirname(agg_purity_dt$dataset)
agg_purity_dt$exprds <- basename(agg_purity_dt$dataset)

fcc_file <- file.path(pipFolder, hicds, exprds, step8fcc_folder, "all_obs_prodSignedRatio.Rdata")
stopifnot(file.exists(fcc_file))
all_fcc <- get(load(fcc_file))

all_fcc <- sort(all_fcc, decreasing = TRUE)

fcc_breaks <- c(0.5, 1)

mycols <- setNames(pal_d3()(3), c("<=0.5", "]0.5, 1[", "1"))

mycols_id <- setNames(names(mycols), mycols)

fcc_cols <- ifelse(all_fcc <= fcc_breaks[1], mycols[1],
                   ifelse(all_fcc > fcc_breaks[1] & all_fcc < fcc_breaks[2], mycols[2],
                          ifelse(all_fcc == fcc_breaks[2], mycols[3], NA)))
stopifnot(!is.na(fcc_cols))

rel_rank <- c(1:length(all_fcc))/length(all_fcc)

myTit <- paste0(hicds_lab, " - ", exprds_lab)
my_xlab <- "Ranked TADs"
my_ylab <- "FCC"

all_FCC_dt <- data.frame(hicds=hicds, exprds=exprds, region=names(all_fcc), relative_rank=rel_rank, FCC = all_fcc, stringsAsFactors=FALSE)
saveFile <- file.path(outFolder, paste0("fig1G_", hicds, "_", exprds, "_all_FCC_dt.Rdata"))
save(all_FCC_dt, file=saveFile, version=2)
cat(paste0("... written:" , saveFile, "\n"))

all_FCC_dt <- merge(all_FCC_dt, agg_purity_dt, by=c("hicds", "exprds", "region"), all=FALSE)
stopifnot(nrow(all_FCC_dt) == length(all_fcc))  

rbPal <- colorRampPalette(c('red','blue'))
all_FCC_dt$purityCols <- rev(rbPal(10))[as.numeric(cut(all_FCC_dt$purityCorr,breaks = 10))]

all_FCC_dt <- all_FCC_dt[order(all_FCC_dt$FCC, decreasing = T),]

ybreaks <- seq(-1, 1, by=0.25)

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_fcc_ranked_with_pie.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
par(bty="L", family=fontFamily, mar=par()$mar+c(0,2,0,0))
plot(
  x=rel_rank, 
  y=all_FCC_dt$FCC, type="h",
  ylim = range(all_FCC_dt$FCC),
  # col = fcc_cols,
  col = all_FCC_dt$purityCols,
  main = myTit,
  xlab="",
  ylab="",
  # xlab=my_xlab,
  # ylab= my_ylab,
	cex.main=plotCex,
  cex.axis=plotCex,
  cex.lab=plotCex,
 axes=FALSE 
)
mtext(side=2, text = paste0(my_ylab), font = 2, cex = plotCex, line=4)
mtext(side=1, text = paste0(my_xlab), font=2, cex=plotCex, line=1)
axis(2, at = ybreaks, lwd=-1, lwd.ticks = 1, cex.axis=plotCex, las=2)
# abline(h = min(all_fcc))
# axis(1, tick = -1, labels = FALSE)
box(type="L")

fcc_dist <- setNames(as.numeric(table(fcc_cols)), names(table(fcc_cols)))
names(fcc_dist) <- mycols_id[names(fcc_dist)]

slices <- fcc_dist
lbls <- names(fcc_dist)
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
# pie(slices,labels = lbls, col=rainbow(length(lbls)),
#     main="Pie Chart of Countries") 

add.pie(z=pct, x=0.8,y=0.5, radius=0.3, col = mycols[names(pct)], labels = paste0(pct, "%"))

legend(
  # "top",
"bottomleft",
       bty="n",
       # inset = c(0,-0.05),
       seg.len=0.5,
       col = rev(mycols),
       lty=1,
       lwd=4,
       legend=rev(names(mycols)),
       # horiz=T
       x.intersp = 0.5,
y.intersp=0.8,
cex=plotCex
       )


foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))
