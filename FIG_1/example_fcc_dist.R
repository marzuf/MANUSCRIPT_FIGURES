
# Rscript example_fcc_dist.R ENCSR489OCU_NCI-H460_40kb TCGAluad_norm_luad

require(doMC)
require(foreach)
require(ggplot2)
require(ggsci)
require(mapplots)

require(ggpubr)

registerDoMC(40)

plotType <- "svg"

source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")


outFolder <- "EXAMPLE_FCC_DIST"
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

ybreaks <- seq(-1, 1, by=0.25)

outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_fcc_ranked_with_pie.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))
par(bty="L", family=fontFamily, mar=par()$mar+c(0,2,0,0))
plot(
  x=rel_rank, 
  y=all_fcc, type="h",
  ylim = range(all_fcc),
  col = fcc_cols,
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
