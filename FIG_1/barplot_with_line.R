# Rscript barplot_with_line.R

require(ggplot2)
require(ggsci)
require(patchwork)
require(foreach)

plotType <- "svg"

source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")

outFolder <- "BARPLOT_WITH_LINE"
dir.create(outFolder, recursive = TRUE)


myHeight <- 5
myWidth <- 7

plotCex <- 1.4

interval_fcc <- c("<=0.5", "]0.5-1[", "1")
lineVar <- c( "]0.75, 1]")

dt2 <- get(load("BARPLOT_WITH_FCC_FRACT//all_dt.Rdata"))
dt2$dataset <- file.path(dt2$hicds, dt2$exprds)
# dt2$intervalFCC <- factor(dt2$intervalFCC, levels=interval_fcc)
stopifnot(!is.na(dt2$intervalFCC))
stopifnot(!duplicated(dt2))


dt3 <- get(load("BARPLOT_WITH_TOPFCC_FRACT///all_dt.Rdata"))
dt3$dataset <- file.path(dt3$hicds, dt3$exprds)

for(curr_ds in unique(dt3$dataset)) {
stopifnot(  sum(dt2$nFCC[dt2$dataset == curr_ds & dt2$intervalFCC %in% c("]0.5, 0.75]", "]0.75, 1]")]) == 
                  sum(dt3$countFCC[dt3$dataset == curr_ds]))
}


dt2$nFCC[dt2$dataset == "Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas" & dt2$intervalFCC %in% c("]0.5, 0.75]", "]0.75, 1]")] == 
dt3$countFCC[dt3$dataset == "Barutcu_MCF-10A_40kb/TCGAbrca_lum_bas"]


dt1 <- get(load("FCC_WAVE_PLOT_NOABS/all_fcc_dt.Rdata"))
dt1$dataset <- file.path(dt1$hicds, dt1$exprds)
dt1 <- dt1[order(dt1$fcc_auc, decreasing = TRUE),]
ds_levels <- as.character(dt1$dataset)
stopifnot(!duplicated(dt1))

dscols <- all_cols[all_cmps[basename(as.character(dt1$dataset))]]

dt2$dataset <- factor(dt2$dataset, levels=ds_levels)
stopifnot(!is.na(dt2$dataset))


plot_dt2 <- dt2[dt2$intervalFCC == lineVar,]
plot_dt2$dataset <- factor(plot_dt2$dataset, levels=ds_levels)
stopifnot(!is.na(plot_dt2$dataset))
plot_dt2 <- plot_dt2[order(plot_dt2$dataset),]
stopifnot(!is.na(plot_dt2))

my_main <- "FCC AUC ratio"
my_main2 <- paste0("Ratio TADs with FCC \u2208 ", lineVar)

linecol <- "brown"


############################## BAR COLS BY DATASET

outFile <- file.path(outFolder, paste0("fcc_barplot_coloredBars.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*1.2))

par(mar=par()$mar+c(2,0,0,2))

par(family=fontFamily)
barp <- barplot(dt1$fcc_auc-1,
                ylab="FCC AUC ratio", cex.lab=1.2, 
                main = my_main,
                # xlab="Datasets",
                cex.main = plotCex,
                # xlab=paste0("Datasets\n(n=", nrow(dt1), ")"),
                xlab="",
                col=dscols, axes=F)
axis(2, at = seq(0, 0.8, by=0.1), labels = seq(0, 0.8, by=0.1)+1)
mtext(1, text=paste0("Datasets\n(n=", nrow(dt1), ")"), line=2, cex=plotCex)

# add the line
par(new = T, family=fontFamily)

plot(x=barp,
     # ylim = c(0,1),
     xlab="", ylab="", lty=1,pch=16,lwd=2,
     y=plot_dt2$countFCC,  type="b", col = linecol,
     axes=FALSE)

axis(side=4, col = linecol, col.ticks = linecol, col.axis=linecol, at = seq(0, 1, by=0.05))
mtext(side = 4, line = 2, text=my_main2, col=linecol,  cex=plotCex)

# legend("bottom", pch=16, col=c(all_cols, linecol), 
#        legend=c(names(all_cols), paste0("Ratio\n", lineVar, " \nTADs")),
#         lty=c(rep(-1, length(all_cols), 1)),
#        cex = c(rep(plotCex, length(all_cols)), -1),
#        inset=c(0,-1), xpd=TRUE,
#        horiz = TRUE,
#        bty="n")
legend("bottom", pch=16, col=c(all_cols), 
       legend=c(names(all_cols)),
       lty=c(rep(-1, length(all_cols))),
       cex = c(rep(plotCex, length(all_cols))),
       inset=c(0,-0.5),
       xpd=TRUE,
       horiz = TRUE,
       bty="n")

foo <- dev.off()

cat(paste0("... written: ", outFile, "\n"))








