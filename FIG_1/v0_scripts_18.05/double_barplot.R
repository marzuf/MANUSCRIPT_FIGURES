# Rscript double_barplot.R
plotType <- "svg"

source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")

outFolder <- "DOUBLE_BARPLOT"
dir.create(outFolder, recursive = TRUE)

myHeight <- 5
myWidth <- 7

plotCex <- 1.4

plotRange <-  "]0.75, 1]"

dt1_init <- get(load("BARPLOT_WITH_FCC_FRACT/all_dt.Rdata"))
dt2 <- get(load("FCC_WAVE_PLOT_NOABS/all_fcc_dt.Rdata"))

dt1 <- dt1_init[dt1_init$intervalFCC == plotRange,]
stopifnot(nrow(dt1) > 0)

mdt1 <-  reshape(dt1, idvar=c("hicds", "exprds"),direction = "wide", timevar="intervalFCC")
colnames(mdt1)[3] <- "ratioTopFCC" 

plot_dt <- merge(mdt1, dt2, by=c("hicds", "exprds"), all=T)
stopifnot(!is.na(plot_dt))

plot_dt <- plot_dt[order(plot_dt$fcc_auc, decreasing = TRUE),]

plot_dt$barcols <- all_cols[all_cmps[plot_dt$exprds]]

plot_dt$fcc_auc_resc <- plot_dt$fcc_auc - 1

# > range(plot_dt$ratioTopFCC)
# [1] 0.1950948 0.3180987
# > range(plot_dt$fcc_auc_resc)
# [1] 0.1605670 0.5318338

toplabs <- seq(from=0.1, to =0.7, by =0.1)

botlabs <- seq(from=0, to =0.4, by =0.05)

outFile <- file.path(outFolder, paste0("double_barplot_fccAUC_fccRange.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))

tit1 <- "FCC AUC ratio (PERMG2T)"
tit1 <- "FCC AUC ratio"
tit1 <- "FCC cumsum curve AUC ratio"

par(mfrow=c(2,1))

# par(cex=0.7, mai=c(0.1,0.1,0.2,0.1), mar=c(0.1, 0.1, 0.1, 0.1))
par(cex=0.7, mai=c(1,4,2,1), mar=c(1,5,2,1)) # c(bottom, left, top, right) w
# par(fig=c(0.05, 0.9, 0.6, 0.9 )) # c(x1, x2, y1, y2)
par(bty="L")
par(family=fontFamily)
barplot(plot_dt$fcc_auc_resc, col = plot_dt$barcols, 
        xlab = "",
        # ylab="",
        ylab="FCC AUC ratio",
        cex.main=plotCex, cex.axis=plotCex, cex.lab=plotCex,
        axes=FALSE)
axis(2, at=c(0, toplabs), labels = c("", format(toplabs+1, digits=3)), las=2, cex=plotCex)
# mtext(side=2,"FCC AUC ratio (PERMG2T)", line=2.75, cex=0.8)
# mtext(side=2,"FCC AUC ratio", line=2.75, cex=0.8)

legend("topright", legend=names(all_cols), text.col = all_cols, bty="n")

mtext(side=3, paste0(tit1))

tit2 <- paste0("Ratio TADs with FCC \u2208 ", plotRange)
      
# par(fig=c(0.05, 0.9, 0.05, 0.4 ), new=T)
barplot(plot_dt$ratioTopFCC, col = plot_dt$barcols,
        xlab = "",
        # ylab = paste0("Fract. of TADs in ", plotRange),
        # ylab = "",
        # ylab=paste0("Ratio TADs with FCC \u2208 ", plotRange),
        ylab=paste0("Ratio of TADs"),
        cex.main=plotCex, cex.axis=plotCex, cex.lab=plotCex,
        axes=FALSE)
axis(2, at=c(0, botlabs), labels=c("", format(botlabs, digits = 2)), las=2, cex=plotCex)
# mtext(side=2,paste0("Ratio TADs with FCC \u2208 ", plotRange), line=2.75, cex=0.8)

# mtext(side=3, paste0("Ratio TADs with high FCC scores"))
mtext(side=3, paste0(tit2))

addCorr(legPos = "topright", x=plot_dt$ratioTopFCC, y= plot_dt$fcc_auc, bty="n")

# mtext(side=1, paste0("Datasets ramlkedall datasets (n=", nrow(plot_dt), ")"))

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))
# fig
# A numerical vector of the form c(x1, x2, y1, y2) which gives the (NDC) coordinates of the figure region in the display region of the device. 
# If you set this, unlike S, you start a new plot, so to add to an existing plot use new = TRUE as well.
# mai
# A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.


# make labels and margins smaller
# par(cex=0.7, mai=c(0.1,0.1,0.2,0.1))
# Temperature <- airquality$Temp
# # define area for the histogram
# par(fig=c(0.1,0.7,0.3,0.9))
# hist(Temperature)
# # define area for the boxplot
# par(fig=c(0.8,1,0,1), new=TRUE)
# boxplot(Temperature)
# # define area for the stripchart
# par(fig=c(0.1,0.67,0.1,0.25), new=TRUE)
# stripchart(Temperature, method="jitter")


