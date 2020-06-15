# NEW VERSION: RANK BY MINIMAL # OF SAMPLES OR BY SUM # OF SAMPLES


plotType <- "svg"

require(ggsci)

source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")
require(ggsci)

# Rscript nSignifGenes_nSignifTADs_minSamp.R

cmp_levels <- sort(cmp_names)

outFolder <- file.path("NSIGNIFGENES_NSIGNIFTADS_MINSAMP")
dir.create(outFolder, recursive = TRUE)

dotpch <- 19
# segcol <-  "#BEBEBE19"
segcol <- "grey"
dotCex <- 1.1


setDir <- ""
mainFolder <- runFolder
settingFolder <- file.path(mainFolder, "PIPELINE", "INPUT_FILES")


inDT <- get(load(file.path(runFolder, "GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata")))
geneDT <- inDT[,c("hicds", "exprds", "entrezID", "adj.P.Val")]
geneDT <- unique(geneDT)
nSignifGenes_dt <- aggregate(adj.P.Val~hicds + exprds, data = geneDT, FUN=function(x) sum(x<=geneSignifThresh))
colnames(nSignifGenes_dt)[colnames(nSignifGenes_dt) == "adj.P.Val"] <- "nSignifGenes"

ratioSignifGenes_dt <- aggregate(adj.P.Val~hicds + exprds, data = geneDT, FUN=function(x) mean(x<=geneSignifThresh))
colnames(ratioSignifGenes_dt)[colnames(ratioSignifGenes_dt) == "adj.P.Val"] <- "ratioSignifGenes"

tadDT <- inDT[,c("hicds", "exprds", "region", "tad_adjCombPval")]
tadDT <- unique(tadDT)
nSignifTADs_dt <- aggregate(tad_adjCombPval~hicds + exprds, data = tadDT, FUN=function(x) sum(x<=tadSignifThresh))
colnames(nSignifTADs_dt)[colnames(nSignifTADs_dt) == "tad_adjCombPval"] <- "nSignifTADs"

nSignif_dt <- merge(nSignifGenes_dt, nSignifTADs_dt, by=c("hicds", "exprds"), all=TRUE)
stopifnot(!is.na(nSignif_dt))

nSignif_dt <- nSignif_dt[order(nSignif_dt$nSignifTADs, decreasing = TRUE),]


ratioSignifTADs_dt <- aggregate(tad_adjCombPval~hicds + exprds, data = tadDT, FUN=function(x) mean(x<=tadSignifThresh))
colnames(ratioSignifTADs_dt)[colnames(ratioSignifTADs_dt) == "tad_adjCombPval"] <- "ratioSignifTADs"

ratioSignif_dt <- merge(ratioSignifGenes_dt, ratioSignifTADs_dt, by=c("hicds", "exprds"), all=TRUE)
stopifnot(!is.na(ratioSignif_dt))

nSignif_dt <- merge(nSignifGenes_dt, nSignifTADs_dt, by=c("hicds", "exprds"), all=TRUE)
stopifnot(!is.na(nSignif_dt))

nSignif_dt <- nSignif_dt[order(nSignif_dt$nSignifTADs, decreasing = TRUE),]

# retrieve the number of samples
ratioSignif_dt$dataset <- file.path(ratioSignif_dt$hicds, ratioSignif_dt$exprds)
ratioSignif_dt$minNbrSample <- sapply(ratioSignif_dt$dataset, function(x) {
  settingFile <- file.path(settingFolder, dirname(x), paste0("run_settings_", basename(x), ".R"))
  stopifnot(file.exists(settingFile))
  source(settingFile)
  samp1 <- get(load(file.path(setDir, sample1_file)))
  samp2 <- get(load(file.path(setDir, sample2_file)))
  min(c(length(samp1), length(samp2)))
})
ratioSignif_dt$dataset <- paste0(ratioSignif_dt$hicds, "\n", ratioSignif_dt$exprds)
ratioSignif_dt <- ratioSignif_dt[order(ratioSignif_dt$minNbrSample),]



# retrieve the number of samples

nSignif_dt$dataset <- file.path(nSignif_dt$hicds, nSignif_dt$exprds)

nSignif_dt$minNbrSample <- sapply(nSignif_dt$dataset, function(x) {
  settingFile <- file.path(settingFolder, dirname(x), paste0("run_settings_", basename(x), ".R"))
  stopifnot(file.exists(settingFile))
  source(settingFile)
  samp1 <- get(load(file.path(setDir, sample1_file)))
  samp2 <- get(load(file.path(setDir, sample2_file)))
  min(c(length(samp1), length(samp2)))
})

nSignif_dt$dataset <- paste0(nSignif_dt$hicds, "\n", nSignif_dt$exprds)


nSignif_dt <- nSignif_dt[order(nSignif_dt$minNbrSample),]

outFile <- file.path(outFolder, "nSignif_dt_minSample.Rdata")
save(nSignif_dt, file = outFile, version=2)


labcols <- all_cols[all_cmps[nSignif_dt$exprds]]

maxTADs <- max(ceiling(nSignif_dt$nSignifTADs/10)*10)
maxGenes <- max(ceiling(nSignif_dt$nSignifGenes/1000)*1000)

nSignif_dt$nSignifTADs_rescaled <- nSignif_dt$nSignifTADs/maxTADs * maxGenes
nSignif_dt$nSignifGenes_rescaled <- nSignif_dt$nSignifGenes/maxGenes * maxTADs

outFile <- file.path(outFolder, paste0("nSignifGenes_nSignifTADs_all_ds_withSymb_geneSignif",geneSignifThresh, "_tadSignif", tadSignifThresh, ".", plotType))
do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth*2))
dev.control(displaylist="enable")
par(bty="U", family=fontFamily)
par(mar = c(5,5,2,5))
plot(
  NULL,
  # x = 1:nrow(nSignif_dt),
  # y = nSignif_dt$nSignifTADs, 
  # col = tad_signif_col,
  xlim = c(1, nrow(nSignif_dt)),
  ylim=c(0, maxTADs),
  ylab = "# signif. TADs",
  xlab = "",
  pch=dotpch,
  cex = dotCex,
  main = paste0("# signif. features"),
  cex.axis=plotCex,
  cex.main = plotCex,
  cex.lab = plotCex,
  col.lab = tad_signif_col,
  axes = FALSE
)


segments(
  x0= 1:nrow(nSignif_dt),
  y0 = nSignif_dt$nSignifGenes_rescaled,
  x1= 1:nrow(nSignif_dt),
  y1=nSignif_dt$nSignifTADs,
  col = segcol
)

points(
  x = 1:nrow(nSignif_dt),
  y = nSignif_dt$nSignifTADs, 
  col = tad_signif_col,
  cex = dotCex,
  pch=dotpch
)

abline(lm(nSignif_dt$nSignifTADs~c(1:nrow(nSignif_dt))), col = tad_signif_col, lty=2)

legend("topleft",
       legend=c(paste0("gene signif.: adj. p-val <= ", geneSignifThresh),
                paste0("TAD signif.: adj. p-val <= ", tadSignifThresh)),
       bty="n")



mtext(side=3, line=-1, text=paste0("all datasets - n= ", length(unique(nSignif_dt$dataset))))
# mtext(side=1, col = labcols, text = nSignif_dt$dataset, at= 1:nrow(nSignif_dt), las=2, cex =0.5)
axis(2, col = tad_signif_col, col.ticks = tad_signif_col, col.axis=tad_signif_col, at=seq(from=0, to = maxTADs, by=10))
axis(1, labels=F, lwd.ticks = -1)
# axis(1, labels=F, at=1:nrow(nSignif_dt))
par(new = T, family=fontFamily)
plot(
  x = 1:nrow(nSignif_dt),
  y = nSignif_dt$nSignifGenes,
  col = gene_signif_col,
  ylab = NA,
  ylim=c(0, maxGenes),
  xlab = NA,
  pch=dotpch,
  cex.axis=plotCex,
  cex.main = plotCex,
  cex.lab = plotCex,
  axes = FALSE
)
axis(side=4, col = gene_signif_col, col.ticks = gene_signif_col, col.axis=gene_signif_col, at = seq(from=0, to=maxGenes, by=1000))
mtext(side = 4, line = 3, '# signif. genes', col=gene_signif_col,  cex=plotCex)

abline(lm(nSignif_dt$nSignifGenes~c(1:nrow(nSignif_dt))), col = gene_signif_col, lty=2)


legend("bottom", 
       # legend = paste0(labsymbol, " ", names(all_cols)),
       legend = paste0(names(all_cols)),
       col=all_cols,
       pch=15,
       # lty=c(1,2),
       horiz=TRUE,
       inset=c(0,-0.12),
       cex = plotCex,
       xpd=TRUE, bty="n"
)
signifPlot <- recordPlot()



invisible(dev.off())
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("nSignifGenes_nSignifTADs_all_ds_withLeg_geneSignif",geneSignifThresh, "_tadSignif", tadSignifThresh, ".", plotType))
do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth*2.5))
par(mar = c(5,5,2,5))
replayPlot(signifPlot) 
mtext(side=1, col = labcols, text = nSignif_dt$dataset, at= 1:nrow(nSignif_dt), las=2, cex = 0.6)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




outFile <- file.path(outFolder, paste0("nSignifGenes_nSignifTADs_all_ds_withNbrSamp_geneSignif",geneSignifThresh, "_tadSignif", tadSignifThresh, ".", plotType))
do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth*2.5))
par(mar = c(5,5,2,5))
replayPlot(signifPlot) 
mtext(side=1, col = labcols, text = paste0(nSignif_dt$minNbrSample, " ", labsymbol), at= 1:nrow(nSignif_dt),  las=2,cex = 0.9, line = 0)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))



dotcols <- all_cols[all_cmps[nSignif_dt$exprds]]

outFile <- file.path(outFolder, paste0("nSignifGenes_vs_nSignifTADs_all_ds_",geneSignifThresh, "_tadSignif", tadSignifThresh, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="L")
plot(
  x = nSignif_dt$nSignifTADs,
  y = nSignif_dt$nSignifGenes,
  col = dotcols,
  ylab = "# signif. genes",
  xlab = "# signif. TADs",
  pch=16,
cex=0.7,
	main = paste0("# signif. features"),
  cex.axis=plotCex,
  cex.main = plotCex,
  cex.lab = plotCex
)
mtext(side=3, text = paste0("all DS - n=", nrow(nSignif_dt), "; gene adj. p-val <= ", geneSignifThresh, " - TAD adj. p-val <= ", tadSignifThresh))
addCorr(x=nSignif_dt$nSignifTADs,y=nSignif_dt$nSignifGenes,bty="n")
#legend("bottomleft", 
#       # legend = paste0(labsymbol, " ", names(all_cols)),
#       legend = paste0(names(all_cols)),
#       col=all_cols,
#       pch=16,
#       cex = plotCex,
# bty="n"
#)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("nSignifGenes_vs_minNbrSample_all_ds_",geneSignifThresh, "_tadSignif", tadSignifThresh, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="L")
plot(
  x = nSignif_dt$minNbrSample,
  y = nSignif_dt$nSignifGenes,
  col = dotcols,
  ylab = "# signif. genes",
  xlab = "Min # of samples",
  pch=16,
cex=0.7,
	main = paste0("# signif. genes vs. sample size"),
  cex.axis=plotCex,
  cex.main = plotCex,
  cex.lab = plotCex
)
mtext(side=3, text = paste0("all DS - n=", nrow(nSignif_dt), "; gene adj. p-val <= ", geneSignifThresh))
addCorr(x=nSignif_dt$minNbrSample,y=nSignif_dt$nSignifGenes,bty="n")
#legend("bottomleft", 
#       # legend = paste0(labsymbol, " ", names(all_cols)),
#       legend = paste0(names(all_cols)),
#       col=all_cols,
#       pch=16,
#       cex = plotCex,
# bty="n"
#)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


outFile <- file.path(outFolder, paste0("nSignifTADs_vs_minNbrSample_all_ds_",geneSignifThresh, "_tadSignif", tadSignifThresh, ".", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
par(bty="L")
plot(
  x = nSignif_dt$minNbrSample,
  y = nSignif_dt$nSignifTADs,
  col = dotcols,
  ylab = "# signif. TADs",
  xlab = "Min # of samples",
  pch=16,
cex=0.7,
	main = paste0("# signif. TADs vs. sample size"),
  cex.axis=plotCex,
  cex.main = plotCex,
  cex.lab = plotCex
)
mtext(side=3, text = paste0("all DS - n=", nrow(nSignif_dt), "; TAD adj. p-val <= ", tadSignifThresh))
addCorr(x=nSignif_dt$minNbrSample,y=nSignif_dt$nSignifGenes,bty="n")
#legend("bottomleft", 
#       # legend = paste0(labsymbol, " ", names(all_cols)),
#       legend = paste0(names(all_cols)),
#       col=all_cols,
#       pch=16,
#       cex = plotCex,
# bty="n"
#)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


# save(nSignif_dt, version=2, file="nSignif_dt.Rdata")
# save(ratioSignif_dt, version=2, file="ratioSignif_dt.Rdata")
# load("nSignif_dt.Rdata")
# load("ratioSignif_dt.Rdata")

nSignif_dt$nSignifGenes_log10 <- log10(nSignif_dt$nSignifGenes)

nSignif_dt$cmpType <- all_cmps[paste0(nSignif_dt$exprds)]
stopifnot(!is.na(nSignif_dt$cmpType))
ratioSignif_dt$cmpType <- all_cmps[paste0(ratioSignif_dt$exprds)]
stopifnot(!is.na(ratioSignif_dt$cmpType))

sampBreaks <- c(30, 80, 130, 180)

myHeightGG <- 5
myWidthGG <- 7

plotTit <- paste0("# signif. features - all datasets - n= ", length(unique(nSignif_dt$dataset)))
subTit <- paste0("TAD signif. adj. p-val <= ", tadSignifThresh, "; gene signif. adj. p-val <= ", geneSignifThresh)

# 
# ggplot(ratioSignif_dt, aes(y=ratioSignifTADs, x=ratioSignifGenes, size=sumSample))+
#   coord_polar(theta = "x")+
# # ggplot(ratioSignif_dt, aes(x=ratioSignifTADs, y=ratioSignifGenes, size=sumSample))+
# #   coord_polar(theta = "x") +
#   ggtitle(plotTit, subtitle = subTit)+
#   scale_y_continuous( breaks = scales::pretty_breaks(n = 10))+
#   scale_x_continuous( breaks = scales::pretty_breaks(n = 10))+
#   labs(size="# samples\n(samp1+samp2)", x = "Ratio signif. TADs", y="Ratio signif. genes")+
#   scale_size_continuous(breaks = sampBreaks)+
#   geom_point() + 
#   my_box_theme +
#   theme(
#     plot.subtitle = element_text(size=12, hjust=0.5, face="italic"),
#     axis.line=element_line(),
#     legend.text=element_text(size=12),
#     legend.title=element_text(size=14)
#   )
#   
  
p_ratio <-  ggplot(ratioSignif_dt, aes(x=ratioSignifTADs, y=ratioSignifGenes, size=minNbrSample, color=cmpType))+
  ggtitle(plotTit, subtitle = subTit)+
  scale_y_continuous( breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous( breaks = scales::pretty_breaks(n = 10))+
  scale_color_manual(values=all_cols, labels=cmp_levels)+
  labs(size="Min(samp1, samp2)\n# samples", x = "Ratio signif. TADs", y="Ratio signif. genes", color="")+
  scale_size_continuous(breaks = sampBreaks)+
  geom_point(alpha=0.6) + 
  my_box_theme +
  theme(
    plot.subtitle = element_text(size=12, hjust=0.5, face="italic"),
    axis.line=element_line(),
    legend.key = element_rect(fill = NA),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14)
    )

outFile <- file.path(outFolder, paste0("ratioGenesSignif_ratioTADsSignif_dotSize_plot.", plotType))
ggsave(p_ratio, filename = outFile, height=myHeightGG, width=myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))


p_nbr <- ggplot(nSignif_dt, aes(x=nSignifTADs, y=nSignifGenes_log10, size=minNbrSample, color=cmpType))+
  ggtitle(plotTit, subtitle = subTit)+
  scale_color_manual(values=all_cols, labels=cmp_levels)+
  geom_point(alpha=0.6) + 
  scale_y_continuous( breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous( breaks = scales::pretty_breaks(n = 10))+
  labs(size="Min(samp1,samp2)\n# samples", x = "# signif. TADs", y="# signif. genes [log10]", color="")+
  scale_size_continuous(breaks = sampBreaks)+
  my_box_theme+
  theme(
    plot.subtitle = element_text(size=12, hjust=0.5, face="italic"),
    legend.key = element_rect(fill = NA),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14),
    axis.line=element_line()
  )
outFile <- file.path(outFolder, paste0("nbrGenesSignif_nbrTADsSignif_dotSize_plot.", plotType))
ggsave(p_nbr, filename = outFile, height=myHeightGG, width=myWidthGG*1.2)
cat(paste0("... written: ", outFile, "\n"))





