options(scipen=100)

# v2 version
# -> one single script for all and for cmpType
# -> if share more than 80% of the genes -> merge conserved regions
# -> count as conserved in one dataset at the exprds level
# -> min conserved region

SSHFS=F

# Rscript purityFilter_density_sameSize.R



script_name <- "purityFilter_density_sameSize.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(GenomicRanges)
require(ggplot2)
require(ggpubr)
require(RColorBrewer)
require(reshape2)
# require(gplots)
registerDoMC(4)

source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")


buildTable <- TRUE

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 5)
myWidth <- 12
axisCex <- 1.4

### HARD-CODED - MAIN SETTINGS
# _final:  discussion 04.08.2020 Giovanni - keep aran CPE data, first vial only

purity_ds <- "aran"
pm <- "CPE"
purity_plot_name <- paste0("Aran - ", pm)

corMet <- "pearson"
transfExpr <- "log10"
corrPurityQtThresh <- 0.05
signif_column <- "adjPvalComb"
signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", signifThresh)


source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")

plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 5)
myWidth <- 10

buildSamp <- FALSE

script0_name <- "0_prepGeneData"

outFolder <- file.path("PURITYFILTER_DENSITY_SAMESIZE",  purity_ds, pm, transfExpr)
dir.create(outFolder, recursive = TRUE)


setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)


purity_file <- file.path(runFolder,"ALLTADS_AND_PURITY_FINAL", purity_ds, pm, transfExpr, "all_ds_corrPurity_dt.Rdata") # here _final INPUT
purityData <- get(load(purity_file))
agg_purity <- aggregate(purityCorr~dataset+region, FUN=mean, data=purityData)

result_file <- file.path(runFolder,"CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)

merge_dt <- merge(agg_purity, resultData, by=c("dataset", "region"))
merge_dt$signif <- merge_dt$adjPvalComb <= signifThresh
merge_dt$signif_lab <- ifelse(merge_dt$signif, paste0("adj. p-val <=", signifThresh), paste0("adj. p-val >", signifThresh) )
purityCorrThresh <- as.numeric(quantile(merge_dt$purityCorr[!merge_dt$signif], probs = corrPurityQtThresh ))

subTit <- paste0(corMet, "'s corr.", " - ", purity_plot_name, " data")
plotTit <- paste0("Purity corr. distribution")
myx_lab <- paste0(transfExpr, " expr. and purity correlation (meanTAD)")
outFile <- file.path(outFolder, paste0("exprPurityCorr_meanTAD_signif_notSignif_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(split(merge_dt$purityCorr, merge_dt$signif_lab),
               plotTit = plotTit, my_xlab = myx_lab)
lines(density(merge_dt$purityCorr), col="green")
abline(v=purityCorrThresh, col="blue")
legend("topleft", lty=1, lwd=2, col=c("green", "blue"), bty="n", legend=c("all", paste0(corrPurityQtThresh, "-qt non-signif. TADs\n(=", round(purityCorrThresh, 2), ")")))
mtext(side=3, text = subTit)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

stopifnot(!duplicated(file.path(merge_dt$hicds, merge_dt$exprds, merge_dt$region)))



# change here subsample by TAD size
require(stringr)

set.seed(20200908)


nSamp <- 1000

merge_dt$nGenes <- str_count(string=merge_dt$region_genes, pattern=",") + 1
stopifnot(merge_dt$nGenes >= 3)

signif_merge_dt <- merge_dt[merge_dt$signif,]
notSignif_merge_dt <- merge_dt[!merge_dt$signif,]

tad_size = unique(signif_merge_dt$nGenes)[1]

if(buildSamp) {
  allSamp_notSignif_merge_dt <- foreach(i = 1:nSamp, .combine='rbind') %dopar% {
    subSamp_notSignif_merge_dt <- foreach(tad_size = unique(signif_merge_dt$nGenes), .combine='rbind') %dopar% {
      tmp_signif_merge_dt <- signif_merge_dt[signif_merge_dt$nGenes == tad_size,]
      tmp_notSignif_merge_dt <- notSignif_merge_dt[notSignif_merge_dt$nGenes == tad_size,]
      
      nSignif <- nrow(tmp_signif_merge_dt)
      
      stopifnot(nrow(tmp_notSignif_merge_dt) >= nSignif)
      
      sample_rows <- sample(x=1:nrow(tmp_notSignif_merge_dt), size=nSignif)
      
      out_dt <- tmp_notSignif_merge_dt[sample_rows,]
      stopifnot(nrow(out_dt) == nrow(tmp_signif_merge_dt))
      out_dt
    }  
    stopifnot(sum(subSamp_notSignif_merge_dt$nGenes) == sum(signif_merge_dt$nGenes))
    subSamp_notSignif_merge_dt
  }
  stopifnot(nrow(allSamp_notSignif_merge_dt) == nSamp * nrow(signif_merge_dt))
  
  samp_merge_dt <- rbind(allSamp_notSignif_merge_dt, signif_merge_dt)
  
  save(samp_merge_dt, file=file.path(outFolder, "samp_merge_dt.Rdata"), version=2)
  
} else {
  samp_merge_dt <- get(load(file.path(outFolder, "samp_merge_dt.Rdata")))
}

purityCorrThresh <- as.numeric(quantile(samp_merge_dt$purityCorr[!samp_merge_dt$signif], probs = corrPurityQtThresh ))

subTit <- paste0(corMet, "'s corr.", " - ", purity_plot_name, " data (nSamp=", nSamp, ")")
plotTit <- paste0("Purity corr. distribution")
myx_lab <- paste0(transfExpr, " expr. and purity correlation (meanTAD)")
outFile <- file.path(outFolder, paste0("exprPurityCorr_meanTAD_signif_notSignif_subSampData_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(split(samp_merge_dt$purityCorr, samp_merge_dt$signif_lab),
               plotTit = plotTit, my_xlab = myx_lab)
lines(density(merge_dt$purityCorr), col="green")
abline(v=purityCorrThresh, col="blue")
legend("topleft", lty=1, lwd=2, col=c("green", "blue"), bty="n", legend=c("all", paste0(corrPurityQtThresh, "-qt non-signif. TADs\n(=", round(purityCorrThresh, 2), ")")))
mtext(side=3, text = subTit)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


notSignif_merge_dt <- merge_dt[!merge_dt$signif,]
notSignif_merge_dt$signif_lab <- "OBS. adj. p-val >0.01"

samp_and_obs_merge_dt <- rbind(samp_merge_dt, notSignif_merge_dt)

purityCorrThreshObs <- as.numeric(quantile(samp_and_obs_merge_dt$purityCorr[samp_and_obs_merge_dt$signif_lab == "OBS. adj. p-val >0.01"], probs = corrPurityQtThresh ))
purityCorrThreshRandom <- as.numeric(quantile(samp_and_obs_merge_dt$purityCorr[samp_and_obs_merge_dt$signif_lab == "adj. p-val >0.01"], probs = corrPurityQtThresh ))

stopifnot(round(purityCorrThreshObs, 2) == round(purityCorrThreshRandom, 2))

subTit <- paste0(corMet, "'s corr.", " - ", purity_plot_name, " data (nSamp=", nSamp, ")")
plotTit <- paste0("Purity corr. distribution")
myx_lab <- paste0(transfExpr, " expr. and purity correlation (meanTAD)")
outFile <- file.path(outFolder, paste0("exprPurityCorr_meanTAD_signif_notSignif_obs_and_subSampData_density.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot_multiDens(split(samp_and_obs_merge_dt$purityCorr, samp_and_obs_merge_dt$signif_lab),
               plotTit = plotTit, my_xlab = myx_lab)
lines(density(merge_dt$purityCorr), col="green")
abline(v=c(purityCorrThreshObs,purityCorrThreshRandom) , col="blue")
legend("topleft", lty=1, lwd=2, col=c("green", "blue"), bty="n", 
       legend=c("all obs.", paste0(corrPurityQtThresh, "-qt non-signif. TADs\n(OBS=", round(purityCorrThreshObs, 2), "; RANDOM=", round(purityCorrThreshRandom, 2), ")")))
mtext(side=3, text = subTit)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))





stop("--ok\n")





merge_dt <- merge(agg_purity, resultData, by=c("dataset", "region"))
merge_dt$signif <- merge_dt$adjPvalComb <= signifThresh
merge_dt$signif_lab <- ifelse(merge_dt$signif, paste0("adj. p-val <=", signifThresh), paste0("adj. p-val >", signifThresh) )
merge_dt_all <- merge_dt
merge_dt_all$signif <- "all"
merge_dt$signif <- ifelse(merge_dt$signif, "signif.", "not signif.")
plot_dt <- rbind(merge_dt, merge_dt_all)
tad_labs <- c( "signif.", "not signif.", "all")
plot_dt$signif <- factor(plot_dt$signif, levels=tad_labs)
my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], tad_labs)

legTitle <- "TADs"

plotTit <- "Distribution of TAD mean corr. with purity"

mySub <- paste0(corMet, "'s corr.", " - ", purity_plot_name, " data - TAD signif. thresh = ", signifThresh)

txt_xpos <- purityCorrThresh 
txt_ypos <- max(unlist(by(plot_dt, plot_dt$signif, function(x) density(x$purityCorr)$y))) - 0.1
txt_lab <- paste0(corrPurityQtThresh, "-qt non-signif. TADs")


all_purity_values_dt <- plot_dt
save(all_purity_values_dt, file=file.path(outFolder, "all_purity_values_dt.Rdata"), version=2)


p3 <- ggdensity(plot_dt,
          x = "purityCorr",
          y = "..density..",
          # combine = TRUE,                  # Combine the 3 plots
          xlab = "TAD mean corr. with purity",
          # add = "median",                  # Add median line.
          rug = FALSE,                      # Add marginal rug
          color = "signif",
          fill = "signif",
          palette = "jco"
) +
  ggtitle(plotTit, subtitle = mySub)+
  scale_color_manual(values=my_cols)+
  scale_fill_manual(values=my_cols)  +
  labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
  guides(color=FALSE)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(
    text = element_text(family=fontFamily),
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
  ) +
  geom_vline(xintercept = purityCorrThresh, linetype=2, colour="darkgrey") +
  annotate("text", x=txt_xpos, y=txt_ypos, label=txt_lab, fontface="italic", hjust=1)

outFile <- file.path(outFolder, paste0("all_signif_notSignif_meanCorr_dist_density.", plotType))
ggsave(p3, file=outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))

p2 <- ggdensity(plot_dt[plot_dt$signif !="all",],
          x = "purityCorr",
          y = "..density..",
          # combine = TRUE,                  # Combine the 3 plots
          xlab = "TAD mean corr. with purity", 
          # add = "median",                  # Add median line. 
          rug = FALSE,                      # Add marginal rug
          color = "signif", 
          fill = "signif",
          palette = "jco"
) +
  scale_color_manual(values=my_cols)+
  scale_fill_manual(values=my_cols)  +
  ggtitle(plotTit, subtitle = mySub)+
  labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") + 
  guides(color=FALSE)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(
    text = element_text(family=fontFamily),
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
  ) +
  geom_vline(xintercept = purityCorrThresh, linetype=2, colour="darkgrey")+
  annotate("text", x=txt_xpos, y=txt_ypos, label=txt_lab, fontface="italic", hjust=1)

outFile <- file.path(outFolder, paste0("signif_notSignif_meanCorr_dist_density.", plotType))
ggsave(p2, file=outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))



p2b <- ggdensity(plot_dt[plot_dt$signif !="not signif.",],
                x = "purityCorr",
                y = "..density..",
                # combine = TRUE,                  # Combine the 3 plots
                xlab = "TAD mean corr. with purity", 
                # add = "median",                  # Add median line. 
                rug = FALSE,                      # Add marginal rug
                color = "signif", 
                fill = "signif",
                palette = "jco"
) +
  scale_color_manual(values=my_cols)+
  ggtitle(plotTit, subtitle = mySub)+
  scale_fill_manual(values=my_cols)  +
  labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") + 
  guides(color=FALSE)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(
    text = element_text(family=fontFamily),
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
  ) +
  geom_vline(xintercept = purityCorrThresh, linetype=2, colour="darkgrey")+
  annotate("text", x=txt_xpos, y=txt_ypos, label=txt_lab, fontface="italic", hjust=1)

outFile <- file.path(outFolder, paste0("signif_all_meanCorr_dist_density.", plotType))
ggsave(p2b, file=outFile, height=myHeight, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))



