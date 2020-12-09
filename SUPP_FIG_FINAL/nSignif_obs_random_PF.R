# SuppFig2a: On these random chromatin domains, DADo invariably returned a lower number of significant hits (Fig. 2b and Supplementary Fig. 2a).

# Rscript nSignif_obs_random_PF.R

# NB: here I  plot the # signif once the immune domains discard
# to have the # signif after runing DADo on filtered data -> look at nbr_signif_v2_purit


require(reshape2)
require(ggpubr)
require(ggsci)

outFolder <- file.path("NSIGNIF_OBS_RANDOM_PF")
dir.create(outFolder)

runFolder <- "../../v2_Yuanlong_Cancer_HiC_data_TAD_DA/"

randomPattern <- "RANDOMMIDPOSDISC"

pthresh <- 0.01

rd_col <- paste0("nSignifNotFlagged_", randomPattern)
obs_col <-  paste0("nSignifNotFlagged_obs")

plotCex <- 1.2

random_dt <- get(load(file.path(runFolder, "SIGNIF_PURITYFLAGGED_FINAL_RANDOMMIDPOS/aran/CPE/log10/all_dt.Rdata")))
random_dt$hicds <- dirname(random_dt$dataset)
random_dt <- random_dt[grepl(randomPattern, random_dt$hicds),]
random_dt$hicds <- gsub(paste0("_", randomPattern), "", random_dt$hicds)
random_dt$exprds <- basename(random_dt$dataset)
random_dt$dataset <- file.path(random_dt$hicds, random_dt$exprds)
random_dt[,rd_col] <- random_dt$nSignif - random_dt$nSignifAndFlagged


obs_dt <- get(load(file.path(runFolder, "SIGNIF_PURITYFLAGGED_FINAL/aran/CPE/log10/all_dt.Rdata")))
obs_dt$hicds <- dirname(obs_dt$dataset)
obs_dt$exprds <- basename(obs_dt$dataset)
obs_dt$dataset <- file.path(obs_dt$hicds, obs_dt$exprds)
obs_dt[,obs_col] <- obs_dt$nSignif - obs_dt$nSignifAndFlagged

range(random_dt[,rd_col])
range(obs_dt[,obs_col] )

m_dt <- merge(random_dt[,c("hicds", "exprds","dataset", rd_col)], obs_dt[,c("hicds", "exprds","dataset", obs_col)], by=c("hicds", "exprds","dataset"), all=TRUE)
stopifnot(!is.na(m_dt))

save(m_dt, version=2, file=file.path(outFolder, "m_dt.Rdata"))

x_lab <- paste0("# of signif. domains - observed data")
y_lab <- paste0("# of signif. domains - random data")

outFile <- file.path(outFolder, "nSignif_obs_random_PF_scatterplot.svg")
svg(outFile, height=5, width=5)
plot(
  x = m_dt[,obs_col],
  y = m_dt[,rd_col],
  xlab=x_lab,
  ylab=y_lab,
  pch=16,
  cex=0.9,
  cex.lab=plotCex,
  cex.axis=plotCex,
  cex.main=plotCex
)
curve(1*x, lty=2, col="red", add=T)
legend("topleft", legend="y=x", bty="n", text.col="red")
foo <- dev.off()



m_dt$dataset <- file.path(m_dt$hicds, m_dt$exprds)
m_dt <- m_dt[order(m_dt[,obs_col], m_dt[,rd_col], decreasing = T),]


obsCol <- pal_nejm()(2)[1]
rdCol <- pal_nejm()(2)[2]
lineCol <- "grey"

all_x <- 1:nrow(m_dt)

outFile <- file.path(outFolder, "nSignif_obs_random_connected_PF_dotplot.svg")
svg(outFile, height=6, width=8)
plot(NULL,
     xlim=range(all_x),
     ylim=range(c(m_dt[,rd_col], m_dt[,obs_col])),
     cex.lab=plotCex,
     cex.axis=plotCex,
     cex.main=plotCex,
     ylab = "# signif. differentially active domains",
     xlab=paste0(length(unique(m_dt$dataset)), " comparisons"),
     axes=FALSE
)
legend("topright", col =c(obsCol, rdCol), pch=16, bty="n", legend=c("observed", "random"))
axis(2)
axis(1, lwd.ticks = 0, labels = FALSE)
segments(x0=all_x, y0=m_dt[,obs_col],
         x1=all_x, y1=m_dt[,rd_col],
         col = lineCol)
points(x = all_x, y = m_dt[,obs_col],col=obsCol, pch=16)
points(x = all_x , y = m_dt[,rd_col],col=rdCol, pch=16)


foo <- dev.off()





plot_dt <- melt(m_dt, id=c("hicds", "exprds", "dataset"))
plot_dt$variable <- as.character(plot_dt$variable)
plot_dt$variable[plot_dt$variable ==rd_col ] <- "random"
plot_dt$variable[plot_dt$variable ==obs_col ] <- "observed"
plot_dt$dataset <- factor(plot_dt$dataset, levels = m_dt$dataset)

stopifnot(!is.na(plot_dt))

plotTit <- ""
subTit <- paste0("")

p <- ggbarplot(plot_dt, x="dataset", y="value", fill="variable", 
               xlab = "", ylab = "# signif. domains", 
               position=position_dodge(0.9))+
  scale_fill_nejm() + 
  # geom_bar(stat='identity', position='dodge')
  labs(fill="Data: ") + 
  ggtitle(plotTit, subtitle=subTit)+
  scale_y_continuous(expand = c(0, 0), limits = c(0,65)) + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    plot.title = element_text(size=16, face = "bold", hjust=0.5),
    plot.subtitle = element_text(size=14, face = "italic", hjust=0.5)
  )

outFile <- file.path(outFolder, paste0("nSignif_obs_random_PF_barplot.svg"))
ggsave(p, filename=outFile, height=5, width=8)


