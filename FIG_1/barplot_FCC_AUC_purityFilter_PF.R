all_dt <- get(load("FCC_WAVE_PLOT_NOABS/all_fcc_dt.Rdata"))
pF_dt <-  get(load("FCC_WAVE_PLOT_NOABS_PURITYFILTER_FINAL_PF//all_fcc_dt.Rdata"))

# Rscript barplot_FCC_AUC_purityFilter_PF.R

plotType <- "svg"

source("../settings.R")
source("../full_dataset_names.R")
source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")

plotCex <- 1.2
plotType <- "svg"
myHeight <- 5
myWidth <- 5

outFolder <- "BARPLOT_FCC_AUC_PURITYFILTER_PF"
dir.create(outFolder, recursive = TRUE)

plot_dt <- merge(pF_dt, all_dt, by=c("hicds", "exprds"), all.x=T, all.y=F, suffix=c("_pF", "_noF"))
stopifnot(!is.na(plot_dt))

plot_dt$dataset <- paste0(hicds_names[plot_dt$hicds], " - ", exprds_names[plot_dt$exprds])
plot_dt <- plot_dt[order(plot_dt$fcc_auc_noF,  decreasing=TRUE),]
plot_dt$dataset <- factor(plot_dt$dataset, levels = plot_dt$dataset)

plot_dt$cmp <- all_cmps[plot_dt$exprds]
stopifnot(!is.na(plot_dt$cmp))

plot_dt$fcc_auc_minus1 <- plot_dt$fcc_auc_noF - 1
plot_dt$fcc_auc_PF_minus1 <- plot_dt$fcc_auc_pF - 1

p1_tit <- "FCC cumsum curve AUC ratio"
p1_sub <- ""
p1_ylab <- "FCC AUC ratio"
p1_legTitle <- ""
plotMargin <- c(0.1,0.1,0.1,0.1)

y_range <- seq(from=0, to=0.55, by=0.05)
y_labs <- y_range+1
y_labs <- format(y_labs, digits=3)

cmp_levels <- sort(cmp_names)

p <- ggplot(plot_dt, aes(x=dataset, y = fcc_auc_minus1, color=cmp, fill=cmp))+ 
  ggtitle(p1_tit, subtitle = p1_sub)+
  scale_y_continuous(breaks = y_range, labels = y_labs, name=p1_ylab)+
  scale_x_discrete(name="")+
  scale_color_manual(values=all_cols, labels=cmp_levels)+
  scale_fill_manual(values=all_cols, labels=cmp_levels)+
  labs(fill=p1_legTitle) +
  my_box_theme+
  geom_point()+
  geom_bar(aes(y=fcc_auc_PF_minus1), stat="identity")+
  guides(color = FALSE)+
  coord_cartesian(clip = 'off', expand=F)  + 
  theme(
    text = element_text(family=fontFamily),
    legend.text = element_text(size=12),
    legend.title = element_text(size=14),
    axis.text.x = element_blank(),
    plot.margin = unit(plotMargin, "lines"))

outFile <- file.path(outFolder, paste0("fcc_auc_obs_vs_PF_barplot.", plotType))
ggsave(p, filename = outFile, height=myHeight, width=myWidth*2)
cat(paste0("... written: ", outFile, "\n"))


auc_obs_pf_cmp_dt <- plot_dt
save(auc_obs_pf_cmp_dt, file = file.path(outFolder, "auc_obs_pf_cmp_dt.Rdata"))

my_x <- plot_dt[,c("fcc_auc_noF")]
my_y <- plot_dt[,c("fcc_auc_pF")]

my_cols <- all_cols[plot_dt$cmp]

outFile <- file.path(outFolder, paste0("fcc_auc_obs_vs_PF_scatter.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth))
plot(
  x=my_x,
  y =my_y,
  xlab="FCC AUC ratio",
  ylab="FCC AUC ratio - with purity filter",
  main = "FCC AUC ratio",
  pch=16, 
  col= my_cols,
  cex.main=plotCex,
  cex.lab=plotCex,
  cex.axis=plotCex,
  cex=1.2
)
curve(1*x, lty=2, col="darkgrey", add=TRUE)
addCorr(x=my_x, y=my_y, bty="n", legPos = "topleft")
mtext(side=3, text=paste0("all TADs vs. purity-filtered TADs"))
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))




