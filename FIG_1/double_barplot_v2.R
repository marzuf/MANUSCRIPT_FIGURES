# Rscript double_barplot_v2.R

require(ggplot2)
require(ggsci)
require(patchwork)

plotType <- "svg"

source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")

outFolder <- "DOUBLE_BARPLOT_V2"
dir.create(outFolder, recursive = TRUE)

ggsci_pal <- "d3"
ggsci_subpal <- ""

mycols <- setNames(pal_d3()(3), c("<=0.5", "]0.5, 1[", "1"))
mycols <- mycols[2:3]

cmp_levels <- sort(cmp_names)

myHeight <- 5
myWidth <- 8

plotCex <- 1.4

interval_fcc <- c("]0.5, 1[", "1")

dt2 <- get(load("BARPLOT_WITH_TOPFCC_FRACT/all_dt.Rdata"))
dt2$dataset <- file.path(dt2$hicds, dt2$exprds)
dt2$intervalFCC <- factor(dt2$intervalFCC, levels=rev(interval_fcc))
stopifnot(!is.na(dt2$intervalFCC))


dt1 <- get(load("FCC_WAVE_PLOT_NOABS/all_fcc_dt.Rdata"))
dt1$dataset <- file.path(dt1$hicds, dt1$exprds)
dt1 <- dt1[order(dt1$fcc_auc, decreasing = TRUE),]
ds_levels <- as.character(dt1$dataset)

dt2$dataset <- factor(dt2$dataset, levels=ds_levels)
stopifnot(!is.na(dt2$dataset))

dt1$dataset <- factor(dt1$dataset, levels=ds_levels)
stopifnot(!is.na(dt1$dataset))

m_dt <- merge(dt1, dt2, by=c("dataset"))

dotcols <- all_cols[all_cmps[basename(as.character(ds_levels))]]
stopifnot(!is.na(dotcols))

dt1$fcc_auc_minus1 <- dt1$fcc_auc - 1

p1_tit <- "FCC AUC ratio"
p1_tit <- ""
p1_tit <- "FCC cumsum curve AUC ratio"
p1_sub <- ""
p1_xlab <- ""
p1_ylab <- "FCC AUC ratio"


p2_tit <- "Ratio top-FCC TADs"
p2_tit <- ""
p2_tit <- paste0("Ratio TADs in FCC range")
p2_sub <- ""
p2_xlab <- ""
p2_ylab <- "Ratio of TADs"

p1_legTitle <- ""
p2_legTitle <- "FCC range"

plotMargin <- c(0.1,0.1,0.1,0.1)

y_range <- seq(from=0, to=0.55, by=0.05)
y_labs <- y_range+1
y_labs <- format(y_labs, digits=3)

dt1$cmp <- all_cmps[basename(as.character(dt1$dataset))]
stopifnot(!is.na(dt1$cmp))


p1_aucRatio_plot <- ggplot(dt1, aes(x=dataset, y = fcc_auc_minus1, color=cmp, fill=cmp))+ 
  ggtitle(p1_tit, subtitle = p1_sub)+
  scale_y_continuous(breaks = y_range, labels = y_labs, name=p1_ylab)+
  scale_x_discrete(name="")+
  scale_color_manual(values=all_cols, labels=cmp_levels)+
  scale_fill_manual(values=all_cols, labels=cmp_levels)+
  labs(fill=p1_legTitle) +
  my_box_theme+
  geom_bar(stat="identity")+
  guides(color = FALSE)+
  coord_cartesian(clip = 'off', expand=F)  + 
  theme(
    text = element_text(family=fontFamily),
    axis.text.x = element_blank(),
        plot.margin = unit(plotMargin, "lines"))


p2_fccFract_plot <- ggplot(dt2, aes(x=dataset, y = ratioFCC, color = intervalFCC, fill = intervalFCC))+ 
  ggtitle(p2_tit, subtitle = p2_sub)+
  geom_bar(stat="identity") + 
  scale_fill_manual(values=mycols)+
  scale_color_manual(values=mycols) +
  # eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")")))+
  # eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")")))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), name=p2_ylab)+
  scale_x_discrete(labels=rep(labsymbol, length(dt2$dataset)), name=p2_xlab)+
  labs(fill=p2_legTitle) +
  my_box_theme+
  guides(color = FALSE)+
  coord_cartesian(clip = 'off', expand=F) + 
  theme(  
    text = element_text(family=fontFamily),
    axis.text.x = element_text(size=10, hjust=0.5, vjust=1, color=dotcols),
        plot.margin = unit(plotMargin, "lines"))


out_p <- p1_aucRatio_plot/p2_fccFract_plot + 
  plot_layout(heights = c(2, 1))

outFile <- file.path(outFolder, paste0("FCC_aucRatio_topFract.", plotType))
ggsave(plot = out_p, filename = outFile, height=myHeightGG, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))


write.table(dt1, file=file.path(outFolder, "fcc_auc_dt.txt"), col.names=T, row.names=F, sep="\t", quote=F)
write.table(dt2, file=file.path(outFolder, "fcc_range_dt.txt"), col.names=T, row.names=F, sep="\t", quote=F)

# library(gtable)
# g1 <- ggplotGrob(p1_aucRatio_plot)
# g2 <- ggplotGrob(p2_fccFract_plot)
# g <- rbind(g1, g2, size = "first")
# g$widths <- unit.pmax(g2$widths, g3$widths)
# grid.newpage()
# grid.draw(g)

# require(gridExtra)
# grid.arrange(p1_aucRatio_plot, p2_fccFract_plot, nrow=2, widths=c(1), heights=c(1,0.5))



