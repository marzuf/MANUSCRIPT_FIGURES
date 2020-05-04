# Rscript double_barplot_v3.R

require(ggplot2)
require(ggsci)
require(patchwork)
require(foreach)

plotType <- "svg"

source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")

outFolder <- "DOUBLE_BARPLOT_V3"
dir.create(outFolder, recursive = TRUE)


myHeight <- 5
myWidth <- 7

plotCex <- 1.4

interval_fcc <- c("<=0.5", "]0.5, 1[", "1")
p2_intervals <- c("]0.5, 1[", "1")

dt2 <- get(load("BARPLOT_WITH_TOPFCC_FRACT_V2//all_dt.Rdata"))
dt2$dataset <- file.path(dt2$hicds, dt2$exprds)
dt2$intervalFCC <- factor(dt2$intervalFCC, levels=interval_fcc)
stopifnot(!is.na(dt2$intervalFCC))
stopifnot(!duplicated(dt2))

dt1 <- get(load("FCC_WAVE_PLOT_NOABS/all_fcc_dt.Rdata"))
dt1$dataset <- file.path(dt1$hicds, dt1$exprds)
dt1 <- dt1[order(dt1$fcc_auc, decreasing = TRUE),]
ds_levels <- as.character(dt1$dataset)
stopifnot(!duplicated(dt1))

dt2$dataset <- factor(dt2$dataset, levels=ds_levels)
stopifnot(!is.na(dt2$dataset))

fcc_fract = interval_fcc[1]
col_dt <- foreach(fcc_fract = interval_fcc, .combine='rbind') %dopar% {
  sub_dt <- dt2[as.character(dt2$intervalFCC) == fcc_fract,]
  sub_dt <- sub_dt[as.numeric(sub_dt$dataset),]
  stopifnot(!is.na(sub_dt))
  sub_dt$cols <- colorRampPalette(c("blue","white", "red" ))(length(sub_dt$ratioFCC))[rank(sub_dt$ratioFCC)]
  sub_dt$ratioFCC_scaled <- as.numeric(scale(sub_dt$ratioFCC))
  sub_dt
}

dt2 <- merge(dt2, col_dt, by=colnames(col_dt)[!colnames(col_dt) %in% c("ratioFCC_scaled", "cols")], all=T)
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


p2_tit <- ""
p2_tit <- "FCC distribution"

p2_sub <- ""
p2_xlab <- ""
p2_ylab <- "FCC range"

p3_tit <- "FCC distribution - enrichment across datasets"
p3_tit <- ""
p3_tit <- "FCC distribution (by range)"
p3_sub <- ""
p3_xlab <- ""
p3_ylab <- "FCC range"

p1_legTitle <- ""
p2_legTitle <- "Ratio of TADs"
p3_legTitle <- "z-score\nratio of TADs\n(dataset resc.)"

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
  scale_color_manual(values=all_cols)+
  scale_fill_manual(values=all_cols)+
  labs(fill=p1_legTitle) +
  my_box_theme+
  geom_bar(stat="identity")+
  guides(color = FALSE)+
  coord_cartesian(clip = 'off', expand=F)  + 
  theme(axis.text.x = element_blank(),
        plot.margin = unit(plotMargin, "lines"))


# p2_fccFract_plot <- ggplot(dt2, aes(x=dataset, y = intervalFCC,  fill = ratioFCC))+ 
#   ggtitle(p2_tit, subtitle = p2_sub)+
#   geom_tile(colour="white") + 
#   scale_x_discrete(labels=rep(labsymbol, length(dt2$dataset)), name=p2_xlab)+
#   scale_fill_gradient2( high="red", low="blue", na.value = "grey", mid ="white", midpoint=mean(plot_dt$density_y))  +
#   labs(fill=p2_legTitle) +
#   my_box_theme+
#   guides(color = FALSE)+
#   coord_cartesian(clip = 'off', expand=F) +
#   theme(  
#     axis.text.x = element_text(size=10, hjust=0.5, vjust=1, color=dotcols),
#     plot.margin = unit(plotMargin, "lines")) 
#   # coord_flip()

p2_dt <- dt2[as.character(dt2$intervalFCC) %in% p2_intervals,]
p2_fccFract_plot <- ggplot(p2_dt,aes(x=dataset, y = intervalFCC,  fill = ratioFCC))+ 
  ggtitle(p2_tit, subtitle = p2_sub)+
  geom_tile(colour="white") + 
  scale_x_discrete(labels=rep(labsymbol, length(dt2$dataset)), name=p2_xlab)+
  scale_y_discrete(name=p2_ylab)+
  scale_fill_gradient2( high="red", low="blue", na.value = "grey", mid ="white", midpoint=mean(p2_dt$ratioFCC))  +
  labs(fill=p2_legTitle) +
  my_box_theme+
  guides(color = FALSE)+
  coord_cartesian(clip = 'off', expand=F) +
  theme(  
    text = element_text(family=fontFamily),
    axis.text.x = element_text(size=10, hjust=0.5, vjust=1, color=dotcols),
    plot.margin = unit(plotMargin, "lines")) 
# coord_flip()
save(dt2, file="dt2_resc.Rdata", version=2)

p3_fccFractResc_plot <- ggplot(dt2, aes(x=dataset, y = intervalFCC,  fill = ratioFCC_scaled))+ 
  ggtitle(p3_tit, subtitle = p2_sub)+
  geom_tile(colour="white") + 
  scale_x_discrete(labels=rep(labsymbol, length(dt2$dataset)), name=p3_xlab)+
  scale_y_discrete(name=p2_ylab)+
  scale_fill_gradient2( high="red", low="blue", na.value = "grey", mid ="white", midpoint=mean(dt2$ratioFCC_scaled))  +
  labs(fill=p3_legTitle) +
  my_box_theme+
  guides(color = FALSE)+
  coord_cartesian(clip = 'off', expand=F) +
  theme(  
    text = element_text(family=fontFamily),
    axis.text.x = element_text(size=10, hjust=0.5, vjust=1, color=dotcols),
    plot.margin = unit(plotMargin, "lines")) 
# coord_flip()


out_p <- p1_aucRatio_plot/p2_fccFract_plot/p3_fccFractResc_plot +
  plot_layout(heights = c(2, 1,1))


outFile <- file.path(outFolder, paste0("FCC_aucRatio_topFract.", plotType))
ggsave(plot = out_p, filename = outFile, height=myHeightGG, width=myWidth)
cat(paste0("... written: ", outFile, "\n"))
