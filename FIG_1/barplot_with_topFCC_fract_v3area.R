
# Rscript barplot_with_topFCC_fract_v3area.R

require(doMC)
require(foreach)
require(ggplot2)
require(ggsci)

require(ggpubr)

registerDoMC(40)

plotType <- "svg"

source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")


outFolder <- "BARPLOT_WITH_TOPFCC_FRACT_V3AREA"
dir.create(outFolder, recursive = TRUE)

buildData <- TRUE

minRange1 <- 0  # <= 0; ]0-0.5] ; ]0.5, 1[ ; 1
minRange2 <- 0.5
maxRange <- 1

fcc_fract_names <- c(paste0("<=", minRange1), 
					paste0("]", minRange1, ", ", minRange2, "]"),
				paste0("]", minRange2, ", ", maxRange, "["),
				 paste0(maxRange))

fract_sort <- fcc_fract_names[length(fcc_fract_names)]

ggsci_pal <- "lancet"
ggsci_subpal <- ""

legTitle <- "FCC ranges:"
fractBarSubTitle <- "AUC ratios:\n"
fractBarTitle <- "Fold-change concordance scores"

plotMargin <- c(1,2,1,1)

auc_ratio_file <- file.path("FCC_WAVE_PLOT_NOABS/all_fcc_dt.Rdata")
stopifnot(file.exists(auc_ratio_file))
x <- get(load(auc_ratio_file))
x$dataset <- paste0(x$hicds, "\n", x$exprds)
x <- x[order(x$fcc_auc, decreasing=TRUE),]
fcc_ds_order <- x$dataset
fcc_auc_sort <- TRUE
x$ds_rank <- 1:nrow(x)

signif_file <- file.path(runFolder, "CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(signif_file))

all_hicds <- all_obs_hicds
all_exprds <- all_obs_exprds

if(buildData){
  all_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
    cat(paste0("... start: ", hicds, "\n"))
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      cat(paste0("... start: ", hicds," - ", exprds,  "\n"))
      fcc_file <- file.path(pipFolder, hicds, exprds, step8fcc_folder, "all_obs_prodSignedRatio.Rdata")
      stopifnot(file.exists(fcc_file))
      all_fcc <- get(load(fcc_file))
      
      nRange0 <- sum(all_fcc <= minRange1)
      nRange1 <- sum(all_fcc > minRange1 & all_fcc <= minRange2)
      nRange2 <- sum(all_fcc > minRange2 & all_fcc < maxRange)
      nRange3 <- sum(all_fcc == 1)
      
stopifnot(nRange0+nRange1+nRange2+nRange3 == length(all_fcc))
	  
	  # data.frame(
	  #   hicds = hicds,
	  #   exprds = exprds,
	  #   FCC = all_fcc[all_fcc > 0.75],
	  #   stringsAsFactors = FALSE
	  # )
      data.frame(
        hicds = hicds,
        exprds = exprds,
        intervalFCC = fcc_fract_names,
        countFCC = c(nRange0, nRange1, nRange2, nRange3),
        ratioFCC = c(nRange0/length(all_fcc), nRange1/length(all_fcc), nRange2/length(all_fcc), nRange3/length(all_fcc)),
        stringsAsFactors = FALSE
      )
    }
    exprds_dt
  }
  outFile <- file.path(outFolder, paste0("all_dt.Rdata"))
  auc_fract_file <- outFile
  save(all_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  inFile <- file.path(outFolder, paste0("all_dt.Rdata"))
  auc_fract_file <- inFile
  all_dt <- get(load(inFile))
  # load("BARPLOT_WITH_FCC_FRACT/all_dt.Rdata")
}
# stop("-ok\n")
all_dt$intervalFCC <- factor(all_dt$intervalFCC, levels = rev(fcc_fract_names))
stopifnot(!is.na(all_dt$intervalFCC))

all_dt$dataset <- paste0(all_dt$hicds, "\n", all_dt$exprds)
all_dt <- all_dt[order(all_dt$countFCC, decreasing = TRUE),]
ds_order <- all_dt$dataset[all_dt$intervalFCC == fract_sort]
exprds_order <- all_dt$exprds[all_dt$intervalFCC == fract_sort]
exprds_order_fract <- exprds_order

if(fcc_auc_sort) {
  all_dt$dataset <- factor(all_dt$dataset, levels=fcc_ds_order)
  stopifnot(!is.na(all_dt$dataset))
  
} else {
  all_dt$dataset <- factor(all_dt$dataset, levels=ds_order)
  stopifnot(!is.na(all_dt$dataset))
}
tmp1 <- all_dt
tmp1 <- tmp1[order(as.numeric(tmp1$dataset)),]
exprds_order <- as.character(tmp1$exprds[tmp1$intervalFCC == fract_sort])

nDS <- length(unique(all_dt$dataset))

myds <- as.character(unique(file.path(all_dt$hicds, all_dt$exprds)))
countCmp <- setNames(as.numeric(table(all_cmps[basename(myds) ])), names(table(all_cmps[basename(myds) ])))

legDT <- data.frame(cmpType = names(all_cols), cmpColor = as.character(all_cols))
legDT <- legDT[rev(order(legDT$cmpType)),]
legDT$count <- countCmp[legDT$cmpType]
legDT$legLabel <- paste0(legDT$cmpType, " (", legDT$count, ")")

mycols <- all_cols[all_cmps[exprds_order]]
# tmp_dt <- aggregate(countFCC~ intervalFCC, data=all_dt, FUN=sum) # none has zero

mycols_fract <- all_cols[all_cmps[exprds_order_fract]]

all_dt$ds_rank <- as.numeric(all_dt$dataset)

fract_plot_with_lab <- ggplot(all_dt, aes(x=ds_rank, y=ratioFCC, fill=intervalFCC, color=intervalFCC)) + 
  geom_area()+
  ggtitle(paste0(fractBarTitle), 
          subtitle = paste0(fractBarSubTitle) )+
          # subtitle = "(all datasets)")+
  scale_x_continuous(name=paste0("(all datasets - n=", nDS, ")"))+
  labs(fill=paste0(legTitle)) + 
  guides(color=FALSE)+
  # coord_cartesian(expand=FALSE)+
  coord_cartesian(clip = 'off', expand=FALSE) +
  scale_y_continuous(name=paste0("Fraction of TADs"),
                     limits = c(0,1), 
                     breaks = seq(from=0, to=1, by=0.1),
                     labels = seq(from=0, to=1, by=0.1))+
  eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")")))+
  eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")")))+
  theme( # Increase size of axis lines
    plot.margin = unit(plotMargin, "lines"), # top, right, bottom, and left 
    plot.title = element_text(hjust = 0.5, face = "bold", size=16, family=fontFamily),
    plot.subtitle = element_text(hjust = 0, vjust=2,face = "italic", size =8, family=fontFamily, lineheight = 1.75),
    panel.grid = element_blank(),
    # panel.grid.major.y = element_line(colour = "grey"),
    # panel.grid.minor.y = element_line(colour = "grey"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x = element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .2, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12, family=fontFamily),
    axis.text.x = element_text(color=mycols, hjust=1,vjust = 0.5, size=7, angle=90, family=fontFamily),
    axis.title.y = element_text(color="black", size=14, family=fontFamily),
    axis.title.x = element_text(color="black", size=14, family=fontFamily),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.key = element_blank(),
    legend.title = element_text(face="bold")
  ) +
  geom_text(data=x, aes(x = x$ds_rank, y=1, 
                        label=sprintf("%.2f", x$fcc_auc)),
            inherit.aes=FALSE, angle=90, size=3, 
            vjust=0.5, hjust=0)+
  theme(
    # legend.position = c(.95, .95),
    # legend.box.just = "right",
    # legend.margin = margin(6, 6, 6, 6),
    legend.justification = c("right", "top")
  )+ 
  geom_text(data=legDT, aes(label = legDT$legLabel, x = 59, y =c(0, 0.05, 0.1)),
            vjust = 0, hjust=0,
            inherit.aes = FALSE, color = legDT$cmpColor)



outFile <- file.path(outFolder, paste0("all_ds_fcc_fract_scores_withLabs_areaplot.", plotType))
ggsave(plot = fract_plot_with_lab, filename = outFile, height=myHeightGG*1.5, width = myWidthGG*2)
cat(paste0("... written: ", outFile, "\n"))

all_dt$labSymb <- labsymbol

fract_plot_with_symb <- ggplot(all_dt, aes(x=ds_rank, y=ratioFCC, fill=intervalFCC, color=intervalFCC)) + 
  geom_area()+
  # coord_cartesian(expand = FALSE) +
  coord_cartesian(clip = 'off', expand=FALSE) +
  ggtitle(paste0(fractBarTitle), 
          subtitle = "AUC ratios:\n")+
          # subtitle = "(all datasets)")+
  labs(fill=legTitle)+
  guides(color=FALSE)+
  eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")")))+
  eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")")))+
  scale_y_continuous(name=paste0("Fraction of TADs"),
                     limits = c(0,1), 
                     breaks = seq(from=0, to=1, by=0.1),
                     labels = seq(from=0, to=1, by=0.1))+
  scale_x_continuous(labels=all_dt$labSymb, breaks =1:nrow(all_dt), name=paste0("(all datasets - n=", nDS, ")"))+
  theme( # Increase size of axis lines
    plot.margin = unit(plotMargin, "lines"), # top, right, bottom, and left 
    plot.title = element_text(hjust = 0.5, face = "bold", size=16, family=fontFamily),
    plot.subtitle = element_text(hjust = 0, vjust=2, face = "italic", size = 8, family=fontFamily, lineheight = 1.75),
    panel.grid = element_blank(),
    # panel.grid.major.y = element_line(colour = "grey"),
    # panel.grid.minor.y = element_line(colour = "grey"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x = element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .2, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12, family=fontFamily),
    axis.text.x = element_text(color=mycols, hjust=0.5,vjust = 0.5, size=12, family=fontFamily),
    axis.title.y = element_text(color="black", size=14, family=fontFamily),
    axis.title.x = element_text(color="black", size=14, family=fontFamily),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.key = element_blank(),
    legend.title = element_text(face="bold", family=fontFamily)
  )+
  geom_text(data=x, aes(x = x$ds_rank, y=1, 
                        label=sprintf("%.2f", x$fcc_auc)),
            inherit.aes=FALSE, angle=90, size=3, vjust=0.5, hjust=0) +
  theme(
    # legend.position = c(.95, .95),
    # legend.box.just = "right",
    # legend.margin = margin(6, 6, 6, 6),
    legend.justification = c("right", "top")
  ) + 
  geom_text(data=legDT, aes(label = paste0(labsymbol, " ", legDT$legLabel), x = 59, y =c(0, 0.05, 0.1)),
            vjust = 0, hjust=0,
            inherit.aes = FALSE, color = legDT$cmpColor)

outFile <- file.path(outFolder, paste0("all_ds_fcc_fract_scores_withSymb_areaplot.", plotType))
ggsave(plot = fract_plot_with_symb, filename = outFile, height=myHeightGG, width = myWidthGG*2)
cat(paste0("... written: ", outFile, "\n"))

save(all_dt, file="all_dt.Rdata", version=2)

tmp <- all_dt[as.character(all_dt$intervalFCC) == fract_sort,]
tmp <- tmp[order(tmp$countFCC, decreasing = TRUE),]

all_dt$dataset <- factor(all_dt$dataset, levels=tmp$dataset)
all_dt$ds_rank <- as.numeric(all_dt$dataset)

fract_plot_with_symb <- ggplot(all_dt, aes(x=ds_rank, y=ratioFCC, fill=intervalFCC, color=intervalFCC)) + 
  # coord_cartesian(expand = FALSE) +
  coord_cartesian(clip = 'off', expand=FALSE) +
  ggtitle(paste0(fractBarTitle), 
          subtitle = "AUC ratios:\n")+
  # subtitle = "(all datasets)")+
  labs(fill=legTitle)+
  guides(color=FALSE)+
  eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")")))+
  eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")")))+
  scale_y_continuous(name=paste0("Fraction of TADs"),
                     limits = c(0,1), 
                     breaks = seq(from=0, to=1, by=0.1),
                     labels = seq(from=0, to=1, by=0.1))+
  scale_x_continuous(labels=all_dt$labSymb, breaks = 1:nrow(all_dt), name=paste0("(all datasets - n=", nDS, ")"))+
  theme( # Increase size of axis lines
    plot.margin = unit(plotMargin, "lines"), # top, right, bottom, and left 
    plot.title = element_text(hjust = 0.5, face = "bold", size=16, family=fontFamily),
    plot.subtitle = element_text(hjust = 0, vjust=2, face = "italic", size = 8, family=fontFamily, lineheight = 1.75),
    panel.grid = element_blank(),
    # panel.grid.major.y = element_line(colour = "grey"),
    # panel.grid.minor.y = element_line(colour = "grey"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x = element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .2, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12, family=fontFamily),
    axis.text.x = element_text(color=mycols_fract, hjust=0.5,vjust = 0.5, size=12, family=fontFamily),
    axis.title.y = element_text(color="black", size=14, family=fontFamily),
    axis.title.x = element_text(color="black", size=14, family=fontFamily),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.key = element_blank(),
    legend.title = element_text(face="bold", family=fontFamily)
  )+
  geom_text(data=x, aes(x = x$ds_rank, y=1, 
                        label=sprintf("%.2f", x$fcc_auc)),
            inherit.aes=FALSE, angle=90, size=3, vjust=0.5, hjust=0) +
  theme(
    # legend.position = c(.95, .95),
    # legend.box.just = "right",
    # legend.margin = margin(6, 6, 6, 6),
    legend.justification = c("right", "top")
  ) + 
  geom_text(data=legDT, aes(label = paste0(labsymbol, " ", legDT$legLabel), x = 59, y =c(0, 0.05, 0.1)),
            vjust = 0, hjust=0,
            inherit.aes = FALSE, color = legDT$cmpColor)



all_dt$ds_rank <- as.numeric(all_dt$dataset)



outFile <- file.path(outFolder, paste0("all_ds_fcc_fract_scores_withSymb_fractSorted_areaplot.", plotType))
ggsave(plot = fract_plot_with_symb, filename = outFile, height=myHeightGG, width = myWidthGG*2)
cat(paste0("... written: ", outFile, "\n"))



fract_plot_with_symb <- ggplot(all_dt, aes(x=ds_rank, y=ratioFCC, fill=intervalFCC, color=intervalFCC)) + 
  geom_area()+
  # coord_cartesian(expand = FALSE) +
  coord_cartesian(clip = 'off', expand=FALSE) +
  ggtitle(paste0(fractBarTitle), 
          subtitle = "AUC ratios:\n")+
  # subtitle = "(all datasets)")+
  labs(fill=legTitle)+
  guides(color=FALSE)+
  eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")")))+
  eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")")))+
  scale_y_continuous(name=paste0("Fraction of TADs"),
                     limits = c(0,1), 
                     breaks = seq(from=0, to=1, by=0.1),
                     labels = seq(from=0, to=1, by=0.1))+
  scale_x_continuous(name=paste0("(all datasets - n=", nDS, ")"))+
  theme( # Increase size of axis lines
    plot.margin = unit(plotMargin, "lines"), # top, right, bottom, and left 
    plot.title = element_text(hjust = 0.5, face = "bold", size=16, family=fontFamily),
    plot.subtitle = element_text(hjust = 0, vjust=2, face = "italic", size = 8, family=fontFamily, lineheight = 1.75),
    panel.grid = element_blank(),
    # panel.grid.major.y = element_line(colour = "grey"),
    # panel.grid.minor.y = element_line(colour = "grey"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x = element_line(size = .2, color = "black"),
    axis.line.y = element_line(size = .2, color = "black"),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12, family=fontFamily),
    axis.text.x = element_text(color=mycols_fract, hjust=1,vjust = 0.5, size=7, family=fontFamily, angle=90),
    axis.title.y = element_text(color="black", size=14, family=fontFamily),
    axis.title.x = element_text(color="black", size=14, family=fontFamily),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.key = element_blank(),
    legend.title = element_text(face="bold", family=fontFamily)
  )+
  geom_text(data=x, aes(x = x$ds_rank, y=1, 
                        label=sprintf("%.2f", x$fcc_auc)),
            inherit.aes=FALSE, angle=90, size=3, vjust=0.5, hjust=0) +
  theme(
    # legend.position = c(.95, .95),
    # legend.box.just = "right",
    # legend.margin = margin(6, 6, 6, 6),
    legend.justification = c("right", "top")
  ) + 
  geom_text(data=legDT, aes(label = paste0(labsymbol, " ", legDT$legLabel), x = 59, y =c(0, 0.05, 0.1)),
            vjust = 0, hjust=0,
            inherit.aes = FALSE, color = legDT$cmpColor)



outFile <- file.path(outFolder, paste0("all_ds_fcc_fract_scores_withSymb_withLabs_fractSorted_areaplot.", plotType))
ggsave(plot = fract_plot_with_symb, filename = outFile, height=myHeightGG*1.5, width = myWidthGG*2)
cat(paste0("... written: ", outFile, "\n"))



