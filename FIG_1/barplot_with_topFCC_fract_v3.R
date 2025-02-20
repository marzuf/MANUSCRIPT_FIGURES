
# Rscript barplot_with_topFCC_fract_v3.R

require(doMC)
require(foreach)
require(ggplot2)
require(ggsci)

require(ggpubr)

registerDoMC(40)

plotType <- "svg"

source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")


outFolder <- "BARPLOT_WITH_TOPFCC_FRACT_V3"
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

fract_plot_with_lab <- ggplot(all_dt, aes(x=dataset, y=ratioFCC, fill=intervalFCC, color=intervalFCC)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle(paste0(fractBarTitle), 
          subtitle = paste0(fractBarSubTitle) )+
          # subtitle = "(all datasets)")+
  scale_x_discrete(name=paste0("(all datasets - n=", nDS, ")"))+
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
	text = element_text(family=fontFamily),
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
  geom_text(data=x, aes(x = x$dataset, y=1, 
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



outFile <- file.path(outFolder, paste0("all_ds_fcc_fract_scores_withLabs_barplot.", plotType))
ggsave(plot = fract_plot_with_lab, filename = outFile, height=myHeightGG*1.5, width = myWidthGG*2)
cat(paste0("... written: ", outFile, "\n"))

all_dt$labSymb <- labsymbol

fract_plot_with_symb <- ggplot(all_dt, aes(x=dataset, y=ratioFCC, fill=intervalFCC, color=intervalFCC)) + 
  geom_bar(position="stack", stat="identity") +
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
  scale_x_discrete(labels=all_dt$labSymb, name=paste0("(all datasets - n=", nDS, ")"))+
  theme( # Increase size of axis lines
	text = element_text(family=fontFamily),
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
  geom_text(data=x, aes(x = x$dataset, y=1, 
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

outFile <- file.path(outFolder, paste0("all_ds_fcc_fract_scores_withSymb_barplot.", plotType))
ggsave(plot = fract_plot_with_symb, filename = outFile, height=myHeightGG, width = myWidthGG*2)
cat(paste0("... written: ", outFile, "\n"))

save(all_dt, file="all_dt.Rdata", version=2)

tmp <- all_dt[as.character(all_dt$intervalFCC) == fract_sort,]
tmp <- tmp[order(tmp$countFCC, decreasing = TRUE),]

all_dt$dataset <- factor(all_dt$dataset, levels=tmp$dataset)

fract_plot_with_symb <- ggplot(all_dt, aes(x=dataset, y=ratioFCC, fill=intervalFCC, color=intervalFCC)) + 
  geom_bar(position="stack", stat="identity") +
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
  scale_x_discrete(labels=all_dt$labSymb, name=paste0("(all datasets - n=", nDS, ")"))+
  theme( # Increase size of axis lines
	text = element_text(family=fontFamily),
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
  geom_text(data=x, aes(x = x$dataset, y=1, 
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



outFile <- file.path(outFolder, paste0("all_ds_fcc_fract_scores_withSymb_fractSorted_barplot.", plotType))
ggsave(plot = fract_plot_with_symb, filename = outFile, height=myHeightGG, width = myWidthGG*2)
cat(paste0("... written: ", outFile, "\n"))



fract_plot_with_symb <- ggplot(all_dt, aes(x=dataset, y=ratioFCC, fill=intervalFCC, color=intervalFCC)) + 
  geom_bar(position="stack", stat="identity") +
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
  scale_x_discrete(name=paste0("(all datasets - n=", nDS, ")"))+
  theme( # Increase size of axis lines
	text = element_text(family=fontFamily),
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
  geom_text(data=x, aes(x = x$dataset, y=1, 
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



outFile <- file.path(outFolder, paste0("all_ds_fcc_fract_scores_withSymb_withLabs_fractSorted_barplot.", plotType))
ggsave(plot = fract_plot_with_symb, filename = outFile, height=myHeightGG*1.5, width = myWidthGG*2)
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
# FRACT and FCC AUC - cont x-axis
######################################################################################

auc_fract_dt <- get(load(auc_fract_file))
auc_ratio_dt <- get(load(auc_ratio_file))

auc_fract_ratio_dt <- merge(auc_fract_dt, auc_ratio_dt, by=c("hicds", "exprds"), all.x=TRUE, all.y=TRUE)
stopifnot(!is.na(auc_fract_ratio_dt))

auc_fract_ratio_dt$intervalFCC <- gsub("FCC ", "", auc_fract_ratio_dt$intervalFCC)
# auc_fract_ratio_dt$intervalFCC <- factor(auc_fract_ratio_dt$intervalFCC, levels=rev(fcc_fract_names))
auc_fract_ratio_dt$intervalFCC <- factor(auc_fract_ratio_dt$intervalFCC, levels=rev(gsub("FCC ", "", fcc_fract_names)))

myPals <-  eval(parse(text=paste0("pal_", ggsci_pal, "(", ggsci_subpal, ")")))(length(unique(auc_fract_ratio_dt$intervalFCC)))

scatPlot <- ggscatter(auc_fract_ratio_dt, 
                      title = paste0("all datasets (n=", length(unique(file.path(auc_fract_ratio_dt$hicds, auc_fract_ratio_dt$exprds))), ")"),
                      x = "fcc_auc", 
                      y = "countFCC",
                      color = "intervalFCC",
                      xlab = "FCC AUC ratio",
                      ylab = "Ratio of TADs",
                      palette = myPals)+
  labs(color=legTitle)+
  geom_smooth(aes(color = intervalFCC),method = "lm", linetype=2, se=F)+
  theme( # Increase size of axis lines
	text = element_text(family=fontFamily),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16, family=fontFamily),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14, family=fontFamily),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12, family=fontFamily),
    axis.text.x = element_text(color="black", hjust=0.5,vjust = 1, size=12, family=fontFamily),
    axis.title.y = element_text(color="black", size=14, family=fontFamily),
    axis.title.x = element_text(color="black", size=14, family=fontFamily),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.title = element_text(face="bold", family=fontFamily)
  )+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))



outFile <- file.path(outFolder, paste0("all_ds_fcc_fract_scores_nbrSignifs_auc_scatterplot_xcont.", plotType))
ggsave(plot = scatPlot, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

######################################################################################
# FRACT and FCC AUC - discrete x-axis
######################################################################################

xax_exp <- 0.01

# auc_fract_dt <- get(load("BARPLOT_WITH_FCC_FRACT/all_dt.Rdata"))
# auc_fract_ratio_dt <- get(load("BARPLOT_FCC_AUC_RATIO//all_dt.Rdata"))

auc_fract_dt <- get(load(auc_fract_file))
auc_ratio_dt <- get(load(auc_ratio_file))

auc_fract_ratio_dt <- merge(auc_fract_dt, auc_ratio_dt, by=c("hicds", "exprds"), all.x=TRUE, all.y=TRUE)
stopifnot(!is.na(auc_fract_ratio_dt))

auc_fract_ratio_dt$intervalFCC <- gsub("FCC ", "", auc_fract_ratio_dt$intervalFCC)
# auc_fract_ratio_dt$intervalFCC <- factor(auc_fract_ratio_dt$intervalFCC, levels=rev(fcc_fract_names))
auc_fract_ratio_dt$intervalFCC <- factor(auc_fract_ratio_dt$intervalFCC, levels=rev(gsub("FCC ", "", fcc_fract_names)))

ds_rank_dt <- auc_fract_ratio_dt[,c("hicds", "exprds", "fcc_auc")]
ds_rank_dt$dataset <- paste0(ds_rank_dt$hicds,"\n", ds_rank_dt$exprds)
ds_rank_dt <- unique(ds_rank_dt)
stopifnot(!duplicated(ds_rank_dt$dataset))
ds_rank_dt <- ds_rank_dt[order(ds_rank_dt$fcc_auc, decreasing=TRUE),]
ds_auc_order <- ds_rank_dt$dataset
mycols_scat <- all_cols[all_cmps[paste0(ds_rank_dt$exprds)]]

auc_fract_ratio_dt$dataset <- paste0(auc_fract_ratio_dt$hicds, "\n", auc_fract_ratio_dt$exprds)
auc_fract_ratio_dt$dataset <- factor(auc_fract_ratio_dt$dataset, levels=ds_auc_order)
stopifnot(!is.na(auc_fract_ratio_dt$dataset))

auc_fract_ratio_dt$ds_rank <- as.numeric(auc_fract_ratio_dt$dataset)
xbreaks <- sort(unique(auc_fract_ratio_dt$ds_rank))

myPals <-  eval(parse(text=paste0("pal_", ggsci_pal, "(", ggsci_subpal, ")")))(length(unique(auc_fract_ratio_dt$intervalFCC)))

scatPlot <- ggscatter(auc_fract_ratio_dt, 
                      title = paste0("all datasets (n=", length(unique(file.path(auc_fract_ratio_dt$hicds, auc_fract_ratio_dt$exprds))), ")"),
                      x = "ds_rank", 
                      y = "countFCC",
                      color = "intervalFCC",
                      xlab = "Datasets (ranked by FCC AUC ratio)",
                      ylab = "Ratio of TADs",
                      palette = myPals) +
  labs(color=legTitle)+
  geom_smooth(aes(color = intervalFCC),method = "lm", linetype=2, se=F)+
  
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_x_continuous(labels = rep(labsymbol, length(xbreaks)), breaks= xbreaks, expand=c(xax_exp,xax_exp)) + 
  
  
  
  theme( # Increase size of axis lines
	text = element_text(family=fontFamily),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16, family=fontFamily),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14, family=fontFamily),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12, family=fontFamily),
    axis.text.x = element_text(color=mycols_scat, hjust=0.5,vjust = 0.5, size=12, family=fontFamily),
    axis.title.y = element_text(color="black", size=14, family=fontFamily),
    axis.title.x = element_text(color="black", size=14, family=fontFamily),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.title = element_text(face="bold", family=fontFamily)
  )

outFile <- file.path(outFolder, paste0("all_ds_fcc_fract_scores_nbrSignifs_auc_scatterplot_xcont.", plotType))
ggsave(plot = scatPlot, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))

scatPlot <- ggscatter(auc_fract_ratio_dt, 
                      title = paste0("all datasets (n=", length(unique(file.path(auc_fract_ratio_dt$hicds, auc_fract_ratio_dt$exprds))), ")"),
                      x = "ds_rank", 
                      y = "countFCC",
                      color = "intervalFCC",
                      xlab = "Datasets (ranked by FCC AUC ratio)",
                      ylab = "Ratio of TADs",
                      palette = myPals) +
  labs(color=legTitle)+
  geom_smooth(aes(color = intervalFCC),method = "lm", linetype=2, se=F)+
  
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_x_continuous(labels = ds_auc_order, breaks=xbreaks, expand=c(xax_exp,xax_exp))+

  theme( # Increase size of axis lines
	text = element_text(family=fontFamily),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16, family=fontFamily),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14, family=fontFamily),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12, family=fontFamily),
    axis.text.x = element_text(color=mycols_scat, hjust=1,vjust = 0.5, size=7, family=fontFamily, angle=90),
    axis.title.y = element_text(color="black", size=14, family=fontFamily),
    axis.title.x = element_text(color="black", size=14, family=fontFamily),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.title = element_text(face="bold", family=fontFamily)
  )

outFile <- file.path(outFolder, paste0("all_ds_fcc_fract_scores_nbrSignifs_auc_scatterplot_xcont_withLabs.", plotType))
ggsave(plot = scatPlot, filename = outFile, height=myHeightGG, width = myWidthGG*1.5)
cat(paste0("... written: ", outFile, "\n"))


######################################################################################
# FRACT and nSignif
######################################################################################



all_dt <- get(load(auc_fract_file))
auc_fract_dt <- all_dt
all_auc_signif_dt <- get(load(signif_file))

signifTADthresh <- 0.01
auc_signif_dt <- aggregate(adjPvalComb~hicds+exprds, data=all_auc_signif_dt, FUN=function(x) sum(x <= signifTADthresh))
colnames(auc_signif_dt)[colnames(auc_signif_dt) == "adjPvalComb"] <- "nSignifTADs"

auc_fract_signif_dt <- merge(auc_fract_dt, auc_signif_dt, by=c("hicds", "exprds"), all.x=TRUE, all.y=TRUE)
stopifnot(!is.na(auc_fract_signif_dt))

auc_fract_signif_dt$intervalFCC <- gsub("FCC ", "", auc_fract_signif_dt$intervalFCC)
# auc_fract_signif_dt$intervalFCC <- factor(auc_fract_signif_dt$intervalFCC, levels=rev(fcc_fract_names))
auc_fract_signif_dt$intervalFCC <- factor(auc_fract_signif_dt$intervalFCC, levels=rev(gsub("FCC ", "", fcc_fract_names)))

myPals <-  eval(parse(text=paste0("pal_", ggsci_pal, "(", ggsci_subpal, ")")))(length(unique(auc_fract_signif_dt$intervalFCC)))

scatPlot <- ggscatter(auc_fract_signif_dt, 
                      title = paste0("all datasets (n=", length(unique(file.path(auc_fract_signif_dt$hicds, auc_fract_signif_dt$exprds))), ")"),
                      x = "nSignifTADs", 
                      y = "countFCC",
                      color = "intervalFCC",
                      xlab = paste0("# signif. TADs (p-val <= ", signifTADthresh, ")"),
                      ylab = "Ratio of TADs",
                      palette = myPals)+
  labs(color=legTitle)+
  geom_smooth(aes(color = intervalFCC),method = "lm", linetype=2, se=F)+
  theme( # Increase size of axis lines
	text = element_text(family=fontFamily),
    plot.title = element_text(hjust = 0.5, face = "bold", size=16, family=fontFamily),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14, family=fontFamily),
    axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12, family=fontFamily),
    axis.text.x = element_text(color="black", hjust=0.5,vjust =1, size=12, family=fontFamily),
    axis.title.y = element_text(color="black", size=14, family=fontFamily),
    axis.title.x = element_text(color="black", size=14, family=fontFamily),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    legend.background =  element_rect(),
    legend.title = element_text(face="bold", family=fontFamily)
  )+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))


outFile <- file.path(outFolder, paste0("all_ds_fcc_fract_scores_fcc_auc_scatterplot_xcont.", plotType))
ggsave(plot = scatPlot, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))




