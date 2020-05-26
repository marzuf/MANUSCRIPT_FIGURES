# Rscript fcc_ratioNbrTADs_pychart_obsPerm.R 

require(doMC)
require(foreach)
require(ggplot2)
require(ggsci)
require(colorRamps)
require(ggcharts)
require(reshape2)
require(ggpubr)

registerDoMC(40)

plotType <- "svg"

source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")

a_t <- "permg2t"

outFolder <- "FCC_RATIONBRTADS_PYCHART_OBSPERM"
dir.create(outFolder, recursive = TRUE)

buildData <- FALSE

fcc_fract <- seq(from=-1, to=1, by=0.1)
# fcc_fract_names <- paste0("FCC > ", fcc_fract[1:(length(fcc_fract)-1)], " and FCC <= ",fcc_fract[2:length(fcc_fract)])
fcc_fract_names <- paste0("FCC \u2208 ]", fcc_fract[1:(length(fcc_fract)-1)], ", ",fcc_fract[2:length(fcc_fract)], "]")
fcc_fract_names <- paste0("]", fcc_fract[1:(length(fcc_fract)-1)], ", ",fcc_fract[2:length(fcc_fract)], "]")

fcc_fract_names[fcc_fract_names == "]-1, -0.9]"] <- "[-1, -0.9]"


# fract_sort <- "FCC > 0.75 and FCC <= 1"
fract_sort <- fcc_fract_names[length(fcc_fract_names)]

ggsci_pal <- "lancet"
ggsci_subpal <- ""

legTitle <- "FCC ranges:"
fractBarSubTitle <- "AUC ratios:\n"
fractBarTitle <- "Fold-change concordance scores"

plotMargin <- c(1,2,1,1)

plotType <- "svg"

myWidthGG <- 9
myHeightGG <- 7
ontFamily <- "Hershey"

auc_ratio_file <- file.path("FCC_WAVE_PLOT_NOABS/all_fcc_dt.Rdata")
stopifnot(file.exists(auc_ratio_file))
x <- get(load(auc_ratio_file))
x$dataset <- file.path(x$hicds,  x$exprds)
x <- x[order(x$fcc_auc, decreasing=TRUE),]
fcc_ds_order <- x$dataset

dsCols <- all_cols[all_cmps[basename(fcc_ds_order)]]
stopifnot(!is.na(dsCols))

signif_file <- file.path(runFolder, "CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(signif_file))

all_hicds <- all_obs_hicds
all_exprds <- all_obs_exprds

# all_hicds=all_hicds[1]

keepPermut <- 1000

obs_dt <- get(load(file.path("FCC_RATIONBRTADS_HEATMAP/all_dt.Rdata")))
# colnames(obs_dt)[4:5] <- paste0(colnames(obs_dt)[4:5], "_obs")
perm_dt <- get(load(file.path("FCC_RATIONBRTADS_HEATMAP_PERMG2T//all_dt.Rdata")))
# colnames(perm_dt)[4:5] <- paste0(colnames(perm_dt)[4:5], "_perm")

obsPerm_dt <- merge(obs_dt, perm_dt, by=c("hicds", "exprds", "intervalFCC"), suffixes = c("_obs", "_perm"), all=TRUE)
stopifnot(!is.na(obsPerm_dt))

obsPerm_dt$intervalFCC <- factor(obsPerm_dt$intervalFCC, levels=fcc_fract_names)


obsPerm_dt$nFCC_perm_mean <- obsPerm_dt$nFCC_perm/keepPermut


obsPerm_dt$ratioFCC_obs_perm <- obsPerm_dt$ratioFCC_obs/obsPerm_dt$ratioFCC_perm

obsPerm_dt$nFCC_obs_perm <- obsPerm_dt$nFCC_obs/obsPerm_dt$nFCC_perm_mean

stopifnot(round(na.omit(obsPerm_dt$nFCC_obs_perm),4) ==  round(na.omit(obsPerm_dt$ratioFCC_obs_perm),4))

obsPerm_dt$dataset <- file.path(obsPerm_dt$hicds, obsPerm_dt$exprds)
nDS <- length(unique(obsPerm_dt$dataset))

obsPerm_dt$ratioFCC_obs_perm_log2 <- log2(obsPerm_dt$ratioFCC_obs_perm)

obsPerm_dt$nFCC_obs_perm_log2 <- log2(obsPerm_dt$nFCC_obs_perm)

cut_dt <- obsPerm_dt[as.character(obsPerm_dt$intervalFCC) %in% unique(as.character(obsPerm_dt$intervalFCC[obsPerm_dt$ratioFCC_obs > 0])),]

cut_dt$dataset <- file.path(cut_dt$hicds, cut_dt$exprds)


plot_cut_dt <- cut_dt[,c("dataset","intervalFCC", "ratioFCC_obs", "ratioFCC_perm")]

agg1 <- aggregate(ratioFCC_obs~intervalFCC, data =plot_cut_dt, FUN=mean)
agg2 <- aggregate(ratioFCC_perm~intervalFCC, data =plot_cut_dt, FUN=mean)
plot_cut_dt2 <- merge(agg1, agg2, by="intervalFCC", all=T)
m_plot_dt <- melt(plot_cut_dt2, id=c("intervalFCC"))
m_plot_dt$intervalFCC <- factor(m_plot_dt$intervalFCC, levels=fcc_fract_names)
stopifnot(!is.na(m_plot_dt))

m_plot_dt$variable <- as.character(m_plot_dt$variable)

var1 <- "obs. data"
var2 <- paste0("perm. data (mean", keepPermut, " permut)")

m_plot_dt$variable[as.character(m_plot_dt$variable) == "ratioFCC_obs"] <- var1
m_plot_dt$variable[as.character(m_plot_dt$variable) == "ratioFCC_perm"] <- var2

m_plot_dt <- m_plot_dt[order(as.numeric(m_plot_dt$intervalFCC)),]
  
# pc_plot <- pyramid_chart(data = m_plot_dt, x = intervalFCC, y = value, group = variable, xlab="Ratio of TADs")


bar_colors <-  c("steelblue3", "orangered")
plot_limit <- abs(max(m_plot_dt$value))
plots <- list()
sides <- c("left", "right")

mytit1 <- "Observed data."
mytit2 <- "Permut. data."

horiz1 <- "Ratio of TADs (mean all datasets)"
horiz2 <- "Ratio of TADs(mean all permut mean)"
horiz1 <- "(mean all datasets)"
horiz2 <- paste0("(mean ", keepPermut, " permut mean)")

plotTit <- "Distribution TAD FCC"
subTit <- "Ratio of TADs"

# draw the first plot

y_scale1<- scale_y_reverse(limits = c(plot_limit, 0), name=paste0(horiz1),  # horizontal axis
                           expand = expand_scale(mult = c(0.05, 0)))

plots[[1]]<- ggplot(m_plot_dt[as.character(m_plot_dt$variable) == var1, ], aes(x = intervalFCC, y = value))+ 
  geom_col(fill = bar_colors[1], width = 0.7) + 
  y_scale1 + coord_flip() + ggcharts:::pyramid_theme(sides[1]) + 
  scale_x_discrete(expand = expand_scale(add = 0.5)) + 
  theme(
    axis.title.x=element_text(),
    axis.text.x=element_text(size=14),
    axis.line.x = element_line()
    )+
  ggtitle(mytit1)


y_scale2 <- scale_y_continuous(limits = c(0, plot_limit), name=paste0(horiz2),
                               expand = expand_scale(mult = c(0, 0.05)))

plots[[2]]<- ggplot(m_plot_dt[as.character(m_plot_dt$variable) == var2, ], aes(x = intervalFCC, y = value))+ 
  geom_col(fill = bar_colors[2], width = 0.7) + scale_x_discrete(expand = expand_scale(add = 0.5)) + 
  y_scale2 + 
  labs(xlab="")+
  theme(axis.title.x=element_text())+
  coord_flip() + ggcharts:::pyramid_theme(sides[2]) + 
  theme(
    axis.title.x=element_text(),
    axis.text.x=element_text(size=14),
    axis.line.x = element_line()
  )+
  ggtitle(mytit2)

pc_plot <- plots[[1]] + plots[[2]] + patchwork::plot_annotation(caption = subTit, 
                                                     title = plotTit, 
                                                     theme = theme(plot.caption = element_text(hjust = 0.5, face="bold", family=fontFamily,
                                                                                                              size = 16),
                                                                   plot.title = element_text(hjust = 0.5, face="bold",family=fontFamily,
                                                                                               size = 18)
                                                                   ))


outFile <- file.path(outFolder, paste0("FCC_score_mean_pyramidchart.", plotType))
ggsave(pc_plot, filename = outFile,  height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


