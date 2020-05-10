# Rscript n_fcc_over_qtile_permg2t_iterate2.R

outFolder <- "N_FCC_OVER_QTILE_PERMG2T_ITERATE2"
dir.create(outFolder)


plotType <- "svg"
myHeight <- 7
myWidth <- 10
plotCex <- 1.4

# source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
# source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")

require(foreach)
require(doMC)
registerDoMC(40)
require(ggplot2)
require(ggpubr)
require(reshape2)

buildData <- FALSE

all_fccThresh <- seq(from = -1, to = 1, by = 0.05)

myHeightGG <- 0.8*myHeightGG
myWidthGG <- 1.2*myWidthGG

keepPermut <- 1000


rd_type <- "PERMG2T"

buildData1 <- FALSE
buildData2 <- FALSE

all_hicds <- all_obs_hicds


if(buildData1){
  hicds = all_hicds[1]
  hicds = all_hicds[2]
  all_fcc_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    hicds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      cat("... start ", hicds, " - ", exprds, " \n")
      
      rd_file <- file.path(pipFolder, hicds, exprds, "8cOnlyFCC_runAllDown/prodSignedRatio_permDT.Rdata")
      rd_fcc_dt <- get(load(rd_file))
      rd_fcc_dt <- rd_fcc_dt[,1:keepPermut]
      rd_fcc <- as.numeric(rd_fcc_dt)
      
      fccFile <- file.path(pipFolder,  hicds, exprds, "8cOnlyFCC_runAllDown/all_obs_prodSignedRatio.Rdata")
      obs_fcc <- get(load(fccFile))
      
      dt1 <- data.frame(
        hicds=hicds,
        exprds=exprds,
        fcc_type="random",
        fcc_value = rd_fcc,
        stringsAsFactors = FALSE
      )
      dt2 <- data.frame(
        hicds=hicds,
        exprds=exprds,
        fcc_type="observed",
        fcc_value = as.numeric(obs_fcc),
        stringsAsFactors = FALSE
      )
      rbind(dt1,dt2)
    }
    hicds_dt
  }
  outFile <- file.path(outFolder, "all_fcc_dt.Rdata")
  save(all_fcc_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
  overall_min_fcc <- min(all_fcc_dt$fcc_value) # -0.61
  outFile <- file.path(outFolder, "overall_min_fcc.Rdata")
  save(overall_min_fcc, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  if(buildData2) {
    cat(paste0("... load all_fcc_dt\n"))
    outFile <- file.path(outFolder, "all_fcc_dt.Rdata")
    all_fcc_dt <- get(load(outFile))
    cat(paste0("... done\n"))  
    all_fcc_dt$dataset <- file.path(all_fcc_dt$hicds, all_fcc_dt$exprds)
    
  } else {
    outFile <- file.path(outFolder, "overall_min_fcc.Rdata")
    overall_min_fcc <- get(load(outFile))
    
  }
  
}

all_fccThresh <- all_fccThresh[all_fccThresh >= overall_min_fcc] 

if(buildData2){
  
  ds = unique(all_fcc_dt$dataset)[1]
  nOverThresh_dt <- foreach(fccThresh = all_fccThresh, .combine='rbind') %dopar% {
    
    ds_dt <- foreach(ds = unique(all_fcc_dt$dataset), .combine='rbind') %do% {
      
      cat("... start - ", fccThresh, " - ", ds, "\n")
      
      curr_dt <- all_fcc_dt[all_fcc_dt$dataset == ds,]
      
      nObs <- sum(curr_dt$fcc_type == "observed")
      nRd <- sum(curr_dt$fcc_type == "random")
      
      nOverThresh_obs <-   sum(curr_dt$fcc_value[curr_dt$fcc_type == "observed"] >= fccThresh )
      nOverThresh_rd <-   sum(curr_dt$fcc_value[curr_dt$fcc_type == "random"] >= fccThresh )
      
      data.frame(
        hicds=dirname(ds),
        exprds=basename(ds),
        fccThresh = fccThresh,
        nOverThresh_obs = nOverThresh_obs,
        nOverThresh_rd = nOverThresh_rd,
        ratioOverThresh_obs = nOverThresh_obs/nObs,
        ratioOverThresh_rd = nOverThresh_rd/nRd,
        stringsAsFactors = FALSE
      )
    }
    ds_dt
  }
  
  outFile <- file.path(outFolder, "nOverThresh_dt.Rdata")
  save(nOverThresh_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  cat(paste0("... load nOverThresh_dt\n"))
  outFile <- file.path(outFolder, "nOverThresh_dt.Rdata")
  nOverThresh_dt <- get(load(outFile))
  cat(paste0("... done\n"))
}


# load("N_FCC_OVER_QTILE_V2_EITHER_ITERATE/nOverThresh_dt.Rdata")

plotTit <- paste0("PERMG2T (keepPermut=", keepPermut, ")")

nOverThresh_dt$obs_over_rd <- nOverThresh_dt$ratioOverThresh_obs/nOverThresh_dt$ratioOverThresh_rd

save_dt <- nOverThresh_dt

# all_fccThresh <- all_fccThresh[all_fccThresh >= overall_min_fcc] 

nOverThresh_dt <- nOverThresh_dt[nOverThresh_dt$fccThresh >= overall_min_fcc,]

nOverThresh_dt$fccThresh <- round(nOverThresh_dt$fccThresh, 2)

nOverThresh_dt$fccThresh <- format(nOverThresh_dt$fccThresh, nsmall=2)

box_p <- ggboxplot(data = nOverThresh_dt,
                   x = "fccThresh", y = "obs_over_rd",
                   xlab="FCC threshold", ylab="Ratio TADs FCC >= thresh. obs/permut",
                   title = paste0(plotTit))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  geom_hline(yintercept=1, linetype=2, color="red")+
  theme(
    text = element_text(family=fontFamily),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=14),
    panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
        panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
        axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
        plot.title = element_text(hjust=0.5, size = 16, face="bold")
  )

outFile <- file.path(outFolder, paste0("ratioObsPermut_overThresh_boxplot.", plotType))
ggsave(box_p, file=outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


nOverThresh_dt$cmp <- all_cmps[paste0(nOverThresh_dt$exprds)]


box_p_cond <- ggboxplot(data = nOverThresh_dt,
                        x = "fccThresh", y = "obs_over_rd", 
                        color = "cmp",fill = "cmp",
                        xlab="FCC threshold", 
                        ylab="Ratio TADs FCC >= thresh. obs/permut",
                        title = paste0(plotTit))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  geom_hline(yintercept=1, linetype=2, color="red")+
  labs(color="", fill="")+
  theme(
    text = element_text(family=fontFamily),
    panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
        panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=14),
    axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
        plot.title = element_text(hjust=0.5, size = 16, face="bold")
  ) +
  scale_color_manual(values = all_cols) + 
  scale_fill_manual(values = all_cols) 

outFile <- file.path(outFolder, paste0("ratioObsPermut_overThresh_boxplot_cmpType.", plotType))
ggsave(box_p_cond, file=outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


mOverThresh_dt <- melt(nOverThresh_dt[,c("hicds", "exprds", "fccThresh", "ratioOverThresh_obs", "ratioOverThresh_rd")], 
                       id=c("hicds", "exprds", "fccThresh"))

mOverThresh_dt$variable <- gsub("ratioOverThresh_", "", mOverThresh_dt$variable)
mOverThresh_dt$variable <- gsub("rd", "random", mOverThresh_dt$variable)
mOverThresh_dt$variable <- gsub("obs", "observed", mOverThresh_dt$variable)


nDS <- length(unique(file.path(mOverThresh_dt$hicds, mOverThresh_dt$exprds)))

myTit <- paste0("Ratio of TADs over FCC thresh. (obs and permut)")
subTit <- paste0("(# datasets = ",  nDS, "; PERMUTG2T - # permut=", keepPermut, ")")

my_xlab <- "FCC threshold"
my_ylab <- "Ratio TADs FCC >= thresh."

require(ggsci)
ggsci_pal <- "lancet"
ggsci_subpal <- ""

all_th <- setNames(rank(unique(as.numeric(as.character(mOverThresh_dt$fccThresh)))), unique(mOverThresh_dt$fccThresh))

mOverThresh_dt$thresh_rank <- all_th[paste0(mOverThresh_dt$fccThresh)]
stopifnot(!is.na(mOverThresh_dt$thresh_rank))

head(mOverThresh_dt)


cor(mOverThresh_dt$thresh_rank, mOverThresh_dt$value) 
cor(as.numeric(as.character(mOverThresh_dt$fccThresh)), mOverThresh_dt$value)
# stopifnot(cor(mOverThresh_dt$thresh_rank, mOverThresh_dt$value) == 
#             cor(as.numeric(as.character(mOverThresh_dt$fccThresh)), mOverThresh_dt$value))

mOverThresh_dt$fccThresh <- as.numeric(as.character(mOverThresh_dt$fccThresh))
format_pval <- function(pval){
  pval <- scales::pvalue(pval, accuracy= 0.0001, add_p = TRUE)
  gsub(pattern = "(=|<)", replacement = " \\1 ", x = pval)
}

x_breaks <- sort(unique(mOverThresh_dt$fccThresh))
x_breaks2 <- x_breaks[seq(1, length(x_breaks), by = 2)]

point_p <- ggplot(mOverThresh_dt,aes(x=fccThresh, y = value, col = variable)) +
  geom_point()+
  ggtitle(myTit, subtitle=subTit)+
  # stat_summary(fun.data=mean_cl_normal) + 
  geom_smooth(method='lm', formula= y~x, se=F)  +
  labs(x = my_xlab, y = my_ylab, color = "")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous(breaks = x_breaks2,
                     labels = format(x_breaks2, nsmall=2))+
  my_box_theme+
  eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")")))+
  theme(
    panel.grid.major.x =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.grid.minor.x =  element_line(colour = "grey", size = 0.5, linetype=1),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size=10),
    axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=0.5),
    axis.line =  element_line()
  )+
  stat_cor(
    # aes(label = ifelse(..p.value.. < 0.05, formatC(..p.value.., format = "e", digits = 2),
    #                    ..p.value..)), 
    aes(label = paste(..r.label.., format_pval(..p..), sep = "~`,`~")),
    # aes(label = ifelse(p < 0.05,sprintf("p = %2.1e", as.numeric(..p.value..)), ..p.value..)),
    # aes(label=sprintf("p = %5.4f", as.numeric(..p.value..))),
    # aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                   
method = "pearson", label.x.npc = "left", label.y.npc = "bottom", p.digits=2)

outFile <- file.path(outFolder, paste0("ratioObsAndPermut_overThresh_pointLine.", plotType))
ggsave(point_p, file=outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


box_p <- ggplot(mOverThresh_dt,aes( x= factor(thresh_rank), y = value, col = variable)) +
  geom_boxplot()+
  ggtitle(myTit, subtitle=subTit)+
  labs(x = my_xlab, y = my_ylab, color = "")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_discrete(breaks = all_th[seq(1, length(all_th), by = 2)], labels=names(all_th)[seq(1, length(all_th), by = 2)])+
  my_box_theme+
  eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")")))+
  theme(
    # axis.text.x = element_text(angle=90),
    panel.grid.major.x =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.grid.minor.x =  element_line(colour = "grey", size = 0.5, linetype=1),
    axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=0.5),
    legend.key = element_rect(fill = NA),
    legend.text = element_text(size=10),
    axis.line =  element_line()
  )
outFile <- file.path(outFolder, paste0("ratioObsAndPermut_overThresh_boxplot.", plotType))
ggsave(box_p, file=outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

#***********************************************************************************
nOverThresh_dt$nOverThresh_rd_mean <- nOverThresh_dt$nOverThresh_rd/keepPermut

cOverThresh_dt <- melt(nOverThresh_dt[,c("hicds", "exprds", "fccThresh", "nOverThresh_obs", "nOverThresh_rd_mean")],
                       id=c("hicds", "exprds", "fccThresh"))
cOverThresh_dt$thresh_rank <- all_th[paste0(cOverThresh_dt$fccThresh)]
stopifnot(!is.na(cOverThresh_dt$thresh_rank))

cOverThresh_dt$variable <- gsub("nOverThresh_", "", cOverThresh_dt$variable)
cOverThresh_dt$variable <- gsub("rd_mean", "random (mean)", cOverThresh_dt$variable)
cOverThresh_dt$variable <- gsub("obs", "observed", cOverThresh_dt$variable)

myTit <- paste0("# of TADs over FCC thresh. (obs and permut)")
subTit <- paste0("(# datasets = ",  nDS, ")")

my_ylab <- "# TADs FCC >= thresh."

box_c <- ggplot(cOverThresh_dt,aes( x= factor(thresh_rank), y = value, col = variable)) +
  geom_boxplot()+
  ggtitle(myTit, subtitle=subTit)+
  labs(x = my_xlab, y = my_ylab, color = "")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_discrete(breaks = all_th[seq(1, length(all_th), by = 2)], labels=names(all_th)[seq(1, length(all_th), by = 2)])+
  my_box_theme+
  eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")")))+
  theme(
    # axis.text.x = element_text(angle=90),
    panel.grid.major.x =  element_line(colour = "grey", size = 0.5, linetype=1),
    panel.grid.minor.x =  element_line(colour = "grey", size = 0.5, linetype=1),
    axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=0.5),
    legend.text = element_text(size=10),
    legend.key = element_rect(fill = NA),
    axis.line =  element_line()
  )

outFile <- file.path(outFolder, paste0("nbrObsAndPermutMean_overThresh_boxplot.", plotType))
ggsave(box_c, file=outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


stop("-ok\n")


