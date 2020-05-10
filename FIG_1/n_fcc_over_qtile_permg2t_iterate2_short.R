# Rscript n_fcc_over_qtile_permg2t_iterate2_short.R

outFolder <- "N_FCC_OVER_QTILE_PERMG2T_ITERATE2_SHORT"
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


buildData <- FALSE

all_fccThresh <- seq(from = -1, to = 1, by = 0.05)


keepPermut <- 10


rd_type <- "PERMG2T"

buildData1 <- TRUE
buildData2 <- TRUE

all_hicds <- all_obs_hicds


if(buildData1){
  hicds = all_hicds[1]
  hicds = all_hicds[2]
  all_fcc_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    hicds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      cat("... start ", hicds, " - ", exprds, " \n")
      
      rd_file <- file.path(pipFolder,  hicds, exprds, "8cOnlyFCC_runAllDown/prodSignedRatio_permDT.Rdata")
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
} else {
  outFile <- file.path(outFolder, "all_fcc_dt.Rdata")
  all_fcc_dt <- get(load(outFile))
}

all_fcc_dt$dataset <- file.path(all_fcc_dt$hicds, all_fcc_dt$exprds)

overall_min_fcc <- min(all_fcc_dt$fcc_value) # -0.61


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
  
  outFile <- file.path(outFolder, "nOverThresh_dt.Rdata")
  nOverThresh_dt <- get(load(outFile))
  
}


# load("N_FCC_OVER_QTILE_V2_EITHER_ITERATE/nOverThresh_dt.Rdata")

plotTit <- paste0("PERMG2T (keepPermut=", keepPermut, ")")


nOverThresh_dt$obs_over_rd <- nOverThresh_dt$ratioOverThresh_obs/nOverThresh_dt$ratioOverThresh_rd

save_dt <- nOverThresh_dt

nOverThresh_dt <- nOverThresh_dt[nOverThresh_dt$fccThresh >= -0.5,]

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
                        xlab="FCC threshold", ylab="Ratio TADs FCC >= thresh. obs/permut",
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



ggplot(data,aes(x.plot, y.plot)) +
  stat_summary(fun.data=mean_cl_normal) + 
  geom_smooth(method='lm', formula= y~x)





