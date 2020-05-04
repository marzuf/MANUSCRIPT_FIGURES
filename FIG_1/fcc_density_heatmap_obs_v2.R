# Rscript fcc_density_heatmap_obs_v2.R 

# each column => dataset; color-coded by density

options(scipen = 100)

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

script_name <- "fcc_density_heatmap_obs_v2.R"
cat("> START ", script_name, "\n")
startTime <- Sys.time()

script0_name <- "0_prepGeneData"

require(foreach)
require(doMC)
require(ggpubr)
require(ggplot2)
registerDoMC(40)

plotType <- "png"
myHeight <- 400
myWidth <- 400
myHeightGG <- 7
myWidthGG <- 14

plotCex <- 1.4

auc_ratio_file <- "FCC_WAVE_PLOT_NOABS/all_fcc_dt.Rdata"
tmp <- get(load(auc_ratio_file))
tmp <- tmp[order(tmp$fcc_auc, decreasing = TRUE),]
ds_levels <- file.path(tmp$hicds, tmp$exprds)


source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")

pipFolder<- file.path("../../v2_Yuanlong_Cancer_HiC_data_TAD_DA/")
stopifnot(dir.exists(pipFolder))

source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

outFolder <- "FCC_DENSITY_HEATMAP_OBS_V2" 
dir.create(outFolder, recursive = TRUE)

all_hicds <- all_obs_hicds
all_exprds <- all_obs_exprds

resolution <- 0.01

buildData <- TRUE

if(buildData){
  hicds = all_hicds[1]
  hicds = all_hicds[2]
  all_result_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    hicds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      cat("... start ", hicds, " - ", exprds, " \n")
      
      fcc_file <- file.path(pipOutFolder, hicds, exprds, "8cOnlyFCC_runAllDown", "all_obs_prodSignedRatio.Rdata")
      if(!file.exists(fcc_file)) {
        fcc_file <- file.path(pipOutFolder, hicds, exprds, "8cOnlyFCConlyObs_runAllDown", "all_obs_prodSignedRatio.Rdata")  
      }
      if(!file.exists(fcc_file)) return(NULL)
      stopifnot(file.exists(fcc_file))
      all_fcc <- as.numeric(get(load(fcc_file)))
      # if(!file.exists(fcc_file)) {
      #   fcc_file <- file.path(pipOutFolder, hicds, exprds, "8cOnlyFCConlyObs_runAllDown", "all_obs_prodSignedRatio.Rdata")  
      # }
      
      data.frame(
        hicds = hicds,
        exprds=exprds,
        FCC = all_fcc,
        stringsAsFactors = FALSE
      )
    }
    hicds_dt
  }
  outFile <- file.path(outFolder, "all_result_dt.Rdata")
  save(all_result_dt, file = outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder, "all_result_dt.Rdata")
  all_result_dt <- get(load(outFile))
}  

min_fcc <- min(all_result_dt$FCC)

all_result_dt$dataset <- file.path(all_result_dt$hicds, all_result_dt$exprds)

all_result_dt$hicds_lab <- "OBSERVED"

all_types <- unique(all_result_dt$hicds_lab)
a_t = all_types[1]
for(a_t in all_types) {
  
  sub_dt <- all_result_dt[all_result_dt$hicds_lab == a_t,]

  # to have something for g2t permut
  # retrieve dataset order
  
  sub_dt$dataset <- factor(sub_dt$dataset, levels=ds_levels)
  stopifnot(!is.na(sub_dt$dataset))
  sub_dt$heatmap_x <- as.numeric(sub_dt$dataset) 
  dsCols <- all_cols[all_cmps[basename(ds_levels)]]
  

  sub_density_dt <- do.call(rbind, by(sub_dt, sub_dt$dataset, function(x) {
    hicds <- unique(x$hicds)
    stopifnot(length(hicds) == 1)
    exprds <- unique(x$exprds)
    stopifnot(length(exprds) == 1)
#    ds_dt <- density(x$FCC)
	ds_dt <- density(x$FCC, from=-1, to=1) # update 03.05.2020 -> limit density curves
    data.frame(
      hicds=hicds,
      exprds=exprds,
      density_x  = ds_dt$x,
      density_y  = ds_dt$y,
      stringsAsFactors = FALSE
    )
  }))
  sub_density_dt$dataset <- file.path(sub_density_dt$hicds, sub_density_dt$exprds)
  

  sub_density_dt$dataset <- factor(sub_density_dt$dataset, levels=ds_levels)
  stopifnot(!is.na(sub_density_dt$dataset))
  sub_density_dt$heatmap_x <- as.numeric(sub_density_dt$dataset)

  nDS <- length(unique(as.character(sub_dt$dataset)))

   new_density_x <- seq(min(sub_density_dt$density_x),max(sub_density_dt$density_x),by=resolution) # update 03.05.2020 -> limit density curves
#  new_density_x <- seq(-1,1,by=resolution)
  # new_density_x <- seq(min(all_result_dt$FCC),max(all_result_dt$FCC),by=resolution)

  plot_dt <- foreach(i = 1:nDS, .combine='rbind') %dopar% {
    sub_dt <- sub_density_dt[sub_density_dt$heatmap_x == i,]
    new_density <- approx(x=sub_dt$density_x,  y=sub_dt$density_y, xout = new_density_x, rule=2)
    data.frame(
      dataset=ds_levels[i],
      density_x = new_density[["x"]],
      density_y = new_density[["y"]],
      stringsAsFactors = FALSE
    )
  }
  plot_dt$dataset <- factor(plot_dt$dataset, levels=ds_levels)
  stopifnot(!is.na(plot_dt$dataset))

  density_plot <- ggplot(plot_dt, aes(x = dataset, y = density_x, fill = density_y))+
    geom_tile() +
  ggtitle(paste0("FCC score distribution - ", a_t),
          subtitle = paste0("all datasets (n=", nDS, ")")) +
    scale_x_discrete(name="Datasets ranked by decreasing AUC FCC ratio", labels = rep(labsymbol, nDS ), expand = c(0, 0))  +
    scale_y_continuous(name="FCC score",
                       breaks = scales::pretty_breaks(n = 20),  expand = c(0, 0))+
    labs(fill = "Density")+
    theme(
	text = element_text(family=fontFamily),
      axis.text.x = element_text(colour = dsCols, size=12),
      axis.text.y= element_text(colour = "black", size=12),
      axis.title.x = element_text(colour = "black", size=14, face="bold"),
      axis.title.y = element_text(colour = "black", size=14, face="bold"),
      plot.title = element_text(hjust=0.5, size=16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size=14, face="italic"),
      panel.background = element_rect(fill = "transparent")
      # legend.background =  element_rect()
    )+
	geom_vline(xintercept=seq(from=1.5, by=1, length.out = nDS-1), linetype=3) +

  density_plot1 <- density_plot + 
    # scale_fill_gradient( high="red", low="blue", na.value = "white")  +
    scale_fill_gradient2( high="red", low="blue", na.value = "grey", mid ="white", midpoint=mean(plot_dt$density_y))  

  density_plot2 <- density_plot + 
					scale_fill_viridis_c(option="A")  

  
  outFile <- file.path(outFolder, paste0("FCC_score_dist_allDS_", a_t, "_densityheatmap.", plotType))
  ggsave(density_plot1, filename = outFile,  height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))

  outFile <- file.path(outFolder, paste0("FCC_score_dist_allDS_", a_t, "_densityheatmap_vPal.", plotType))
  ggsave(density_plot2, filename = outFile,  height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))

  
  outFile <- file.path(outFolder, paste0("density_plot_", a_t, ".Rdata"))
  save(plot_dt, file = outFile)
  cat(paste0("... written: ", outFile, "\n"))



















  
  sub_dt <- all_result_dt[all_result_dt$hicds_lab == a_t,]
  sub_dt$FCC[sub_dt$hicds == dirname(ds_levels[1]) & sub_dt$exprds == basename(ds_levels[1])] <- sub_dt$FCC[sub_dt$hicds == dirname(ds_levels[1]) & sub_dt$exprds == basename(ds_levels[1])] - 0.5
    
  sub_density_dt <- do.call(rbind, by(sub_dt, sub_dt$dataset, function(x) {
    hicds <- unique(x$hicds)
    stopifnot(length(hicds) == 1)
    exprds <- unique(x$exprds)
    stopifnot(length(exprds) == 1)
#    ds_dt <- density(x$FCC)
    ds_dt <- density(x$FCC, from=-1, to=1) # TRIAL 4520
    data.frame(
      hicds=hicds,
      exprds=exprds,
      density_x  = ds_dt$x,
      density_y  = ds_dt$y,
      stringsAsFactors = FALSE
    )
  }))
  
  sub_density_dt$dataset <- file.path(sub_density_dt$hicds, sub_density_dt$exprds)
  
  
  sub_density_dt$dataset <- factor(sub_density_dt$dataset, levels=ds_levels)
  stopifnot(!is.na(sub_density_dt$dataset))
  sub_density_dt$heatmap_x <- as.numeric(sub_density_dt$dataset) 
  
  
  nDS <- length(unique(as.character(sub_dt$dataset)))
  
   new_density_x <- seq(min(sub_density_dt$density_x),max(sub_density_dt$density_x),by=resolution) # update 03.05.2020 -> limit density curves
#  new_density_x <- seq(-1,1,by=resolution)
  # new_density_x <- seq(min(all_result_dt$FCC),max(all_result_dt$FCC),by=resolution)
  
  
  plot_dt <- foreach(i = 1:nDS, .combine='rbind') %dopar% {
    sub_dt <- sub_density_dt[sub_density_dt$heatmap_x == i,]
    new_density <- approx(x=sub_dt$density_x,  y=sub_dt$density_y, xout = new_density_x, rule=2)
    data.frame(
      dataset=ds_levels[i],
      density_x = new_density[["x"]],
      density_y = new_density[["y"]],
      stringsAsFactors = FALSE
    )
  }
  plot_dt$dataset <- factor(plot_dt$dataset, levels=ds_levels)
  stopifnot(!is.na(plot_dt$dataset))
  
  density_plot <- ggplot(plot_dt, aes(x = dataset, y = density_x, fill = density_y))+ 
    geom_tile() +
    ggtitle(paste0("CHECKPLOT - FCC score distribution - ", a_t),   
            subtitle = paste0("all datasets (n=", nDS, ")")) +
    scale_x_discrete(name="Datasets ranked by decreasing AUC FCC ratio", labels = rep(labsymbol, nDS ), expand = c(0, 0))  +
    scale_y_continuous(name="FCC score",
                       breaks = scales::pretty_breaks(n = 20),  expand = c(0, 0))+
    labs(fill = "Density")+
    # scale_fill_gradient( high="red", low="blue", na.value = "white" )  +
    scale_fill_gradient2( high="red", low="blue", na.value = "grey", mid ="white", midpoint=mean(plot_dt$density_y))  +
    theme(
	text = element_text(family=fontFamily),
      axis.text.x = element_text(colour = dsCols, size=12),
      axis.text.y= element_text(colour = "black", size=12),
      axis.title.x = element_text(colour = "black", size=14, face="bold"),
      axis.title.y = element_text(colour = "black", size=14, face="bold"),
      plot.title = element_text(hjust=0.5, size=16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size=14, face="italic"),
      panel.background = element_rect(fill = "transparent")
      # legend.background =  element_rect()
    )
  
  density_plot <- density_plot + geom_vline(xintercept=seq(from=1.5, by=1, length.out = nDS-1), linetype=3)
  
  outFile <- file.path(outFolder, paste0("FCC_score_dist_allDS_", a_t, "_densityheatmap_checkplot.", plotType))
  ggsave(density_plot, filename = outFile,  height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  

  
  
  
}



















