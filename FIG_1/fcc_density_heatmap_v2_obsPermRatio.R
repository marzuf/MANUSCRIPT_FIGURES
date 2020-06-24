# Rscript fcc_density_heatmap_v2_obsPermRatio.R

require(ggplot2)

plotType <- "png"

source("../settings.R")

outFolder <- file.path("FCC_DENSITY_HEATMAP_V2_OBSPERMRATIO")
dir.create(outFolder, recursive = TRUE)

auc_ratio_file <- "FCC_WAVE_PLOT_NOABS/all_fcc_dt.Rdata"
tmp <- get(load(auc_ratio_file))
tmp <- tmp[order(tmp$fcc_auc, decreasing = TRUE),]
ds_levels <- file.path(tmp$hicds, tmp$exprds)

dsCols <- all_cols[all_cmps[basename(ds_levels)]]

obs_dt <- get(load("FCC_DENSITY_HEATMAP_OBS_V2/density_plot_OBSERVED.Rdata"))
colnames(obs_dt)[colnames(obs_dt) == "density_y"] <- "density_y_OBS"

perm_dt <- get(load("FCC_DENSITY_HEATMAP_PERMDT_V2/density_plot_PERMG2T.Rdata"))
colnames(perm_dt)[colnames(perm_dt) == "density_y"] <- "density_y_PERMUT"

cmp_dt <- merge(obs_dt, perm_dt, by=c("dataset", "density_x"), all=TRUE)
stopifnot(!is.na(cmp_dt))

cmp_dt$density_y_Ratio <- log2(cmp_dt$density_y_OBS/cmp_dt$density_y_PERMUT)

nDS <- length(unique(cmp_dt$dataset))


myWidthGG=myWidthGG*1.1

subTit <- paste0("all datasets (n=", nDS, ")")

cmp_dt$dataset <- factor(cmp_dt$dataset, levels = ds_levels)

density_plot_s <- ggplot(cmp_dt, aes(x = dataset, y = density_x, fill = density_y_Ratio))+
  geom_tile() +
  ggtitle(paste0("FCC score distribution - OBS/PERMUT [log2]"),
          subtitle = subTit) +
  scale_x_discrete(name="Datasets ranked by decreasing AUC FCC ratio", labels = rep(labsymbol, nDS ), expand = c(0, 0))  +
  scale_y_continuous(name="FCC score",
                     breaks = scales::pretty_breaks(n = 20),  expand = c(0, 0))+
  labs(fill = "Density\nratio [log2]")+
                          theme(
	text = element_text(family=fontFamily),
                            axis.text.x = element_text(colour = dsCols, size=8),
                            axis.text.y= element_text(colour = "black", size=12),
                            axis.title.x = element_text(colour = "black", size=14, face="bold"),
                            axis.title.y = element_text(colour = "black", size=14, face="bold"),
                            plot.title = element_text(hjust=0.5, size=16, face="bold"),
                            plot.subtitle = element_text(hjust=0.5, size=14, face="italic"),
                            panel.background = element_rect(fill = "transparent")
                            # legend.background =  element_rect()
                          )  
density_plot <- density_plot_s + geom_vline(xintercept=seq(from=1.5, by=1, length.out = nDS-1), linetype=3)
                        
density_plot1 <- density_plot  +
  scale_fill_gradient( high="red", low="blue", na.value = "white" ) 
# +
#  scale_fill_gradient2( high="red", low="blue", na.value = "grey", mid ="white", midpoint=mean(cmp_dt$density_y_Ratio)) 


cat(myHeightGG, "\n")
cat(myWidthGG, "\n")

  density_plot2 <- density_plot + 
					scale_fill_viridis_c(option="A")  
  
  density_plot3 <- density_plot + 
    geom_vline(xintercept=seq(from=1.5, by=1, length.out = nDS-1), linetype=1)+
    scale_fill_gradientn(colours=colorRamps::matlab.like2(20))
    

outFile <- file.path(outFolder, paste0("FCC_score_ratio_dist_allDS_obs_perm_densityheatmap.", plotType))
ggsave(density_plot1, filename = outFile,  height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("FCC_score_ratio_dist_allDS_obs_perm_densityheatmap_vPal.", plotType))
ggsave(density_plot2, filename = outFile,  height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("FCC_score_ratio_dist_allDS_obs_perm_densityheatmap_vMatlab.", plotType))
ggsave(density_plot3, filename = outFile,  height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


all_dt <- get(load("FCC_DENSITY_HEATMAP_OBS_V2/all_result_dt.Rdata"))
min_fcc <- min(all_dt$FCC)
stopifnot(!is.na(min_fcc))


cmp_dt_cut <- cmp_dt[cmp_dt$density_x >= min_fcc,]

density_plot_s <- ggplot(cmp_dt_cut, aes(x = dataset, y = density_x, fill = density_y_Ratio))+
  geom_tile() +
  ggtitle(paste0("FCC score distribution - OBS/PERMUT [log2]"),
          subtitle = subTit) +
  scale_x_discrete(name="Datasets ranked by decreasing AUC FCC ratio", labels = rep(labsymbol, nDS ), expand = c(0, 0))  +
  scale_y_continuous(name="FCC score",
                     breaks = scales::pretty_breaks(n = 20),  expand = c(0, 0))+
  labs(fill = "Density\nratio [log2]")+
                          theme(
	text = element_text(family=fontFamily),
                            axis.text.x = element_text(colour = dsCols, size=8),
                            axis.text.y= element_text(colour = "black", size=12),
                            axis.title.x = element_text(colour = "black", size=14, face="bold"),
                            axis.title.y = element_text(colour = "black", size=14, face="bold"),
                            plot.title = element_text(hjust=0.5, size=16, face="bold"),
                            plot.subtitle = element_text(hjust=0.5, size=14, face="italic"),
                            panel.background = element_rect(fill = "transparent")
                            # legend.background =  element_rect()
                          )  
density_plot <- density_plot_s + geom_vline(xintercept=seq(from=1.5, by=1, length.out = nDS-1), linetype=3)
                        
density_plot1 <- density_plot  +
  scale_fill_gradient( high="red", low="blue", na.value = "white" ) 
# +
#  scale_fill_gradient2( high="red", low="blue", na.value = "grey", mid ="white", midpoint=mean(cmp_dt$density_y_Ratio)) 


  density_plot2 <- density_plot + 
					scale_fill_viridis_c(option="A")  

  density_plot3 <- density_plot_s + 
    geom_vline(xintercept=seq(from=1.5, by=1, length.out = nDS-1), linetype=1)+
    scale_fill_gradientn(colours=colorRamps::matlab.like2(20))
    
  
outFile <- file.path(outFolder, paste0("FCC_score_ratio_dist_allDS_obs_perm_densityheatmap_cutMin.", plotType))
ggsave(density_plot1, filename = outFile,  height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("FCC_score_ratio_dist_allDS_obs_perm_densityheatmap_vPal_cutMin.", plotType))
ggsave(density_plot2, filename = outFile,  height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, paste0("FCC_score_ratio_dist_allDS_obs_perm_densityheatmap_vMatlab_cutMin.", plotType))
ggsave(density_plot3, filename = outFile,  height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


#> load("FCC_DENSITY_HEATMAP_OBS_V2/all_result_dt.Rdata")
#> head(all_result_dt
#+ )
#                 hicds           exprds       FCC
#1 Barutcu_MCF-10A_40kb TCGAbrca_lum_bas 0.3017254
#2 Barutcu_MCF-10A_40kb TCGAbrca_lum_bas 0.4187001
#3 Barutcu_MCF-10A_40kb TCGAbrca_lum_bas 1.0000000
#4 Barutcu_MCF-10A_40kb TCGAbrca_lum_bas 1.0000000
#5 Barutcu_MCF-10A_40kb TCGAbrca_lum_bas 0.4083328
#6 Barutcu_MCF-10A_40kb TCGAbrca_lum_bas 0.2549637
#> min(all_result_dt$FCC)
#[1] -0.4820128
#> load("FCC_DENSITY_HEATMAP_PERMDT_V2/all_result_dt.Rdata")                                                                                                                                              
#> min(all_result_dt$FCC)
#[1] -0.6138406



