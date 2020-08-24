options(scipen=100)

SSHFS=F

# Rscript sameFam_signifTADs.R

script_name <- "sameFam_signifTADs.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(ggplot2)
require(ggpubr)
require(ggsci)
require(RColorBrewer)
require(reshape2)
# require(gplots)
registerDoMC(40)

jitterCol <- "blue"

hicds <- "ENCSR489OCU_NCI-H460_40kb"
exprds <- "TCGAluad_norm_luad"
family <- "hgnc"
familyType <- "hgnc_family_short"
# cannot take AUC_COEXPRDIST_WITHFAM_SORTNODUP/ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad_hgnc/hgnc_family_short/allData_dt.Rdata because crop dist


signifThresh <- 0.01


buildTable <- F

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 8)
myWidth <- myHeight
axisCex <- 1.4

myWidthGG <- 12
myHeightGG <- 12


outFolder <- "SAMEFAM_SIGNIFTADS"
dir.create(outFolder, recursive = TRUE)

family <- "hgnc"
familyType <- "hgnc_family_short"
all_plot_vars <- c("nbrFam", "ratioFam")
all_col_vars <- c("ratioUniqueFams", "nbrUniqueFams")

source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../settings.R")
source("../full_dataset_names.R")


final_dt <- get(load(file.path(runFolder, "CREATE_FINAL_TABLE/all_result_dt.Rdata")))



hicds="Rao_HCT-116_2017_40kb"
exprds="TCGAcoad_msi_mss"

if(buildTable) {


  all_sameFam_signif_dt <- foreach(hicds=all_obs_hicds, .combine='rbind') %dopar% {
    exprds_dt <- foreach(exprds=all_obs_exprds[[paste0(hicds)]], .combine='rbind') %dopar% {
      
      ds_dt <- final_dt[final_dt$hicds == hicds & final_dt$exprds == exprds,]
      stopifnot(nrow(ds_dt) > 0)
      signifTADs <- ds_dt$region[ds_dt$adjPvalComb <= signifThresh]
      
      
      gene_file <- file.path(pipFolder, hicds, exprds, "0_prepGeneData", "pipeline_geneList.Rdata")
      geneList <- get(load(gene_file))
      
      
      
      fam_dt <- get(load(file.path(runFolder, "PREP_GENE_FAMILIES_TAD_DATA", hicds, paste0(family, "_entrezID_family_TAD_DT.Rdata"))))
      
      pip_fam_dt <- fam_dt[fam_dt$entrezID %in% geneList,]
      stopifnot(pip_fam_dt$region %in% ds_dt$region)
      
      agg_nbrFam_dt <- aggregate(as.formula(paste0(familyType, " ~ region")), data=pip_fam_dt, FUN=function(x)length(unique(x)))
      colnames(agg_nbrFam_dt)[2] <- "nbrUniqueFams"
      
      agg_ratioFam_dt <- aggregate(as.formula(paste0(familyType, " ~ region")), data=pip_fam_dt, FUN=function(x)length(unique(x))/length(x))
      colnames(agg_ratioFam_dt)[2] <- "ratioUniqueFams"
      
      
      agg_nbrFam_dt$tad_signif <- ifelse(agg_nbrFam_dt$region %in% signifTADs, "signif.", "not signif.")
      agg_ratioFam_dt$tad_signif <- ifelse(agg_ratioFam_dt$region %in% signifTADs, "signif.", "not signif.")
      
      
      if(hicds=="ENCSR489OCU_NCI-H460_40kb") {
        
        for(plot_var in all_plot_vars){
          
          
          plot_dt <- eval(parse(text=paste0("agg_", plot_var, "_dt")))
          tad_labs <- c("signif.", "not signif.")
          my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], tad_labs)
          
          col_var <- colnames(plot_dt)[2]
          
          legTitle <- "TADs"
          
          plotTit <- paste0(col_var, " by TAD")
          
          mySub <- paste0(hicds, " - ",exprds, " - ", familyType, " data - TAD signif. thresh = ", signifThresh)
          
          p3 <- ggdensity(plot_dt,
                          x = col_var,
                          y = "..density..",
                          # combine = TRUE,                  # Combine the 3 plots
                          xlab = paste0(plot_var),
                          # add = "median",                  # Add median line.
                          rug = FALSE,                      # Add marginal rug
                          color = "tad_signif",
                          fill = "tad_signif",
                          palette = "jco"
          ) +
            ggtitle(plotTit, subtitle = mySub)+
            scale_color_manual(values=my_cols)+
            scale_fill_manual(values=my_cols)  +
            labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
            guides(color=FALSE)+
            scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
            scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
            theme(
              text = element_text(family=fontFamily),
              panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
              panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
              panel.background = element_rect(fill = "transparent"),
              panel.grid.major.x =  element_blank(),
              panel.grid.minor.x =  element_blank(),
              axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
              axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
              axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
              axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
              plot.title = element_text(hjust=0.5, size = 16, face="bold"),
              plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
              legend.title = element_text(face="bold")
            ) 
          
          outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", plot_var, "_signif_notSignif_density.", plotType))
          ggsave(p3, file=outFile, height=myHeight, width=myWidth)
          cat(paste0("... written: ", outFile, "\n"))

          outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", plot_var, "_signif_notSignif_boxplot.", plotType))
          do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth))
          par(mar = par()$mar + c(3,3,0,0))
          boxplot(as.formula(paste0(col_var, " ~ ", "tad_signif")), outline=FALSE,
                  data = plot_dt, main = plotTit,
                  xlab="", ylab="", cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, las=2)
          
          stripchart(as.formula(paste0(col_var, " ~ ", "tad_signif")), vertical = TRUE, data = plot_dt,
                     method = "jitter", add = TRUE, pch = 20, col = jitterCol)
          # legend("topright", legend=legText, bty="n", cex=0.9)
          
          mtext(side=2, text=paste0(col_var), cex=plotCex, line=5)
          mtext(side=3, text = paste0(hicds, " - ",exprds, " - signif. tresh = ", signifThresh), font = 3)
          foo <- dev.off()
          cat(paste0("... written: ", outFile, "\n"))
          
        } # end iterating plot_vars
      } # end if hicds  == NCI460
        
        out_dt <- merge(agg_ratioFam_dt, agg_nbrFam_dt,by="region", all=TRUE)
        out_dt$dataset <- file.path(hicds, exprds)
        stopifnot(!is.na(out_dt))
        stopifnot(out_dt$tad_signif.x == out_dt$tad_signif.y)        
        out_dt$tad_signif <- out_dt$tad_signif.x
        out_dt$tad_signif.x <- out_dt$tad_signif.y <- NULL
        # out_dt$tad_signif2 <- ifelse(out_dt$region %in% signifTADs, "signif.", "not signif.")
        # stopifnot(out_dt$tad_signif == out_dt$tad_signif2)        
        out_dt  
        
    }
    exprds_dt
  }
  outFile <- file.path(outFolder, "all_sameFam_signif_dt.Rdata")
  save(all_sameFam_signif_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  
} else {
  outFile <- file.path(outFolder, "all_sameFam_signif_dt.Rdata")
  all_sameFam_signif_dt <- get(load(outFile))
}

nDS <- length(unique(all_sameFam_signif_dt$dataset))
col_var=all_col_vars[1]
for(col_var in all_col_vars){
  
  
  tad_labs <- c("signif.", "not signif.")
  my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], tad_labs)
  

  legTitle <- "TADs"
  
  plotTit <- paste0(col_var, " by TAD")
  
  mySub <- paste0("# DS = ",nDS, " - ", familyType, " data - TAD signif. thresh = ", signifThresh)
  
  p3 <- ggdensity(all_sameFam_signif_dt,
                  x = col_var,
                  y = "..density..",
                  # combine = TRUE,                  # Combine the 3 plots
                  xlab = paste0(col_var),
                  # add = "median",                  # Add median line.
                  rug = FALSE,                      # Add marginal rug
                  color = "tad_signif",
                  fill = "tad_signif",
                  palette = "jco"
  ) +
    ggtitle(plotTit, subtitle = mySub)+
    scale_color_manual(values=my_cols)+
    scale_fill_manual(values=my_cols)  +
    labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
    guides(color=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme(
      text = element_text(family=fontFamily),
      panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major.x =  element_blank(),
      panel.grid.minor.x =  element_blank(),
      axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5),
      axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5),
      plot.title = element_text(hjust=0.5, size = 16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
      legend.title = element_text(face="bold")
    ) 
  
  outFile <- file.path(outFolder, paste0( "allDS_", col_var, "_signif_notSignif_density.", plotType))
  ggsave(p3, file=outFile, height=myHeight, width=myWidth)
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFolder, paste0("allDS_", col_var, "_signif_notSignif_boxplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth))
  par(mar = par()$mar + c(3,3,0,0))
  boxplot(as.formula(paste0(col_var, " ~ ", "tad_signif")), outline=FALSE,
          data = all_sameFam_signif_dt, main = plotTit,
          xlab="", ylab="", cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, las=2)
  
  stripchart(as.formula(paste0(col_var, " ~ ", "tad_signif")), vertical = TRUE, data = all_sameFam_signif_dt,
             method = "jitter", add = TRUE, pch = 20, col = jitterCol)
  # legend("topright", legend=legText, bty="n", cex=0.9)
  
  mtext(side=2, text=paste0(col_var), cex=plotCex, line=5)
  mtext(side=3, text = paste0("# DS = ",nDS,  " - signif. tresh = ", signifThresh), font = 3)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
}
