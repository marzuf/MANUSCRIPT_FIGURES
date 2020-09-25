options(scipen=100)

SSHFS=F

# Rscript sameFam_signifTADs_v3.R

script_name <- "sameFam_signifTADs_v3.R"

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

ratioAnnotFilter <- 0.5



buildTable <- F

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 8)
myWidth <- myHeight
axisCex <- 1.4

myWidthGG <- 7
myHeightGG <- 7


outFolder <- "SAMEFAM_SIGNIFTADS_V3"
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

curr_theme <-     theme(
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

script0_name <- "0_prepGeneData"


custom_p <- function(p) {
  p2 <- p +  labs(fill="")+ scale_fill_nejm() + theme(plot.title=element_text(size=14,face="bold", hjust=0.5 ))
  return(p2)
}


if(buildTable) {


  all_sameFam_signif_dt <- foreach(hicds=all_obs_hicds, .combine='rbind') %dopar% {
    exprds_dt <- foreach(exprds=all_obs_exprds[[paste0(hicds)]], .combine='rbind') %dopar% {
      
      
      gene2tad_dt <- read.table(file.path(runFolder, hicds, "genes2tad", "all_genes_positions.txt"), col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = FALSE, header=F)
      gene2tad_dt$entrezID <- as.character(gene2tad_dt$entrezID)
      
      ds_dt <- final_dt[final_dt$hicds == hicds & final_dt$exprds == exprds,]
      stopifnot(nrow(ds_dt) > 0)
      signifTADs <- ds_dt$region[ds_dt$adjPvalComb <= signifThresh]
      
      tad_adjPvals <- setNames(ds_dt$adjPvalComb,ds_dt$region)
      
      gene_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
      geneList <- get(load(gene_file))
      
      stopifnot(geneList %in% gene2tad_dt$entrezID)
      gene2tad_dt <- gene2tad_dt[gene2tad_dt$entrezID %in% geneList,]
      
      nGenesByTAD <- setNames(as.numeric(table(gene2tad_dt$region)), names(table(gene2tad_dt$region)))
      
      fam_dt <- get(load(file.path(runFolder, "PREP_GENE_FAMILIES_TAD_DATA", hicds, paste0(family, "_entrezID_family_TAD_DT.Rdata"))))
      
      pip_fam_dt <- fam_dt[fam_dt$entrezID %in% geneList,]
      stopifnot(pip_fam_dt$region %in% ds_dt$region)
      
      agg_nbrAnnot_dt <- aggregate(as.formula(paste0(familyType, " ~ region")), data=pip_fam_dt, FUN=function(x)length(x))
      colnames(agg_nbrAnnot_dt)[2] <- "nbrAnnot"
      
      
      agg_nbrFam_dt <- aggregate(as.formula(paste0(familyType, " ~ region")), data=pip_fam_dt, FUN=function(x)length(unique(x)))
      colnames(agg_nbrFam_dt)[2] <- "nbrUniqueFams"
      
      agg_ratioFam_dt <- aggregate(as.formula(paste0(familyType, " ~ region")), data=pip_fam_dt, FUN=function(x)length(unique(x))/length(x))
      colnames(agg_ratioFam_dt)[2] <- "ratioUniqueFams"
      
      
      
      out_dt <- merge(agg_nbrAnnot_dt, merge(agg_ratioFam_dt, agg_nbrFam_dt,by="region", all=TRUE),by="region", all=TRUE)
      out_dt$dataset <- file.path(hicds, exprds)
      out_dt$tad_signif <- ifelse(out_dt$region %in% signifTADs, "signif.", "not signif.")
      stopifnot(!is.na(out_dt))
      out_dt$region <- as.character(out_dt$region)
      stopifnot(out_dt$region %in% names(tad_adjPvals))
      out_dt$tad_adjPval <- tad_adjPvals[out_dt$region]
      
      stopifnot(out_dt$region %in% names(nGenesByTAD))
      out_dt$nbrGenes <- nGenesByTAD[out_dt$region]
      
      out_dt$ratioAnnoted <- out_dt$nbrAnnot/out_dt$nbrGenes
      stopifnot(out_dt$ratioAnnoted >= 0 & out_dt$ratioAnnoted <= 1)
      
      
      kept_dt <- out_dt[out_dt$ratioAnnoted >= ratioAnnotFilter,]
      
      cat(nrow(out_dt), " -> ", nrow(kept_dt), "\n")
      
      kept_dt$nbrUniqueFams_lab <- ifelse(kept_dt$nbrUniqueFams == 1, "1 unique fam.", 
                                      ifelse(kept_dt$nbrUniqueFams == 2, "2 unique fam.",
                                             ifelse(kept_dt$nbrUniqueFams == 3, "3 unique fam.", 
                                                  ifelse(kept_dt$nbrUniqueFams > 3, ">= 4 unique fam.", NA))))
      
      stopifnot(!is.na(kept_dt$nbrUniqueFams_lab))
      
      agg_dt <- aggregate(region ~ dataset + tad_signif + nbrUniqueFams_lab, data=kept_dt, FUN = length)
      colnames(agg_dt)[colnames(agg_dt) == "region"] <- "nbr"
      agg_dt$nbr_log10 <- log10(agg_dt$nbr)
      
      nSignif_count <- setNames(as.numeric(table(kept_dt$tad_signif)), names(table(kept_dt$tad_signif)))
      agg_dt$ratio_nbr <- agg_dt$nbr/nSignif_count[agg_dt$tad_signif]
      
      
      subTit <- paste0(names(nSignif_count), " : ", as.numeric(nSignif_count), collapse = "; ")
      subTit <- paste0("TADs with ratioAnnot >=  ",ratioAnnotFilter, " - ",  subTit)
      
      
      agg_dt$nbrUniqueFams_lab <- factor(agg_dt$nbrUniqueFams_lab, levels=rev(c("1 unique fam.", "2 unique fam.","3 unique fam." , ">= 4 unique fam.")))
      stopifnot(!is.na(agg_dt$nbrUniqueFams_lab))
      
      if(hicds=="ENCSR489OCU_NCI-H460_40kb" & exprds=="TCGAluad_norm_luad") {
        
        
        p1 <- ggbarplot(data=agg_dt, y="nbr", x="tad_signif",
                        fill="nbrUniqueFams_lab", 
                        xlab="TADs", 
                        title = paste0(hicds, " - ", exprds),
                        subtitle=subTit) 
        p1 <- custom_p(p1)
        outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_nbr.svg"))
        ggsave(p1, filename=outFile, height=7, width=7)
        
        p2 <- ggbarplot(data=agg_dt, y="nbr_log10", x="tad_signif", fill="nbrUniqueFams_lab", xlab="TADs", 
                        title = paste0(hicds, " - ", exprds),
                        subtitle=subTit) 
        p2 <- custom_p(p2)
        outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_nbr_log.svg"))
        ggsave(p2, filename=outFile, height=7, width=7)
        
        p3 <- ggbarplot(data=agg_dt, y="ratio_nbr", x="tad_signif", fill="nbrUniqueFams_lab", xlab="TADs", 
                        title = paste0(hicds, " - ", exprds),
                        subtitle=subTit) 
        p3 <- custom_p(p3)
        outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_nbr_ratio.svg"))
        ggsave(p3, filename=outFile, height=7, width=7)
        
 
      } # end if hicds  == NCI460
        
      kept_dt  
        
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

agg_dt <- aggregate(region ~  tad_signif + nbrUniqueFams_lab, data=all_sameFam_signif_dt, FUN = length)
colnames(agg_dt)[colnames(agg_dt) == "region"] <- "nbr"
agg_dt$nbr_log10 <- log10(agg_dt$nbr)

nSignif_count <- setNames(as.numeric(table(all_sameFam_signif_dt$tad_signif)), names(table(all_sameFam_signif_dt$tad_signif)))
agg_dt$ratio_nbr <- agg_dt$nbr/nSignif_count[agg_dt$tad_signif]


agg_dt$nbrUniqueFams_lab <- factor(agg_dt$nbrUniqueFams_lab, levels=rev(c("1 unique fam.", "2 unique fam.", "3 unique fam.",">= 4 unique fam.")))
stopifnot(!is.na(agg_dt$nbrUniqueFams_lab))


subTit <- paste0(names(nSignif_count), " : ", as.numeric(nSignif_count), collapse = "; ")
subTit <- paste0("TADs with ratioAnnot >=  ",ratioAnnotFilter, " - ",  subTit)


  outFile <- file.path(outFolder, "agg_dt.Rdata")
  save(agg_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
  

nDS <- length(unique(all_sameFam_signif_dt$dataset))
  
  
  p1 <- ggbarplot(data=agg_dt, y="nbr", x="tad_signif",
                  fill="nbrUniqueFams_lab", 
                  xlab="TADs", 
                  title = paste0("all datasets - n=", nDS),
                  subtitle=subTit)
  p1 <- custom_p(p1)
  
  outFile <- file.path(outFolder, "all_nbr.svg")
  ggsave(p1, filename=outFile, height=7, width=7)
  
  p2 <- ggbarplot(data=agg_dt, y="nbr_log10", x="tad_signif", fill="nbrUniqueFams_lab", xlab="TADs", 
                  title =paste0("all datasets - n=", nDS),
                  subtitle=subTit) 
  p2 <- custom_p(p2)
  outFile <- file.path(outFolder, "all_nbr_log.svg")
  ggsave(p2, filename=outFile, height=7, width=7)
  
  
  p3 <- ggbarplot(data=agg_dt, y="ratio_nbr", x="tad_signif", fill="nbrUniqueFams_lab", xlab="TADs", 
                  title =paste0("all datasets - n=", nDS),
                  subtitle=subTit) 
  p3 <- custom_p(p3)
  outFile <- file.path(outFolder, "all_nbr_ratio.svg")
  ggsave(p3, filename=outFile, height=7, width=7)



# nDS <- length(unique(all_sameFam_signif_dt$dataset))
# col_var=all_col_vars[1]
# 
# atLeastAnnot <- c(0,3)
# 
# init_all_sameFam_signif_dt <- all_sameFam_signif_dt
# 
# for(min_annot in  atLeastAnnot) {
#   
#   all_sameFam_signif_dt <- init_all_sameFam_signif_dt
#   all_sameFam_signif_dt <-  all_sameFam_signif_dt[all_sameFam_signif_dt$nbrAnnot >= min_annot,]
#   
#   for(col_var in all_col_vars){
#     
#     
#     tad_labs <- c("signif.", "not signif.")
#     my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], tad_labs)
#     
#     
#     legTitle <- "TADs"
#     
#     plotTit <- paste0(col_var, " by TAD (# annot >= ", min_annot, ")")
#     
#     mySub <- paste0("# DS = ",nDS, " - ", familyType, " data - TAD signif. thresh = ", signifThresh)
#     
#     p3 <- ggdensity(all_sameFam_signif_dt,
#                     x = col_var,
#                     y = "..density..",
#                     # combine = TRUE,                  # Combine the 3 plots
#                     xlab = paste0(col_var),
#                     # add = "median",                  # Add median line.
#                     rug = FALSE,                      # Add marginal rug
#                     color = "tad_signif",
#                     fill = "tad_signif",
#                     palette = "jco"
#     ) +
#       ggtitle(plotTit, subtitle = mySub)+
#       scale_color_manual(values=my_cols)+
#       scale_fill_manual(values=my_cols)  +
#       labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
#       guides(color=FALSE)+
#       scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
#       scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
#       curr_theme
#     
#     outFile <- file.path(outFolder, paste0( "allDS_", col_var, "_signif_notSignif_minAnnot", min_annot, "_density.", plotType))
#     ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
#     cat(paste0("... written: ", outFile, "\n"))
#     
#     outFile <- file.path(outFolder, paste0("allDS_", col_var, "_signif_notSignif_minAnnot", min_annot, "_boxplot.", plotType))
#     do.call(plotType, list(outFile, height=myHeight*1.2, width=myWidth))
#     par(mar = par()$mar + c(3,3,0,0))
#     boxplot(as.formula(paste0(col_var, " ~ ", "tad_signif")), outline=FALSE,
#             data = all_sameFam_signif_dt, main = plotTit,
#             xlab="", ylab="", cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, las=2)
#     
#     stripchart(as.formula(paste0(col_var, " ~ ", "tad_signif")), vertical = TRUE, data = all_sameFam_signif_dt,
#                method = "jitter", add = TRUE, pch = 20, col = jitterCol)
#     # legend("topright", legend=legText, bty="n", cex=0.9)
#     
#     mtext(side=2, text=paste0(col_var), cex=plotCex, line=5)
#     mtext(side=3, text = paste0("# DS = ",nDS,  " - signif. tresh = ", signifThresh), font = 3)
#     foo <- dev.off()
#     cat(paste0("... written: ", outFile, "\n"))
#     
#     
#     myx <- all_sameFam_signif_dt[,paste0(col_var)]
#     myy <- -log10(all_sameFam_signif_dt[,"tad_adjPval"])
#     
#     
#     outFile <- file.path(outFolder, paste0("allDS_", col_var, "_vs_adjPval_minAnnot", min_annot, "_densplot.", plotType))
#     do.call(plotType, list(outFile, height=myHeight, width=myWidth))
#     # par(mar = par()$mar + c(3,3,0,0))
#     densplot(x=myx, y=myy,
#              main = plotTit,
#              xlab=paste0(col_var), ylab=paste0("TAD adj. p-val [-log10]"), cex.main=plotCex, cex.lab=plotCex, cex.axis=plotCex, las=2)
#     mtext(side=2, text=paste0(col_var), cex=plotCex, line=5)
#     mtext(side=3, text = paste0("# DS = ",nDS), font = 3)
#     addCorr(x=myx,y=myy,legPos="topright", bty="n")
#     foo <- dev.off()
#     cat(paste0("... written: ", outFile, "\n"))
#     
#     
#   }
#   
#   
#   
# }
