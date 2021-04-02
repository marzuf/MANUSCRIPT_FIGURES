setDir <- "/media/electron"
# setDir <- ""

require(ggpubr)
require(ggsci)
require(foreach)
require(doMC)
registerDoMC(40)
require(patchwork)

# Rscript interactome_vPixel_all.R

buildTable <- F

outFolder <- "INTERACTOME_VPIXEL_ALL"
dir.create(outFolder, recursive = TRUE)

signif_level <- 0.05

# if a value, for ex, region_info[[1]]$diffInt_pmat$LNCaP_minus_RWPE1[2,3] equals 0.01, then it means: at the region 1, 
# the [2,3] pixel, the contact intensity (measured by hic-dc p-value, i.e., values in region_info[[1]]$hic_dc_mat) of
# LNCaP is higher than RWPE1, with a p-value of 0.01; if region_info[[1]]$diffInt_pmat$LNCaP_minus_RWPE1[2,3] equals -0.01,
# then it means, the contact intensity of LNCaP is lower than RWPE1, with a p-value of 0.01

# 22Rv1minusRWPE1 -> higher intensity in 22Rv1;

ref_hicds <- "RWPE1"
matching_hicds <- "22Rv1"

matchingDir <- paste0(matching_hicds, "_minus_", ref_hicds)

legText <- c(paste0(">0 = higher in ", matching_hicds),
             paste0("<0 = higher in ", ref_hicds))

dt <- get(load(file.path(setDir,
                         paste0("/mnt/ndata/Yuanlong/1.Projects/19.With_Marie/1.Data/GSE118514_", ref_hicds, "_40kb_TADs_for_hicdc_withSignif.output.v3.Rdata"))))
i=1
# chr12_TAD194	54160001	54440000
# chr7_TAD424	116080001	116320000
# chr17_TAD174	46720001	46880000
tads_to_plot <- list(
  c("chr12", 54160001, 54440000),
  c("chr7", 116080001, 116320000),
  c("chr17", 46720001, 46880000)
)

all_chromos <-  paste0("chr", gsub("^\\s+|\\s+$", "", unlist(lapply(strsplit(names(dt), split=":"), function(x) x[[1]]))))
stopifnot(!is.na(all_chromos))
all_starts <-  as.numeric(gsub("^\\s+|\\s+$", "", unlist(lapply(strsplit(names(dt), split=":"), function(x) x[[2]]))))
stopifnot(!is.na(all_starts))
all_ends <-  as.numeric(gsub("^\\s+|\\s+$", "", unlist(lapply(strsplit(names(dt), split=":"), function(x) x[[3]]))))
stopifnot(!is.na(all_ends))
all_dirs <- as.character(gsub("^\\s+|\\s+$", "", unlist(lapply(strsplit(names(dt), split=":"), function(x) x[[4]]))))
stopifnot(!is.na(all_dirs))
stopifnot(all_dirs %in% c("signif.up", "signif.down", "notsignif.up", "notsignif.down"))

i=1

if(buildTable) {
  
  plot_dt <- foreach(i = seq_along(dt), .combine='rbind') %dopar% {
    
    
    stopifnot(length(i) == 1)
    pval_mat <-   dt[[i]][["diffInt_pmat"]][[matchingDir]]
    stopifnot(isSymmetric(as.matrix(pval_mat)))
    pval_mat <- as.numeric(as.matrix(pval_mat)[lower.tri(as.matrix(pval_mat), diag = TRUE)])
    abs_pval_mat_log10 <- -log10(abs(pval_mat))
    # retrieve back the sign
    pval_mat_log10 <- sign(pval_mat) * abs_pval_mat_log10
    stopifnot(sum(pval_mat_log10 > 0) == sum(pval_mat_log10 > 0))
    
    data.frame(
      tad_dir = all_dirs[i],
      p_values = pval_mat,
      p_values_log10 = pval_mat_log10,
      stringsAsFactors = FALSE
    )
    
  }
  outFile <- file.path(outFolder, "plot_dt.Rdata")
  save(plot_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "plot_dt.Rdata")
  plot_dt <- get(load(outFile))
}

plot_dt$tad_dir_short <- plot_dt$tad_dir
plot_dt$tad_dir_short[grepl("notsignif",plot_dt$tad_dir_short )] <- "not signif." 

plotTit <- paste0(ref_hicds, " TADs")
subTit <- paste0(paste0(">0 = higher in ", matching_hicds),"; ",
             paste0("<0 = higher in ", ref_hicds))

plot_var="p_values_log10"
for(plot_var in c("p_values", "p_values_log10")) {
  
  p_dens <- ggdensity(plot_dt,
            x = paste0(plot_var),
            y = "..density..",
            # combine = TRUE,                  # Combine the 3 plots
            xlab = paste0(plot_var),
            # add = "median",                  # Add median line.
            rug = FALSE,                      # Add marginal rug
            color = "tad_dir_short",
            fill = "tad_dir_short",
            palette = "d3"
  ) +
    labs(x=paste0(matchingDir, " ", plot_var), y="Density", fill="", color="" )+
    ggtitle(plotTit, subtitle=subTit)+
    geom_vline(xintercept=0, linetype=2, color="darkgrey")
  
  if(plot_var == "p_values") {
    signif_lev <- signif_level
  }else if(plot_var == "p_values_log10") {
    signif_lev <- -log10(signif_level)
  }else{
    stop("error\n")
  }
  p_dens <- p_dens + 
    geom_vline(xintercept=c(-signif_lev,signif_lev), linetype=1, color="orange")
  
  outFile <- file.path(outFolder, paste0(plot_var, "_density_by_tadDir.png"))
  ggsave(p_dens, file=outFile, height=5.5, width = 6.5)
  cat(paste0("... written: ", outFile, "\n"))
  
  if(plot_var == "p_values") {
    p_cut <- p_dens +  scale_x_continuous(limits=c(-signif_lev,signif_lev))
    
    widthGG <- 6.5
    
    
  }else if(plot_var == "p_values_log10") {
  
    density_data <- ggplot_build(p_dens)$data[[1]]
    sub_data_left_ymax <- max(density_data[density_data$x <= -signif_lev, "y"])
    sub_data_right_ymax <- max(density_data[density_data$x >= signif_lev, "y"])
    
    p_left <- p_dens + coord_cartesian(xlim=c(NA, -signif_lev), ylim=c(NA, sub_data_left_ymax))
    p_right <- p_dens + coord_cartesian(xlim=c(signif_lev, NA), ylim=c(NA,sub_data_right_ymax))
    
    p_cut <- (p_left + p_right) +
      plot_layout(guides = "collect") & 
      theme(legend.position = 'bottom')
    
    widthGG <- 7.5
    
  }else{
    stop("error\n")
  }
  
  outFile <- file.path(outFolder, paste0(plot_var, "_density_by_tadDir_vCut.png"))
  ggsave(p_cut, file=outFile, height=5.5, width = widthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  
}




# plot_dt$facet <- as.numeric(abs(plot_dt[,plot_var]) >= signif_lev)
# plot_dt$facet <- plot_dt$facet * sign(plot_dt[,plot_var])
# stopifnot(sum(plot_dt$facet == 1) == sum(plot_dt[,plot_var] >= signif_lev))
# stopifnot(sum(plot_dt$facet == -1) == sum(plot_dt[,plot_var] <= -signif_lev))
# plot_dt2 <- plot_dt[plot_dt$facet != 0,]
# 
# ggdensity(plot_dt2,
#           x = paste0(plot_var),
#           y = "..density..",
#           # combine = TRUE,                  # Combine the 3 plots
#           xlab = paste0(plot_var),
#           # add = "median",                  # Add median line.
#           rug = FALSE,                      # Add marginal rug
#           color = "tad_dir_short",
#           fill = "tad_dir_short",
#           palette = "d3"
# ) +
#   facet_grid(. ~ facet, scales="free", space="free") +
#   labs(x=paste0(matchingDir, " ", plot_var), y="Density", fill="", color="" )+
#   ggtitle(plotTit, subtitle=subTit)+
  # geom_vline(xintercept=0, linetype=2, color="darkgrey")






# all_signifUp <- which(all_dirs == "signif.up")
# all_signifDown <- which(all_dirs == "signif.down")
# all_notSignif <- grepl("notsignif", all_dirs)
# 
# 
# lapply1
# 
# 
# i_tad=1
# for(i_tad in seq_along(tads_to_plot)) {
#   i <- which(
#     all_chromos == tads_to_plot[[i_tad]][1] & 
#       all_starts == tads_to_plot[[i_tad]][2] & 
#       all_ends == tads_to_plot[[i_tad]][3] 
#   )
#   stopifnot(length(i) == 1)
#   pval_mat <-   dt[[i]][["diffInt_pmat"]][[matchingDir]]
#   stopifnot(isSymmetric(as.matrix(pval_mat)))
#   pval_mat <- as.numeric(as.matrix(pval_mat)[lower.tri(as.matrix(pval_mat), diag = TRUE)])
#   # save(pval_mat, file="pval_mat.Rdata", version=2)
#   plottit <- paste0(ref_hicds, " -  ", names(dt)[i])
#   outFile <- file.path(outFolder, paste0(ref_hicds, "_vs_", matching_hicds, 
#                                          paste0("_chr", gsub(" +", "", gsub("\\.", "_", gsub(":", "_", gsub(": ", "_", names(dt)[i]))))), ".png"))
#   do.call("png", list(file=outFile, height=400, width=500))
#   # plot(density(as.numeric(pval_mat)), main=plottit)
#   plot(density(pval_mat), main=plottit,
#        xlab  = paste0(matchingDir),
#        cex.axis=1.2, cex.lab=1.2)
#   mtext(side=3, text=paste0(ref_hicds, "_vs_", matching_hicds, " (norm vs. tumor; up=tumor>norm; down=norm>tumor)"))
#   abline(v=0, lty=2, col="darkgrey")
#   abline(v=c(-0.05,0.05), lty=1, col="orange")
#   legend("topleft", legend=legText, bty="n")
#   foo <- dev.off()
#   cat(paste0("... written: ", outFile, "\n"))
#   
#   abs_pval_mat_log10 <- -log10(abs(pval_mat))
#   # retrieve back the sign
#   pval_mat_log10 <- sign(pval_mat) * abs_pval_mat_log10
#   
#   stopifnot(sum(pval_mat > 0) == sum(pval_mat_log10 > 0))
#   
#   outFile <- file.path(outFolder, paste0(ref_hicds, "_vs_", matching_hicds, 
#                                          paste0("_chr", gsub(" +", "", gsub("\\.", "_", gsub(":", "_", gsub(": ", "_", names(dt)[i]))))), "_log10.png"))
#   do.call("png", list(file=outFile, height=400, width=500))
#   plot(density(pval_mat_log10), main=plottit,
#        xlab  = paste0(matchingDir, " [log10]"),
#        cex.axis=1.2, cex.lab=1.2)
#   mtext(side=3, text=paste0(ref_hicds, "_vs_", matching_hicds, " (norm vs. tumor; up=tumor>norm; down=norm>tumor)"))
#   abline(v=0, lty=2, col="darkgrey")
#   abline(v=c(-log10(0.05),log10(0.05)), lty=1, col="orange")
#   legend("topleft", legend=legText, bty="n")
#   foo <- dev.off()
#   cat(paste0("... written: ", outFile, "\n"))
#   
# }
# 
