# setDir <- "/media/electron"
setDir <- ""

require(ggpubr)
require(ggsci)
require(foreach)
require(doMC)
registerDoMC(40)
require(patchwork)

# Rscript interactome_vPixel_all_v2.R

buildTable <- F

outFolder <- "INTERACTOME_VPIXEL_ALL_V2"
dir.create(outFolder, recursive = TRUE)

signif_level <- 0.05

# if a value, for ex, region_info[[1]]$diffInt_pmat$LNCaP_minus_RWPE1[2,3] equals 0.01, then it means: at the region 1, 
# the [2,3] pixel, the contact intensity (measured by hic-dc p-value, i.e., values in region_info[[1]]$hic_dc_mat) of
# LNCaP is higher than RWPE1, with a p-value of 0.01; if region_info[[1]]$diffInt_pmat$LNCaP_minus_RWPE1[2,3] equals -0.01,
# then it means, the contact intensity of LNCaP is lower than RWPE1, with a p-value of 0.01

# 22Rv1minusRWPE1 -> higher intensity in 22Rv1;

ref_hicds <- "RWPE1"
matching_hicds <- "22Rv1"

if(ref_hicds == "22Rv1" | ref_hicds =="LNCaP") {
  ds_dir2 <- matching_hicds
  ds_dir1 <- ref_hicds
  # matchingDir <- paste0(matching_hicds, "_minus_", ref_hicds)
} else {
  ds_dir1 <- matching_hicds
  ds_dir2 <- ref_hicds
  # matchingDir <- paste0(ref_hicds, "_minus_", matching_hicds)  
  # legText <- c(paste0(">0 = higher in ", matching_hicds),
  #              paste0("<0 = higher in ", ref_hicds))
}
matchingDir <- paste0(ds_dir1, "_minus_", ds_dir2)  
legText <- c(paste0(">0 = higher in ", ds_dir1),
             paste0("<0 = higher in ", ds_dir2))
left_text <- paste0("higher in ", ds_dir2)
right_text <- paste0("higher in ", ds_dir1)

all_files <- list.files(file.path(setDir,
                                  paste0("/mnt/ndata/Yuanlong/1.Projects/19.With_Marie/1.Data")), pattern="_40kb_TADs_for_hicdc_withSignif.output.v3.Rdata$",
                        full.names=TRUE)

my_file <- all_files[grepl(ref_hicds, basename(all_files))]
stopifnot(length(my_file) == 1)
dt <- get(load(my_file))
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
    if(is.null(pval_mat)) return(NULL) # for LNCaP -> some are null
    stopifnot(isSymmetric(as.matrix(pval_mat)))
    pval_mat <- as.numeric(as.matrix(pval_mat)[lower.tri(as.matrix(pval_mat), diag = TRUE)])
    abs_pval_mat_log10 <- -log10(abs(pval_mat))
    # retrieve back the sign
    pval_mat_log10 <- sign(pval_mat) * abs_pval_mat_log10
    stopifnot(sum(pval_mat_log10 > 0) == sum(pval_mat_log10 > 0))
    
    data.frame(
      tad_name = names(dt)[i],
      tad_dir = all_dirs[i],
      p_values = pval_mat,
      p_values_log10 = pval_mat_log10,
      stringsAsFactors = FALSE
    )
    
  }
  outFile <- file.path(outFolder, paste0("ref", ref_hicds, "_", matchingDir, "_plot_dt.Rdata"))
  save(plot_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, paste0("ref", ref_hicds, "_", matchingDir, "_plot_dt.Rdata"))
  plot_dt <- get(load(outFile))
}
################################################################################################## 
################################################# prep data

plot_dt$tad_dir_short <- plot_dt$tad_dir
plot_dt$tad_dir_short[grepl("notsignif",plot_dt$tad_dir_short )] <- "not signif." 

all_dirs_short <- all_dirs
all_dirs_short[grepl("notsignif",all_dirs_short )] <- "not signif." 
all_cmps <- setNames(c("","tumor; ", "normal; "),c("not signif.", "signif.up", "signif.down") )
n_vals <- setNames(as.numeric(table(all_dirs_short)), names(table(all_dirs_short)))
plot_dt$tad_dir_short <- paste0(plot_dt$tad_dir_short, "\n(", all_cmps[as.character(plot_dt$tad_dir_short)], 
                                "n=", n_vals[as.character(plot_dt$tad_dir_short)], ")")
stopifnot(!is.na(plot_dt$tad_dir_short))

plot_dt$tad_dir_short2 <- plot_dt$tad_dir
plot_dt$tad_dir_short2[grepl("notsignif",plot_dt$tad_dir_short2 )] <- "not signif." 
plot_dt$tad_dir_short2[grepl("^signif\\.",plot_dt$tad_dir_short2 )] <- "signif." 

all_dirs_short2 <- all_dirs
all_dirs_short2[grepl("notsignif",all_dirs_short2 )] <- "not signif." 
all_dirs_short2[grepl("^signif\\.",all_dirs_short2 )] <- "signif." 
n_vals2 <- setNames(as.numeric(table(all_dirs_short2)), names(table(all_dirs_short2)))
plot_dt$tad_dir_short2 <- paste0(plot_dt$tad_dir_short2, 
                                "\n(n=", n_vals2[as.character(plot_dt$tad_dir_short2)], ")")
stopifnot(!is.na(plot_dt$tad_dir_short2))

plotTit <- paste0(ref_hicds, " CDs")
subTit <- paste0(legText, collapse=";")

save(plot_dt, file="plot_dt.Rdata")

################################################################################################## 
################################################# density plots

plot_dt$tad_dir_short <- gsub("signif\\.up", "signif. up",plot_dt$tad_dir_short)
plot_dt$tad_dir_short <- gsub("signif\\.down", "signif. down",plot_dt$tad_dir_short)

mycols <- pal_jama()(5)[4:5]

plot_dt <- plot_dt[! grepl("notsignif", plot_dt$tad_dir),]

labtext_plotvar <- c(p_values="p-val.",
                     p_values_log10="p-val. [-log10]")

plot_var="p_values_log10"
for(plot_var in c("p_values", "p_values_log10")) {
  
  stopifnot(plot_var %in% names(labtext_plotvar))
  
  ##*********************************** full distribution
  
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
    scale_fill_manual(values=mycols)+
    scale_color_manual(values=mycols)+
    labs(x=paste0(matchingDir, " ", plot_var), y="Density", fill="", color="" )+
    ggtitle(plotTit, subtitle=subTit)+
    geom_vline(xintercept=0, linetype=2, color="darkgrey") +
    theme(  plot.subtitle=element_text(size=14, face="italic", hjust=0.5),
            plot.title=element_text(size=16, face="bold", hjust=0.5)) 
  
  if(plot_var == "p_values") {
    signif_lev <- signif_level
  }else if(plot_var == "p_values_log10") {
    signif_lev <- -log10(signif_level)
  }else{
    stop("error\n")
  }
  p_dens <- p_dens + 
    geom_vline(xintercept=c(-signif_lev,signif_lev), linetype=1, color="orange")
  
  # outFile <- file.path(outFolder, paste0("ref", ref_hicds, "_", plot_var, "_", matchingDir, "_density_by_tadDir.png"))
  # ggsave(p_dens, file=outFile, height=5.5, width = 6.5)
  # cat(paste0("... written: ", outFile, "\n"))
  
  

  
  
  
  ##*********************************** log 10
  
  density_data <- ggplot_build(p_dens)$data[[1]]
  
  
  if(plot_var == "p_values") {
    sub_data_ymax <- max(density_data[abs(density_data$x) <= signif_lev, "y"])
    p_cut <- p_dens +  coord_cartesian(xlim=c(-signif_lev,signif_lev), ylim=c(0, sub_data_ymax))
    
    widthGG <- 7.5
    
  }else if(plot_var == "p_values_log10") {
  
    sub_data_left_ymax <- max(density_data[density_data$x <= -signif_lev, "y"])
    sub_data_right_ymax <- max(density_data[density_data$x >= signif_lev, "y"])
    
    sub_data_left_xmin <- min(density_data[density_data$x <= -signif_lev, "x"])
    sub_data_right_xmax <- max(density_data[density_data$x >= signif_lev, "x"])

        # on electron, does not work with NA !!!
    p_left <- p_dens + 
      coord_cartesian(xlim=c(sub_data_left_xmin, -signif_lev), ylim=c(0, sub_data_left_ymax)) + 
      labs(x=paste0(left_text, " ",labtext_plotvar[plot_var] ))+
      # labs(x=paste0())+
      theme(plot.title = element_blank(),
            plot.subtitle = element_blank())
    p_right <- p_dens + 
      coord_cartesian(xlim=c(signif_lev, sub_data_right_xmax), ylim=c(0,sub_data_right_ymax))+
      # labs(x="left_x")+
      labs(x=paste0(right_text, " ",labtext_plotvar[plot_var] ))+
      theme(plot.title = element_blank(),
            plot.subtitle = element_blank())
    p_cut <- (p_left + p_right) +
      plot_layout(guides = "collect") & 
      plot_annotation(title=plotTit, subtitle=subTit, theme=theme(  plot.subtitle=element_text(size=14, face="italic", hjust=0.5),
                                                                    plot.title=element_text(size=16, face="bold", hjust=0.5))) &
      theme(legend.position = 'bottom', legend.text=element_text(size=10))
          
    
    widthGG <- 7.5
    
  }else{
    stop("error\n")
  }
  
  outFile <- file.path(outFolder, paste0("ref", ref_hicds, "_", plot_var, "_", matchingDir, "_density_by_tadDir_vCut.svg"))
  ggsave(p_cut, file=outFile, height=5.5, width = widthGG)
  cat(paste0("... written: ", outFile, "\n"))
  

    
}


