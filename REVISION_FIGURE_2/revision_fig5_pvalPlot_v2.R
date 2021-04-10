setDir <- "/media/electron"
setDir <- ""

options(scipen = 999)


# Rscript revision_fig5_pvalPlot_v2.R

outFolder <- "REVISION_FIG5_PVALPLOT_V2"
dir.create(outFolder, recursive = TRUE)


plotType <- "pdf"

if(plotType=="png") {
  myHeight <- 400
  myWidth <- 400
  
}else{
  myHeight <- 6
  myWidth <- 6
}


source("my_plot_matrix_v2.R")

plot_around <- c(0,0)


# if a value, for ex, region_info[[1]]$diffInt_pmat$LNCaP_minus_RWPE1[2,3] equals 0.01, then it means: at the region 1, 
# the [2,3] pixel, the contact intensity (measured by hic-dc p-value, i.e., values in region_info[[1]]$hic_dc_mat) of
# LNCaP is higher than RWPE1, with a p-value of 0.01; if region_info[[1]]$diffInt_pmat$LNCaP_minus_RWPE1[2,3] equals -0.01,
# then it means, the contact intensity of LNCaP is lower than RWPE1, with a p-value of 0.01

# 22Rv1minusRWPE1 -> higher intensity in 22Rv1;


binSize <- 20000

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
tads_to_plot <- list(
  c("chr12", 54160001, 54440000)
)
all_chromos <-  paste0("chr", gsub("^\\s+|\\s+$", "", unlist(lapply(strsplit(names(dt), split=":"), function(x) x[[1]]))))
stopifnot(!is.na(all_chromos))
all_starts <-  as.numeric(gsub("^\\s+|\\s+$", "", unlist(lapply(strsplit(names(dt), split=":"), function(x) x[[2]]))))
stopifnot(!is.na(all_starts))
all_ends <-  as.numeric(gsub("^\\s+|\\s+$", "", unlist(lapply(strsplit(names(dt), split=":"), function(x) x[[3]]))))
stopifnot(!is.na(all_ends))


i_tad=1
for(i_tad in seq_along(tads_to_plot)) {
  
  tad_chromo <- tads_to_plot[[i_tad]][1] 
  tad_start <- as.numeric(tads_to_plot[[i_tad]][2] )
  tad_end <- as.numeric(tads_to_plot[[i_tad]][3] )
  
  stopifnot(!is.na(tad_start))
  stopifnot(!is.na(tad_end))
  
  # Yuanlong practice: 0-based
  # 1 \t 40'000 should be converted to 1 \t 1
  # 40'001 \t 120'000 should be converted to 1 \t 1
  bin_start <- (tad_start - 1)/binSize  + 1
  bin_end <- tad_end/binSize
  
  
  i <- which(
    all_chromos == tad_chromo & 
      all_starts ==  as.numeric(tad_start),
    all_ends == tad_end 
  )
  stopifnot(length(i) == 1)
  pval_mat <-   dt[[i]][["diffInt_pmat"]][[matchingDir]]
  
  stopifnot(colnames(pval_mat) == bin_start:bin_end)  
  pval_mat <- as.matrix(pval_mat)
  
  tmp <- colnames(pval_mat)
  
  colnames(pval_mat) <- seq(from=tad_start-1, length.out=ncol(pval_mat), by=binSize)  # was 2709-2722, now is "54160000" "54420000" 
  
  stopifnot(colnames(pval_mat) == as.character((as.numeric(tmp)-1)*binSize))
  
  stopifnot(as.character(colnames(pval_mat)[ncol(pval_mat)]) == as.character(tad_end-binSize))
  rownames(pval_mat) <- seq(from=tad_start-1, length.out=ncol(pval_mat), by=binSize)
  
  ### to emulate hic data, replace colnames of yuanlong
  
  stopifnot(as.character(rownames(pval_mat)[ncol(pval_mat)]) == as.character(tad_end-binSize))
  
  pval_mat_log10 <- -log10(abs(pval_mat))
  signed_pval_mat_log10 <- pval_mat_log10*sign(pval_mat)
  
  stopifnot(colnames(signed_pval_mat_log10) == colnames(pval_mat))
  
  
  # for the plot function, colnames shoudl match
  # bin_start <- coord[1]/resolution
  # bin_end <- coord[2]/resolution -  1
  # lims <- seq(floor(bin_start), ceiling(bin_end)) * resolution
  
  
  outFile <- file.path(outFolder, paste0("ref", ref_hicds, "_", matchingDir, "_", tad_chromo, "_", tad_start, "_", tad_end, "_", binSize/1000, "kb_pvalplot.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))

  plotTit <- paste0(matchingDir)
  matLim <- ceiling(max(abs(pval_mat)))
  subTit <- paste0(ref_hicds, " - ", tad_chromo, ":", tad_start, "-", tad_end, "")
  
  my_plot_matrix(mat = pval_mat, 
                          # legCoords=  c(1,-1,1.5,10), #  xl,yb,xr,yt  c(1,1,1.5,2),
                         legCoords=  c(1,-0.5,1.5,2.5),
                         legMargins=c(0,0,0,0),
                 plotWithLeg=TRUE,
                 # tad_coord = c(tad_start,tad_end-binSize),
                 tad_coord = c(tad_start,tad_end),
                 minColRange=-matLim, maxColRange=matLim,
                 subTit=subTit, legTit =paste0("diff. interact.\np-values:"),
                 color = colorRampPalette(c("blue", "white", "red"))(100), 
                 resolution = binSize, 
                 transformation=identity,
                 bins_around=c(plot_around, plot_around),
                 main=plotTit) 
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  outFile <- file.path(outFolder, paste0("ref", ref_hicds, "_", matchingDir, "_", tad_chromo, "_", tad_start, "_", tad_end, "_", binSize/1000, "kb_pvalplot_log10.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  
  matLim <- ceiling(max(abs(signed_pval_mat_log10)))
  
  my_plot_matrix(mat = signed_pval_mat_log10, 
                 # legCoords=  c(1,-1,1.5,10), #  xl,yb,xr,yt  c(1,1,1.5,2),
                 legCoords=  c(1,-0.5,1.5,2.5),
                 legMargins=c(0,0,0,0),
                 plotWithLeg=TRUE,
                 # tad_coord = c(tad_start,tad_end-binSize),
                 tad_coord = c(tad_start,tad_end),
                 minColRange=-matLim, maxColRange=matLim,
                 subTit=subTit, legTit =paste0("diff. interact.\np-values\n[-log10]:"),
                 color = colorRampPalette(c("blue", "white", "red"))(100), 
                 resolution = binSize, 
                 transformation=identity,
                 bins_around=c(plot_around, plot_around),
                 main=plotTit) 
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  

    
}