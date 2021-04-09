setDir <- "/media/electron"
setDir <- ""

# Rscript tad_matrix_plot.R

outFolder <- "TAD_MATRIX_PLOT"
dir.create(outFolder, recursive = TRUE)

binAround <- 2

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
tads_to_plot <- list(
  c("chr12", 54160001, 54440000)
)
all_chromos <-  paste0("chr", gsub("^\\s+|\\s+$", "", unlist(lapply(strsplit(names(dt), split=":"), function(x) x[[1]]))))
stopifnot(!is.na(all_chromos))
all_starts <-  as.numeric(gsub("^\\s+|\\s+$", "", unlist(lapply(strsplit(names(dt), split=":"), function(x) x[[2]]))))
stopifnot(!is.na(all_starts))
all_ends <-  as.numeric(gsub("^\\s+|\\s+$", "", unlist(lapply(strsplit(names(dt), split=":"), function(x) x[[3]]))))
stopifnot(!is.na(all_ends))

binSize <- 20000
mat_file <- file.path("extract_hic", "RWPE1_chr12_obs_KR_20kb_matrix.txt")
mat_dt <- read.delim(mat_file, stringsAsFactors = FALSE)


i_tad=1
for(i_tad in seq_along(tads_to_plot)) {
  
  tad_chromo <- tads_to_plot[[i_tad]][1] 
  tad_start <- as.numeric(tads_to_plot[[i_tad]][2] )
  tad_end <- as.numeric(tads_to_plot[[i_tad]][3] )
  
  stopifnot(!is.na(tad_start))
  stopifnot(!is.na(tad_end))
  
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
  
  plot_start <- (bin_start - binAround)
  plot_end <- (bin_end + binAround)
  
  plottit <- paste0(ref_hicds, " -  ", names(dt)[i])
  
  n <- (plot_end-plot_start+1)
  
  outFile <- file.path(outFolder, paste0( 
                                         paste0("RWPE1mat_chr", gsub(" +", "", gsub("\\.", "_", gsub(":", "_", gsub(": ", "_", names(dt)[i]))))), ".png"))
  do.call("png", list(file=outFile, height=400, width=400))
  image(1:n,1:n,
    as.matrix(log10(mat_dt[plot_start:plot_end,plot_start:plot_end])), 
    main=plottit, axes=F,
    xlab="", ylab="")
  mtext(side=3, text=paste0("binAround=", binAround))
  k=1
  t_hat <- c(binAround, binAround+bin_end-bin_start+1) + 1
  
  for (k in 1:(length(t_hat)-1)) {
    lines(c(t_hat[k],t_hat[k]),c(t_hat[k],(t_hat[(k+1)]-1)), col="darkgreen", lwd=5)
    lines(c(t_hat[(k+1)]-1,t_hat[(k+1)]-1),c(t_hat[k],t_hat[(k+1)]-1), col="darkgreen", lwd=5)
    lines(c(t_hat[k],t_hat[(k+1)]-1),c(t_hat[k],t_hat[k]), col="darkgreen", lwd=5)
    lines(c(t_hat[k],t_hat[(k+1)]-1),c(t_hat[(k+1)]-1,t_hat[(k+1)]-1), col="darkgreen", lwd=5)
  }
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  outFile <- file.path(outFolder, paste0(ref_hicds, "_vs_", matching_hicds, 
                                         paste0("RWPE1mat_chr", gsub(" +", "", gsub("\\.", "_", gsub(":", "_", gsub(": ", "_", names(dt)[i]))))), "_pvalmat.png"))
  do.call("png", list(file=outFile, height=400, width=400))
  
  image(as.matrix(pval_mat), 
        axes=F,
        xlab="", ylab="", main=plottit)
  mtext(side=3, text=paste0("p-val mat"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  outFile <- file.path(outFolder, paste0(ref_hicds, "_vs_", matching_hicds, 
                                         paste0("RWPE1mat_chr", gsub(" +", "", gsub("\\.", "_", gsub(":", "_", gsub(": ", "_", names(dt)[i]))))), "_pvalmat_log10.png"))
  do.call("png", list(file=outFile, height=400, width=400))
  
  image(as.matrix(log10(abs(pval_mat)) * sign(pval_mat)),
        axes=F,
        xlab="", ylab="", main=plottit)
  mtext(side=3, text=paste0("p-val mat [log10]"))
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  stopifnot(isSymmetric(as.matrix(pval_mat)))

  
}
