setDir <- "/media/electron"
setDir <- ""

# Rscript interactome_vPixel.R

outFolder <- "INTERACTOME_VPIXEL"
dir.create(outFolder, recursive = TRUE)

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

i_tad=1
for(i_tad in seq_along(tads_to_plot)) {
  i <- which(
    all_chromos == tads_to_plot[[i_tad]][1] & 
      all_starts == tads_to_plot[[i_tad]][2] & 
      all_ends == tads_to_plot[[i_tad]][3] 
      )
  stopifnot(length(i) == 1)
  pval_mat <-   dt[[i]][["diffInt_pmat"]][[matchingDir]]
  stopifnot(isSymmetric(as.matrix(pval_mat)))
  pval_mat <- as.numeric(as.matrix(pval_mat)[lower.tri(as.matrix(pval_mat), diag = TRUE)])
  # save(pval_mat, file="pval_mat.Rdata", version=2)
  plottit <- paste0(ref_hicds, " -  ", names(dt)[i])
  outFile <- file.path(outFolder, paste0(ref_hicds, "_vs_", matching_hicds, 
                                         paste0("_chr", gsub(" +", "", gsub("\\.", "_", gsub(":", "_", gsub(": ", "_", names(dt)[i]))))), ".png"))
  do.call("png", list(file=outFile, height=400, width=500))
  # plot(density(as.numeric(pval_mat)), main=plottit)
  plot(density(pval_mat), main=plottit,
       xlab  = paste0(matchingDir),
       cex.axis=1.2, cex.lab=1.2)
  mtext(side=3, text=paste0(ref_hicds, "_vs_", matching_hicds, " (norm vs. tumor; up=tumor>norm; down=norm>tumor)"))
  abline(v=0, lty=2, col="darkgrey")
  abline(v=c(-0.05,0.05), lty=1, col="orange")
  legend("topleft", legend=legText, bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  abs_pval_mat_log10 <- -log10(abs(pval_mat))
  # retrieve back the sign
  pval_mat_log10 <- sign(pval_mat) * abs_pval_mat_log10
  
  stopifnot(sum(pval_mat_log10 > 0) == sum(pval_mat > 0))
  
  outFile <- file.path(outFolder, paste0(ref_hicds, "_vs_", matching_hicds, 
                                         paste0("_chr", gsub(" +", "", gsub("\\.", "_", gsub(":", "_", gsub(": ", "_", names(dt)[i]))))), "_log10.png"))
  do.call("png", list(file=outFile, height=400, width=500))
  plot(density(pval_mat_log10), main=plottit,
       xlab  = paste0(matchingDir, " [log10]"),
       cex.axis=1.2, cex.lab=1.2)
  mtext(side=3, text=paste0(ref_hicds, "_vs_", matching_hicds, " (norm vs. tumor; up=tumor>norm; down=norm>tumor)"))
  abline(v=0, lty=2, col="darkgrey")
  abline(v=c(-log10(0.05),log10(0.05)), lty=1, col="orange")
  legend("topleft", legend=legText, bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}

