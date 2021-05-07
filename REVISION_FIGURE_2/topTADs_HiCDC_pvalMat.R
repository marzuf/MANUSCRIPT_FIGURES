

# Rscript topTADs_HiCDC_pvalMat.R

setDir <- "/media/electron"
setDir <- ""

binSize <- 20000
options(scipen=100)


plotType <- "png"
myHeight <- 400
myWidth <- 400

plot_around <- 0
nTop <- 10

source("my_plot_matrix_v2.R")


outFolder <-file.path("TOPTADS_HICDC_PVALMAT")
dir.create(outFolder,recursive = TRUE)

all_files <- list.files(file.path(setDir,
                                  paste0("/mnt/ndata/Yuanlong/1.Projects/19.With_Marie/1.Data")),
                        pattern="_40kb_TADs_for_hicdc_withSignif.output.v3.Rdata$",
                        full.names=TRUE)

runFolder <- "../../v2_Yuanlong_Cancer_HiC_data_TAD_DA"
final_table_file <- file.path(runFolder, "CREATE_FINAL_TABLE/all_result_dt.Rdata")

result_dt <- get(load(final_table_file))
result_dt$chromo <- gsub("(chr.+)_TAD.+", "\\1", result_dt$region)

ref_hicds = "22Rv1"
ref_hicds <- commandArgs(trailingOnly=TRUE)
match_hicds <- setdiff(c("22Rv1", "RWPE1"), ref_hicds)


all_files <- list.files(file.path(setDir,
                                  paste0("/mnt/ndata/Yuanlong/1.Projects/19.With_Marie/1.Data")), pattern="_40kb_TADs_for_hicdc_withSignif.output.v3.Rdata$",
                        full.names=TRUE)
ref_file <- all_files[grepl(ref_hicds, basename(all_files))]
stopifnot(length(ref_file) == 1)
ref_data <- get(load(ref_file))

all_chromos <-  paste0("chr", gsub("^\\s+|\\s+$", "", unlist(lapply(strsplit(names(ref_data), split=":"), function(x) x[[1]]))))
stopifnot(!is.na(all_chromos))
all_starts <-  as.numeric(gsub("^\\s+|\\s+$", "", unlist(lapply(strsplit(names(ref_data), split=":"), function(x) x[[2]]))))
stopifnot(!is.na(all_starts))
all_ends <-  as.numeric(gsub("^\\s+|\\s+$", "", unlist(lapply(strsplit(names(ref_data), split=":"), function(x) x[[3]]))))
stopifnot(!is.na(all_ends))


# match_file <- all_files[grepl(match_hicds, basename(all_files))]
# stopifnot(length(match_file) == 1)
# match_data <- get(load(match_file))


ds_dt <- result_dt[grepl(ref_hicds, result_dt$hicds),]
stopifnot(nrow(ds_dt) > 0)
stopifnot(length(unique(ds_dt$hicds)) == 1)
ds_dt <- ds_dt[order(ds_dt$adjPvalComb),]

i_top  = 1

for(i_top in 1:nTop) {
  
  
  
  
  # STEP 1 : retrieve the top 10
  tad_region <- ds_dt$region[i_top]
  tad_chromo <- ds_dt$chr[i_top]
  tad_start <- ds_dt$start[i_top]
  tad_end <- ds_dt$end[i_top]
  
  stopifnot(tad_end > tad_start)
  stopifnot(tad_end %% binSize == 0)
  stopifnot( (tad_start-1) %% binSize == 0)
  
  
  stopifnot(is.numeric(tad_start))
  stopifnot(is.numeric(tad_end))
  
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
  
  # matFile <- file.path(outFolder, paste0("ref", ref_hicds, "_", tad_chromo, "_", tad_start, "_", tad_end, "_", binSize/1000, 
  #                                        "kb_22Rv1_hicdc_mat.Rdata"))
  ref_hicdcmat <- ref_data[[i]][["hic_dc_mat"]][[paste0(ref_hicds)]]
  stopifnot(isSymmetric(as.matrix(ref_hicdcmat)))
  
  cat(paste0("tad_start = ", tad_start, "\n"))
  cat(paste0("tad_end = ", tad_end, "\n"))
  
  cat(paste0("(tad_end-tad_start+1)/binSize = ", (tad_end-tad_start+1)/binSize, "\n"))
  
  
  cat(paste0("bin_start = ", bin_start, "\n"))
  cat(paste0("bin_end = ", bin_end, "\n"))
  cat(paste0("range(colnames) = ", colnames(ref_hicdcmat)[1], " - " , colnames(ref_hicdcmat)[ncol(ref_hicdcmat)] , "\n"))
  cat(paste0("ncol(ref_hicdcmat) = ", ncol(ref_hicdcmat), "\n"))
  
  tmp <- as.numeric(as.character(colnames(ref_hicdcmat)))
  stopifnot(tmp == bin_start:bin_end)
  
  stopifnot(as.character(bin_start) == colnames(ref_hicdcmat)[1])
  stopifnot(as.character(bin_end) == colnames(ref_hicdcmat)[ncol(ref_hicdcmat)])
  stopifnot(ncol(ref_hicdcmat) == (tad_end-tad_start+1)/binSize)
  
  # save(ref_hicdcmat, file = matFile)
  # cat(paste0("... written:" , matFile, "\n"))
  
  # matFile <- file.path(outFolder, paste0("ref", ref_hicds, "_", tad_chromo, "_", tad_start, "_", tad_end, "_", binSize/1000, 
  #                                        "kb_RWPE1_hicdc_mat.Rdata"))
  match_hicdcmat <- ref_data[[i]][["hic_dc_mat"]][[paste0(match_hicds)]]
  stopifnot(isSymmetric(as.matrix(match_hicdcmat)))
  # save(hicdcmat_RWPE1, file = matFile)
  
  
  ref_htc_mat <- ref_hicdcmat
  nBins <- bin_end - bin_start + 1
  stopifnot(nBins == ncol(ref_hicdcmat))
  
  
  # colnames(ref_htc_mat) <- seq(from=(tad_start-1)0, by=binSize, length.out=ncol(ref_htc_mat))
  # rownames(ref_htc_mat) <- seq(from=0, by=binSize, length.out=nrow(ref_htc_mat))
  # # replace the colnames to have coord instead of bin
  # stopifnot(colnames(ref_htc_mat) == bin_start:bin_end)  
  # ref_htc_mat <- as.matrix(ref_htc_mat)
  # tmp <- colnames(ref_htc_mat)
  
  
  colnames(ref_htc_mat) <- seq(from=tad_start-1, length.out=ncol(ref_htc_mat), by=binSize)  # was 2709-2722, now is "54160000" "54420000" 
  stopifnot(colnames(ref_htc_mat) == as.character((as.numeric(tmp)-1)*binSize))
  
  stopifnot(as.character(colnames(ref_htc_mat)[ncol(ref_htc_mat)]) == as.character(tad_end-binSize))
  rownames(ref_htc_mat) <- seq(from=tad_start-1, length.out=ncol(ref_htc_mat), by=binSize)
  
  ### to emulate hic data, replace colnames of yuanlong
  stopifnot(as.character(rownames(ref_htc_mat)[ncol(ref_htc_mat)]) == as.character(tad_end-binSize))
  
  # ref_htc_mat[1,1] <- ref_htc_mat[4,4] <- NA
  
  plotTit <- paste0("Top ", i_top, " ", ref_hicds, ": ", tad_region, " ", tad_start, " - ", tad_end)
  
  
  ##################################### 1 ref only
  plotSubtit <- paste0(ref_hicds, " HiCDC pval [-log10] - ", binSize/1000, " kb (size=", nBins, ")")
  
  
  outFile <- file.path(outFolder, 
                       paste0("ref", ref_hicds, "top", i_top, "_", tad_region, "_", ref_hicds, "_hicdcPvalMat.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  my_plot_matrix(mat = ref_htc_mat, 
                 tad_coord = c(tad_start,tad_end),
                 resolution = binSize, 
                 bins_around=c(plot_around, plot_around),
                 # saveMatFile = matFile,
                 # saveMatcolFile =matcolFile,
                 # na.col = "green",
                 plotWithLeg=TRUE,
                 subTit=plotSubtit,
                 main=plotTit) 
  # mtext(side=3, text=plotSubtit)
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  ##################################### 2 match only
  
  match_htc_mat <- match_hicdcmat
  nBins <- bin_end - bin_start + 1
  stopifnot(nBins == ncol(match_htc_mat))
  
  colnames(match_htc_mat) <- seq(from=tad_start-1, length.out=ncol(match_htc_mat), by=binSize)  # was 2709-2722, now is "54160000" "54420000" 
  stopifnot(colnames(match_htc_mat) == as.character((as.numeric(tmp)-1)*binSize))
  
  stopifnot(as.character(colnames(match_htc_mat)[ncol(match_htc_mat)]) == as.character(tad_end-binSize))
  rownames(match_htc_mat) <- seq(from=tad_start-1, length.out=ncol(match_htc_mat), by=binSize)
  
  ### to emulate hic data, replace colnames of yuanlong
  stopifnot(as.character(rownames(match_htc_mat)[ncol(match_htc_mat)]) == as.character(tad_end-binSize))
  
  # match_htc_mat[1,1] <- match_htc_mat[4,4] <- NA
  plotTit <- paste0("Top ", i_top, " ", ref_hicds, ": ", tad_region, " ", tad_start, " - ", tad_end)

  # outFile <- file.path(#outFolder, 
  #                      paste0("match_htc_mat.Rdata"))
  # save(match_htc_mat, file=outFile)
  # 
  # outFile <- file.path(#outFolder,
  #                      paste0("ref_htc_mat.Rdata"))
  # save(ref_htc_mat, file=outFile)
  # 
  plotSubtit <- paste0(match_hicds, " HiCDC pval [-log10] - ", binSize/1000, " kb (size=", nBins, ")")
  
  outFile <- file.path(outFolder, 
                       paste0("ref", ref_hicds, "_top", i_top, "_", tad_region, "_", match_hicds, "_hicdcPvalMat.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  my_plot_matrix(mat = match_htc_mat, 
                 tad_coord = c(tad_start,tad_end),
                 resolution = binSize, 
                 bins_around=c(plot_around, plot_around),
                 # saveMatFile = matFile,
                 # saveMatcolFile =matcolFile,
                 # na.col = "green",
                 plotWithLeg=TRUE,
                 subTit=plotSubtit,
                 main=plotTit) 
  # mtext(side=3, text=plotSubtit)
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
  ##################################### 3 combined
  
  ##### combined
  
  stopifnot(colnames(ref_htc_mat) == colnames(match_htc_mat))
  
  combined_mat <- ref_htc_mat
  combined_mat[upper.tri(combined_mat, diag=F)] <- match_htc_mat[upper.tri(match_htc_mat)]
  diag(combined_mat) <- NA
  stopifnot(combined_mat[lower.tri(combined_mat, diag=F)] == ref_htc_mat[lower.tri(ref_htc_mat, diag=F)])
  stopifnot(combined_mat[upper.tri(combined_mat, diag=F)] == match_htc_mat[upper.tri(match_htc_mat, diag=F)])
  
  plotSubtit <- paste0("L: ", ref_hicds, "; U: ", match_hicds, " HiCDC pval [-log10] - ", binSize/1000, " kb (size=", nBins, ")")

  outFile <- file.path(outFolder, 
                       paste0("ref", ref_hicds, "_top",  i_top, "_ref", ref_hicds, "_", tad_region, "_combined_hicdcPvalMat.", plotType))
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  my_plot_matrix(mat = combined_mat, 
                 tad_coord = c(tad_start,tad_end),
                 resolution = binSize, 
                 bins_around=c(plot_around, plot_around),
                 # saveMatFile = matFile,
                 # saveMatcolFile =matcolFile,
                 na.col = "darkgrey",
                 checkSim=FALSE,
                 plotWithLeg=TRUE,
                 subTit=plotSubtit,
                 color = colorRampPalette(c("blue", "red"))(100), 
                 
                 legTit = paste0("HiCDC p-vals\n[-log10]"),
                 legCoords=  c(1,0,1.5,5),
                 
                 main=plotTit) 
  # mtext(side=3, text=plotSubtit)
  
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  outFile <- file.path(outFolder, 
                       paste0("ref", ref_hicds, "_top", i_top, "_", tad_region, "_", ref_hicds, "_hicdcPvalMat.Rdata"))
  save(ref_hicdcmat, file = outFile)
  cat(paste0("... ", outFile, "\n"))
  
  outFile <- file.path(outFolder, 
                       paste0("ref", ref_hicds, "_top", i_top, "_", tad_region, "_", match_hicds, "_hicdcPvalMat.Rdata"))
  save(match_hicdcmat, file = outFile)
  cat(paste0("... ", outFile, "\n"))
  
  
}
ds_dt$adjPvalComb_rd <- formatC(ds_dt$adjPvalComb, format = "e", digits = 2)
ds_dt$meanLogFC_rd <- round(ds_dt$meanLogFC,4)

rownames(ds_dt) <- 1:nrow(ds_dt)
write.table(ds_dt[1:nTop, c("hicds", "region", "region_genes", "adjPvalComb_rd")], sep="\t", row.names=T, quote=F)
write.table(ds_dt[1:nTop, c("hicds", "region","start", "end","meanLogFC_rd" ,"region_genes", "adjPvalComb_rd")], sep="\t", row.names=T, quote=F)