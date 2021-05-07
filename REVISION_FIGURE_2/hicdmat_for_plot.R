setDir <- "/media/electron"
setDir <- ""

options(scipen = 999)


# Rscript hicdcmat_for_plot.R

outFolder <- "HICDCMAT_FOR_PLOT"
dir.create(outFolder, recursive = TRUE)


binSize <- 20000

ref_hicds <- "RWPE1"


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
# tads_to_plot <- list(
#   c("chr12", 54160001, 54440000)
# )
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
  
  matFile <- file.path(outFolder, paste0("ref", ref_hicds, "_", tad_chromo, "_", tad_start, "_", tad_end, "_", binSize/1000, 
                                         "kb_22Rv1_hicdc_mat.Rdata"))
  
  hicdcmat_22rv1 <- dt[[i]][["hic_dc_mat"]][["22Rv1"]]
  stopifnot(isSymmetric(as.matrix(hicdcmat_22rv1)))
  save(hicdcmat_22rv1, file = matFile)
  cat(paste0("... written:" , matFile, "\n"))
  
  matFile <- file.path(outFolder, paste0("ref", ref_hicds, "_", tad_chromo, "_", tad_start, "_", tad_end, "_", binSize/1000, 
                                         "kb_RWPE1_hicdc_mat.Rdata"))
  hicdcmat_RWPE1 <- dt[[i]][["hic_dc_mat"]][["RWPE1"]]
  stopifnot(isSymmetric(as.matrix(hicdcmat_RWPE1)))
  save(hicdcmat_RWPE1, file = matFile)
  
  cat(paste0("... written:" , matFile, "\n"))
}
