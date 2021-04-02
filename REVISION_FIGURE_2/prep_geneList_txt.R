# Rscript prep_geneList_txt.R K562_40kb TCGAlaml_wt_mutFLT3
# Rscript prep_geneList_txt.R GSE118514_RWPE1_40kb TCGAprad_norm_prad

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 2)

hicds <- args[1]
exprds <- args[2]

inFile <- file.path(file.path("/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA",
                              "PIPELINE/OUTPUT_FOLDER",
                              hicds, exprds,
                              "0_prepGeneData", "pipeline_geneList.Rdata"))

outFile <- file.path(file.path("/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA",
                              "PIPELINE/OUTPUT_FOLDER",
                              hicds, exprds,
                              "0_prepGeneData", "pipeline_geneList.txt"))

inValues <- get(load(inFile))

out_dt <- data.frame(
  name = names(inValues),
  value = as.character(inValues),
  stringsAsFactors = FALSE
)

write.table(out_dt, file = outFile, col.names=T, row.names=F, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))
