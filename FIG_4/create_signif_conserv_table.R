
# Rscript create_signif_conserv_table.R

library(stringr)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) <= 1)
if(length(args) == 1) {
  cmptype <- args[1]
  cmpPrefix <- paste0(cmptype, "_")
} else {
  cmptype <- ""
  cmpPrefix <- paste0(cmptype)
}


inFolder <- file.path("../../v2_Yuanlong_Cancer_HiC_data_TAD_DA/TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2", cmptype)

# x=get(load("GO_SIGNIF_ACROSS_HICDS_v2/conserved_signif_enrich.Rdata"))
# x2=get(load("GO_SIGNIF_ACROSS_HICDS_v2/go_all_conserved_signif_tads_genes.Rdata"))
x3=get(load(file.path(inFolder, "plot_matching_dt_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))
# x4=get(load(file.path(inFolder, "../../v2_Yuanlong_Cancer_HiC_data_TAD_DA/TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2/all_signif_matching_dt_adjPvalComb_0.01.Rdata"))
# x5=get(load(file.path(inFolder, "../../v2_Yuanlong_Cancer_HiC_data_TAD_DA/TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2/"))
# x6=get(load(file.path(inFolder, "/conserved_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata"))
# x7=get(load(file.path(inFolder, "signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata"))

dt <- get(load(file.path(inFolder, paste0(cmpPrefix, "conserved_regions_with_genes_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata"))))
# [1] "conserved_region"       "corresp_tads"           "intersect_genes_symbol"
# [4] "intersect_genes_entrez"

dt$nConserv <- sapply(dt$corresp_tads, function(x) {1+str_count(x, ",")})
dt <- dt[order(dt$nConserv, decreasing = TRUE),]
stopifnot(names(which.max(colSums(x3))) == dt$conserved_region[1])
stopifnot(max(colSums(x3)) == dt$nConserv[1])


outFolder <- file.path("CREATE_SIGNIF_CONSERV_TABLE", cmptype)
dir.create(outFolder, recursive = TRUE)
outFile <- file.path(outFolder, paste0(cmpPrefix, "conserv_regions_table.txt"))
write.table(dt, file=outFile, row.names=F, col.names=TRUE, sep="\t", quote=FALSE)
