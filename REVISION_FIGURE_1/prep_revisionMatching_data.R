
# Rscript prep_revisionMatching_data.R

outFolder <- "PREP_REVISIONMATCHING_DATA"
dir.create(outFolder, recursive=TRUE)

runFolder <- "../../v2_Yuanlong_Cancer_HiC_data_TAD_DA"

all_pairs <- c(
  file.path("LI_40kb","GSE105381_HepG2_40kb", "TCGAlihc_norm_lihc"),
  file.path("LG1_40kb" ,"ENCSR444WCZ_A549_40kb", "TCGAluad_norm_luad"),
  file.path("LG2_40kb" ,"ENCSR444WCZ_A549_40kb" ,"TCGAluad_norm_luad"),
  file.path("LG1_40kb", "ENCSR489OCU_NCI-H460_40kb", "TCGAluad_norm_luad"), 
  file.path("LG2_40kb", "ENCSR489OCU_NCI-H460_40kb", "TCGAluad_norm_luad"), 
  file.path("GSE118514_RWPE1_40kb", "ENCSR346DCU_LNCaP_40kb", "TCGAprad_norm_prad"),
  file.path("GSE118514_RWPE1_40kb", "GSE118514_22Rv1_40kb", "TCGAprad_norm_prad")
)
all_normal_ds <- as.character(sapply(all_pairs, function(x) dirname(dirname(x))))
all_tumor_ds <-  as.character(sapply(all_pairs, function(x) basename(dirname(x))))


tadSignifThresh <- 0.01

matchingFile <- file.path(runFolder, "REVISION_RANKDIFF_ACTIVDIFF", "matching_data.Rdata")
matching_data <- get(load(file=matchingFile))

ds1_matching_dt <- do.call(rbind, lapply(matching_data, function(x)x[["norm_matching_pval_tadRank_dt"]]))
unique(ds1_matching_dt$ref_hicds)
# [1] "LI_40kb"              "LG1_40kb"             "LG2_40kb"             "GSE118514_RWPE1_40kb"
ds2_matching_dt <- do.call(rbind, lapply(matching_data, function(x)x[["tumor_matching_pval_tadRank_dt"]]))
unique(ds2_matching_dt$ref_hicds)
# [1] "GSE105381_HepG2_40kb"      "ENCSR444WCZ_A549_40kb"     "ENCSR489OCU_NCI-H460_40kb"
# [4] "ENCSR346DCU_LNCaP_40kb"    "GSE118514_22Rv1_40kb"     



inFolder <- file.path(runFolder, "REVISION_CHANGES_CPTMTLABELS_ALLDS")
outFile <- file.path(inFolder, "tad2cptmt_dt.Rdata")
all_tad2cptmt_dt <- get(load(outFile))
#all_tad2cptmt_dt <- all_tad2cptmt_dt[all_tad2cptmt_dt$region_ID %in% tokeep_tads,]

tad2cptmts_eight <- setNames(all_tad2cptmt_dt$tad_eightCptmtLab, all_tad2cptmt_dt$region_ID)
tad2cptmts <- setNames(all_tad2cptmt_dt$tad_binaryCptmtLab, all_tad2cptmt_dt$region_ID)
regionID_pvals <- setNames(all_tad2cptmt_dt$adjPvalComb, all_tad2cptmt_dt$region_ID)

matching_withRank_dt <- rbind(ds1_matching_dt, ds2_matching_dt)
rownames(matching_withRank_dt) <- NULL
matching_withRank_dt$ref_chromo <- gsub("(chr.+)_.+", "\\1", matching_withRank_dt$refID )
matching_withRank_dt$matching_chromo <- gsub("(chr.+)_.+", "\\1", matching_withRank_dt$matchingID_maxOverlapBp )
stopifnot( matching_withRank_dt$matching_chromo==matching_withRank_dt$ref_chromo )

stopifnot(matching_withRank_dt$matching_exprds == matching_withRank_dt$ref_exprds )

matching_withRank_dt$ref_region_ID <- file.path(matching_withRank_dt$ref_hicds,
                                                matching_withRank_dt$ref_exprds,
                                                matching_withRank_dt$refID)

matching_withRank_dt$matching_region_ID <- file.path(matching_withRank_dt$matching_hicds,
                                                     matching_withRank_dt$matching_exprds,
                                                     matching_withRank_dt$matchingID_maxOverlapBp)

stopifnot(matching_withRank_dt$ref_region_ID %in% names(regionID_pvals))
matching_withRank_dt$ref_region_pval <- regionID_pvals[paste0(matching_withRank_dt$ref_region_ID)]
stopifnot(!is.na(matching_withRank_dt$ref_region_pval))
stopifnot(round(matching_withRank_dt$ref_region_pval, 6) == round(matching_withRank_dt$adjPval, 6))

stopifnot(matching_withRank_dt$ref_region_ID %in% names(tad2cptmts))
matching_withRank_dt$ref_region_cptmt <- tad2cptmts[paste0(matching_withRank_dt$ref_region_ID)]
stopifnot(!is.na(matching_withRank_dt$ref_region_cptmt))



stopifnot(matching_withRank_dt$ref_region_ID %in% names(tad2cptmts_eight))
matching_withRank_dt$ref_region_cptmt_eight <- tad2cptmts_eight[paste0(matching_withRank_dt$ref_region_ID)]
stopifnot(!is.na(matching_withRank_dt$ref_region_cptmt_eight))


stopifnot(matching_withRank_dt$matching_region_ID %in% names(regionID_pvals))
matching_withRank_dt$matching_region_pval <- regionID_pvals[paste0(matching_withRank_dt$matching_region_ID)]
stopifnot(!is.na(matching_withRank_dt$matching_region_pval))

stopifnot(matching_withRank_dt$matching_region_ID %in% names(tad2cptmts))
matching_withRank_dt$matching_region_cptmt <- tad2cptmts[paste0(matching_withRank_dt$matching_region_ID)]
stopifnot(!is.na(matching_withRank_dt$matching_region_cptmt))




stopifnot(matching_withRank_dt$matching_region_ID %in% names(tad2cptmts_eight))
matching_withRank_dt$matching_region_cptmt_eight <- tad2cptmts_eight[paste0(matching_withRank_dt$matching_region_ID)]
stopifnot(!is.na(matching_withRank_dt$matching_region_cptmt_eight))

matching_withRank_dt$ref_tadSignif <- ifelse(matching_withRank_dt$adjPval <= tadSignifThresh, "signif.", "not signif.")

stopifnot(matching_withRank_dt$ref_hicds %in% c(all_tumor_ds, all_normal_ds))
stopifnot(matching_withRank_dt$matching_hicds %in% c(all_tumor_ds, all_normal_ds))

stopifnot(0 == sum(matching_withRank_dt$ref_hicds %in% all_tumor_ds & matching_withRank_dt$matching_hicds %in% all_tumor_ds))
stopifnot(0 == sum(matching_withRank_dt$ref_hicds %in% all_normal_ds & matching_withRank_dt$matching_hicds %in% all_normal_ds))

matching_withRank_dt$normCptmt <- ifelse(matching_withRank_dt$ref_hicds %in% all_normal_ds, matching_withRank_dt$ref_region_cptmt,
                                         ifelse(matching_withRank_dt$matching_hicds %in% all_normal_ds, matching_withRank_dt$matching_region_cptmt,NA))
stopifnot(!is.na(matching_withRank_dt$normCptmt))

matching_withRank_dt$tumorCptmt <- ifelse(matching_withRank_dt$ref_hicds %in% all_tumor_ds, matching_withRank_dt$ref_region_cptmt,
                                         ifelse(matching_withRank_dt$matching_hicds %in% all_tumor_ds, matching_withRank_dt$matching_region_cptmt,NA))
stopifnot(!is.na(matching_withRank_dt$tumorCptmt))


matching_withRank_dt$normCptmt_eight <- ifelse(matching_withRank_dt$ref_hicds %in% all_normal_ds, matching_withRank_dt$ref_region_cptmt_eight,
                                         ifelse(matching_withRank_dt$matching_hicds %in% all_normal_ds, matching_withRank_dt$matching_region_cptmt_eight,NA))
stopifnot(!is.na(matching_withRank_dt$normCptmt_eight))

matching_withRank_dt$tumorCptmt_eight <- ifelse(matching_withRank_dt$ref_hicds %in% all_tumor_ds, matching_withRank_dt$ref_region_cptmt_eight,
                                          ifelse(matching_withRank_dt$matching_hicds %in% all_tumor_ds, matching_withRank_dt$matching_region_cptmt_eight,NA))
stopifnot(!is.na(matching_withRank_dt$tumorCptmt_eight))




matching_withRank_dt$norm2tumor_cptmtChange <- paste0(matching_withRank_dt$normCptmt, "->", matching_withRank_dt$tumorCptmt)
matching_withRank_dt$tumorMinusNorm_meanLog2FC <- matching_withRank_dt$tumorMeanFC - matching_withRank_dt$normMeanFC

matching_withRank_dt$norm2tumor_cptmtChange_eight <- paste0(matching_withRank_dt$normCptmt_eight, "->", matching_withRank_dt$tumorCptmt_eight)
matching_withRank_dt$norm2tumor_cptmtChange <- paste0(matching_withRank_dt$normCptmt, "->", matching_withRank_dt$tumorCptmt)

################################################################## RANK DIFF AND COMPARTMENT CHANGE
matching_withRank_dt$sameCptmt <- matching_withRank_dt$ref_region_cptmt == matching_withRank_dt$matching_region_cptmt
matching_withRank_dt$sameCptmt_eight <- matching_withRank_dt$ref_region_cptmt_eight == matching_withRank_dt$matching_region_cptmt_eight

outFile <- file.path(outFolder, "matching_withRank_dt.Rdata")
save(matching_withRank_dt, file=outFile, version=2)
