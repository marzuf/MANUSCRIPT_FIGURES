require(foreach)
require(doMC)
registerDoMC(40)
require(ggsci)
require(ggpubr)
require(ggplot2)

# Rscript revision_changes_cptmtLabels_v3_notPF.R

### FOCUS ON THE B->B compartments

setDir <- "/media/electron"
setDir <- ""

runFolder <- "../../v2_Yuanlong_Cancer_HiC_data_TAD_DA"

plotType <- "svg"
myHeightGG <- 5
myWidthGG <- 7

source(file.path(runFolder, "revision_settings.R"))

outFolder <- "REVISION_CHANGES_CPTMTLABELS_V3_NOTPF"
dir.create(outFolder)

buildTable <- TRUE

hierarchyFolder <- file.path(setDir, 
"/mnt/ndata/Yuanlong/1.Projects/2.PROFILE/2.Results_review/2.Sup_tab_used_for_2nd_submission",
"Supplementary_Data_1_domain_hierarchies_127HiCmaps")
  
matchingFile <- file.path(runFolder, "REVISION_RANKDIFF_ACTIVDIFF", "matching_data.Rdata")

tadSignifThresh <- 0.01

# 
# head /mnt/ndata/Yuanlong/1.Projects/2.PROFILE/2.Results_review/2.Sup_tab_used_for_2nd_submission/Supplementary_Data_1_domain_hierarchies_127HiCmaps/LNCaP_prostate_cancer_binsize=40kb.bed
# #fields: chr, pos_start, pos_end, the full compartment label describing the position in the hierachy, normalized compartment domain rank, ., pos_star
# t, pos_end, color based on eight compartment model (A.1.1 - B.2.2), 1, eight compartment score (A.1.1 to B.2.2: 8 to 1 devided by 8)
# chr1    560001  600000  A.2.1.2.1.2.2.1.1.1.2.2 0.554371002132196       .       560001  600000  #FF9191 1       0.75
# chr1    720001  1080000 A.2.1.2.1.2.2.1.1.1.2.2 0.554371002132196       .       720001  1080000 #FF9191 1       0.75

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

cptmt_files <-  c(
  "LI_40kb"="LI_liver_binsize=40kb.bed",
    "GSE105381_HepG2_40kb"="HepG2_liver_cancer_binsize=40kb.bed",
"LG1_40kb" ="Lung_tissue_1_binsize=40kb.bed",
  "ENCSR444WCZ_A549_40kb"="A549_alveolar_basal_epithelial_cancer_binsize=40kb.bed",
"LG2_40kb"="Lung_tissue_2_binsize=40kb.bed",
  "ENCSR489OCU_NCI-H460_40kb"="NCI_H460_large_cell_lung_cancer_binsize=40kb.bed",
"GSE118514_RWPE1_40kb"="RWPE1_prostate_epithelium_binsize=40kb.bed",
"ENCSR346DCU_LNCaP_40kb"="LNCaP_prostate_cancer_binsize=40kb.bed",
   "GSE118514_22Rv1_40kb"="22Rv1_prostate_cancer_binsize=40kb.bed"
)


final_table_file <- file.path(runFolder, "CREATE_FINAL_TABLE/all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_dt <- get(load(final_table_file))
final_table_DT <- final_dt
final_table_DT$regionID <- file.path(final_table_DT$hicds, final_table_DT$exprds, final_table_DT$region)
stopifnot(!duplicated(final_table_DT$regionID))
regionID_pvals <- setNames(final_table_DT$adjPvalComb, final_table_DT$regionID)



#### #### #### #### #### #### #### #### #### #### #### #### 
#### retrieve the purity tagged tads



purity_ds <- "aran"
pm <- "CPE"
purity_plot_name <- paste0("Aran - ", pm)

corMet <- "pearson"
transfExpr <- "log10"
corrPurityQtThresh <- 0.05
signif_column <- "adjPvalComb"
signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", signifThresh)


purity_file <- file.path(runFolder,"ALLTADS_AND_PURITY_FINAL", purity_ds, pm, transfExpr, "all_ds_corrPurity_dt.Rdata") # here _final INPUT
purityData <- get(load(purity_file))
agg_purity <- aggregate(purityCorr~dataset+region, FUN=mean, data=purityData)

result_file <- file.path(runFolder,"CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)

merge_dt <- merge(agg_purity, resultData, by=c("dataset", "region"))  # !!! WILL DISCARD DATA WITHOUT PURITY SCORE !!!
merge_dt$region_ID <- file.path(merge_dt$dataset, merge_dt$region)
merge_dt$signif <- merge_dt$adjPvalComb <= signifThresh
merge_dt$signif_lab <- ifelse(merge_dt$signif, paste0("adj. p-val <=", signifThresh), paste0("adj. p-val >", signifThresh) )
purityCorrThresh <- as.numeric(quantile(merge_dt$purityCorr[!merge_dt$signif], probs = corrPurityQtThresh ))

# purity_flagged_tads <- merge_dt$region_ID[merge_dt$purityCorr <= purityCorrThresh]
tokeep_tads <- merge_dt$region_ID[merge_dt$purityCorr > purityCorrThresh]
cat(paste0( "purityCorrThresh = ", round(purityCorrThresh, 4), "\n"))
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 




stopifnot(names(cptmt_files) %in% final_dt$hicds)
# sub_final_dt <- final_dt[final_dt$hicds %in% names(cptmt_files),]
### CORRECTED 27.02 -> NEED TO KEEP ONLY THE normal vs. tumor !!!
sub_final_dt <- final_dt[final_dt$hicds %in% names(cptmt_files) &
                           final_dt$exprds %in% basename(all_pairs),]
stopifnot(nrow(sub_final_dt) > 0)

hicds = names(cptmt_files)[1]

if(buildTable) {
  tad2cptmt_final_dt <- foreach(hicds = names(cptmt_files), .combine='rbind')%do% {
    
    
    cptmt_dt <- read.delim(file.path(hierarchyFolder, cptmt_files[paste0(hicds)]), header=F, skip=1, stringsAsFactors = FALSE,
                           col.names=c( "chr", "pos_start", "pos_end", "cptmt_label", "normalized_rank",".",
                                        "pos_start", "pos_end", "color", "1", "cptmt_score"))
    stopifnot(cptmt_dt$pos_start.1 == cptmt_dt$pos_start)
    stopifnot(cptmt_dt$pos_end.1 == cptmt_dt$pos_end)
    
    hicds_final_dt <- sub_final_dt[sub_final_dt$hicds == hicds,]  
    stopifnot(nrow(hicds_final_dt) > 0)
    hicds_final_dt$chromo <- gsub("(chr.+)_TAD.+", "\\1", hicds_final_dt$region)
    stopifnot(hicds_final_dt$chromo %in% paste0("chr", 1:22))
    
    hicds_final_dt$midPos <- (hicds_final_dt$start+hicds_final_dt$end)/2
    
    hicds_final_dt$startCptmtLabel <- hicds_final_dt$startCptmtScore <- hicds_final_dt$startCptmtNormRank <- NA
    hicds_final_dt$endCptmtLabel <- hicds_final_dt$endCptmtScore <- hicds_final_dt$endCptmtNormRank <- NA
    hicds_final_dt$midCptmtLabel <- hicds_final_dt$midCptmtScore <- hicds_final_dt$midCptmtNormRank <- NA
    hicds_final_dt$tadCptmtLabel <- hicds_final_dt$tadCptmtScore <- hicds_final_dt$tadCptmtNormRank <- NA
    hicds_final_dt$i_startCptmt <- hicds_final_dt$i_endCptmt <- hicds_final_dt$i_midCptmt <- NA
    i=1
    i=106
    hicds_out_dt <- foreach(i = 1:nrow(hicds_final_dt), .combine='rbind') %dopar% {
      
      chromo <- hicds_final_dt$chromo[i]
      startPos <- hicds_final_dt$start[i]
      endPos <- hicds_final_dt$end[i]
      midPos <- hicds_final_dt$midPos[i]
      
      # !!! there might be gaps !!!
      # 875 chr1 248760001 248920000 B.1.1.2.2.2.2.2.1.1.2.2       0.3274232 .   248760001 248920000
      # 876 chr1 249040001 249240000 B.1.1.2.2.2.2.2.1.1.2.2       0.3274232 .   249040001 249240000
      
      
      # find where start, end, mid fall in
      
      startCptmt <- which(cptmt_dt$chr == chromo & 
                            cptmt_dt$pos_start <= startPos &
                            cptmt_dt$pos_end >= startPos)
      if(length(startCptmt) == 0) startCptmt <- "gap" 
      stopifnot(length(startCptmt) == 1)
      
      midCptmt <- which(cptmt_dt$chr == chromo & 
                          cptmt_dt$pos_start <= midPos &
                          cptmt_dt$pos_end >= midPos)
      if(length(midCptmt) == 0) midCptmt <- "gap" 
      stopifnot(length(midCptmt) == 1)
      
      endCptmt <- which(cptmt_dt$chr == chromo & 
                          cptmt_dt$pos_start <= endPos &
                          cptmt_dt$pos_end >= endPos)
      if(length(endCptmt) == 0) endCptmt <- "gap" 
      stopifnot(length(endCptmt) == 1)
      
      # if start cptmt = mid cptmt -> assign start
      # if end cptmt = mid cptmt -> assign end
      # assign mid otherwise
      tadCptmt <- ifelse(startCptmt == midCptmt, startCptmt,
                         ifelse(endCptmt == midCptmt, endCptmt,midCptmt))
      # for the start
      if(startCptmt == "gap" ) {
        hicds_final_dt$startCptmtLabel[i] <- hicds_final_dt$startCptmtScore[i] <- hicds_final_dt$startCptmtNormRank[i] <- "gap"
        hicds_final_dt$i_startCptmt[i] <- "gap"
      }else {
        hicds_final_dt$startCptmtLabel[i] <- cptmt_dt$cptmt_label[startCptmt] 
        hicds_final_dt$startCptmtScore[i] <- cptmt_dt$cptmt_score[startCptmt] 
        hicds_final_dt$startCptmtNormRank[i] <- cptmt_dt$normalized_rank[startCptmt]
        hicds_final_dt$i_startCptmt[i] <- startCptmt
      }
      # for the end
      if(endCptmt == "gap" ) {
        hicds_final_dt$endCptmtLabel[i] <- hicds_final_dt$endCptmtScore[i] <- hicds_final_dt$endCptmtNormRank[i] <- "gap"
        hicds_final_dt$i_endCptmt[i] <- "gap"
      } else {
        hicds_final_dt$endCptmtLabel[i] <-  cptmt_dt$cptmt_label[endCptmt] 
        hicds_final_dt$endCptmtScore[i] <- cptmt_dt$cptmt_score[endCptmt] 
        hicds_final_dt$endCptmtNormRank[i] <- cptmt_dt$normalized_rank[endCptmt]
        hicds_final_dt$i_endCptmt[i] <- endCptmt
      }
      # for the mid
      if(midCptmt == "gap" ) {
        hicds_final_dt$midCptmtLabel[i] <- hicds_final_dt$midCptmtScore[i] <- hicds_final_dt$midCptmtNormRank[i] <- "gap"
        hicds_final_dt$i_midCptmt[i] <- "gap"
      } else {
        hicds_final_dt$midCptmtLabel[i] <-  cptmt_dt$cptmt_label[midCptmt] 
        hicds_final_dt$midCptmtScore[i] <- cptmt_dt$cptmt_score[midCptmt] 
        hicds_final_dt$midCptmtNormRank[i] <- cptmt_dt$normalized_rank[midCptmt]
        hicds_final_dt$i_midCptmt[i] <- midCptmt
      }
      # for the TAD
      if(tadCptmt == "gap" ) {
        hicds_final_dt$tadCptmtLabel[i] <-  hicds_final_dt$tadCptmtScore[i] <- hicds_final_dt$tadCptmtNormRank[i] <-"gap"
      } else {
        hicds_final_dt$tadCptmtLabel[i] <-  cptmt_dt$cptmt_label[tadCptmt] 
        hicds_final_dt$tadCptmtScore[i] <- cptmt_dt$cptmt_score[tadCptmt] 
        hicds_final_dt$tadCptmtNormRank[i] <- cptmt_dt$normalized_rank[tadCptmt]
      }
      stopifnot(!is.na(hicds_final_dt[i,]))
      hicds_final_dt[i,]
    }
    hicds_out_dt
    
  }
  outFile <- file.path(outFolder, "tad2cptmt_final_dt.Rdata")
  save(tad2cptmt_final_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "tad2cptmt_final_dt.Rdata")
  tad2cptmt_final_dt <- get(load(outFile))
}
# tad2cptmt_final_dt=get(load("REVISION_CHANGES_CPTMTLABELS/tad2cptmt_final_dt.Rdata"))
mean(tad2cptmt_final_dt$tadCptmtLabel == tad2cptmt_final_dt$startCptmtLabel)
# [1] 0.9790391
mean(tad2cptmt_final_dt$tadCptmtLabel == tad2cptmt_final_dt$midCptmtLabel)
# [1] 1
mean(tad2cptmt_final_dt$midCptmtLabel == tad2cptmt_final_dt$startCptmtLabel)
# [1] 0.9790391
mean(tad2cptmt_final_dt$endCptmtLabel == tad2cptmt_final_dt$startCptmtLabel)
# [1] 0.9790391
mean(tad2cptmt_final_dt$endCptmtLabel == tad2cptmt_final_dt$startCptmtLabel)
# [1] 1
mean(tad2cptmt_final_dt$tadCptmtLabel == "gap")
# [1] 0.02096086
mean(tad2cptmt_final_dt$i_endCptmt == tad2cptmt_final_dt$i_startCptmt)
# [1] 0.9390172
sum(tad2cptmt_final_dt$startCptmtLabel == "gap")
# [1] 0
sum(tad2cptmt_final_dt$endCptmtLabel == "gap")
# [1] 0

stopifnot(sum(tad2cptmt_final_dt$startCptmtLabel == "gap") == 0)
stopifnot(sum(tad2cptmt_final_dt$endCptmtLabel == "gap") == 0)
stopifnot(tad2cptmt_final_dt$endCptmtLabel == tad2cptmt_final_dt$startCptmtLabel)

### UPDATE HERE THE TADS ARE THE COMPARTMENTS... ONLY MIDCPTMT MIGHT BE IN GAP
# retaint the start (= end ) compartment
tad2cptmt_final_dt$tadCptmtLabel <- tad2cptmt_final_dt$startCptmtLabel


# from A.1.1 to B.2.2
tad2cptmt_final_dt$tad_eightCptmtLab <- substr(tad2cptmt_final_dt$tadCptmtLabel, start=1, stop=5)

tad2cptmt_final_dt$tad_binaryCptmtLab <- substr(tad2cptmt_final_dt$tadCptmtLabel, start=1, stop=1)
tad2cptmt_final_dt$region_ID <- file.path(tad2cptmt_final_dt$hicds, tad2cptmt_final_dt$exprds, tad2cptmt_final_dt$region)
stopifnot(!duplicated(tad2cptmt_final_dt$region_ID))
tad2cptmts <- setNames(tad2cptmt_final_dt$tad_binaryCptmtLab, tad2cptmt_final_dt$region_ID)
tad2cptmts_full <- setNames(tad2cptmt_final_dt$tadCptmtLabel, tad2cptmt_final_dt$region_ID)
tad2cptmts_eight <- setNames(tad2cptmt_final_dt$tad_eightCptmtLab, tad2cptmt_final_dt$region_ID)

# remove the TADs in gaps
tad2cptmt_dt <- tad2cptmt_final_dt[tad2cptmt_final_dt$tadCptmtLabel != "gap",]
stopifnot(tad2cptmt_dt$tad_binaryCptmtLab %in% c("A", "B"))
### UPDATE HERE THE TADS ARE THE COMPARTMENTS... ONLY MIDCPTMT MIGHT BE IN GAP
# START CPTMT = ENDCPTMT, AND NEITHER START NOR END ARE GAP


######################## PLOT ALL DATASETS
# this can be done for all datasets -> cptmt assignment

all_ds_dt <- tad2cptmt_final_dt
all_ds_dt$tadSignif <- ifelse(all_ds_dt$adjPvalComb <= tadSignifThresh, "signif.", "not signif.")

my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], unique(all_ds_dt$tadSignif))
legTitle <- ""

stopifnot(!duplicated(all_ds_dt$region_ID))

all_cmps <- unique(file.path(all_ds_dt$hicds, all_ds_dt$exprds))

mySub <- paste0("# DS = ", length(all_cmps), "; # TADs = ", nrow(all_ds_dt), 
                " (signif.: ", sum(all_ds_dt$adjPvalComb <= tadSignifThresh), ")")

all_ds_dt$region_cptmt <- all_ds_dt$tad_binaryCptmtLab
all_ds_dt$region_cptmt_eight <- all_ds_dt$tad_eightCptmtLab














# matching_data <- get(load("REVISION_RANKDIFF_ACTIVDIFF/matching_data.Rdata"))
matching_data <- get(load(file=matchingFile))

ds1_matching_dt <- do.call(rbind, lapply(matching_data, function(x)x[["norm_matching_pval_tadRank_dt"]]))
unique(ds1_matching_dt$ref_hicds)
# [1] "LI_40kb"              "LG1_40kb"             "LG2_40kb"             "GSE118514_RWPE1_40kb"
ds2_matching_dt <- do.call(rbind, lapply(matching_data, function(x)x[["tumor_matching_pval_tadRank_dt"]]))
unique(ds2_matching_dt$ref_hicds)
# [1] "GSE105381_HepG2_40kb"      "ENCSR444WCZ_A549_40kb"     "ENCSR489OCU_NCI-H460_40kb"
# [4] "ENCSR346DCU_LNCaP_40kb"    "GSE118514_22Rv1_40kb"     


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

stopifnot(matching_withRank_dt$ref_region_ID %in% names(tad2cptmts_full))
matching_withRank_dt$ref_region_cptmt_full <- tad2cptmts_full[paste0(matching_withRank_dt$ref_region_ID)]
stopifnot(!is.na(matching_withRank_dt$ref_region_cptmt_full))

stopifnot(matching_withRank_dt$ref_region_ID %in% names(tad2cptmts_eight))
matching_withRank_dt$ref_region_cptmt_eight <- tad2cptmts_eight[paste0(matching_withRank_dt$ref_region_ID)]
stopifnot(!is.na(matching_withRank_dt$ref_region_cptmt_eight))


stopifnot(matching_withRank_dt$matching_region_ID %in% names(regionID_pvals))
matching_withRank_dt$matching_region_pval <- regionID_pvals[paste0(matching_withRank_dt$matching_region_ID)]
stopifnot(!is.na(matching_withRank_dt$matching_region_pval))

stopifnot(matching_withRank_dt$matching_region_ID %in% names(tad2cptmts))
matching_withRank_dt$matching_region_cptmt <- tad2cptmts[paste0(matching_withRank_dt$matching_region_ID)]
stopifnot(!is.na(matching_withRank_dt$matching_region_cptmt))

stopifnot(matching_withRank_dt$matching_region_ID %in% names(tad2cptmts_full))
matching_withRank_dt$matching_region_cptmt_full <- tad2cptmts_full[paste0(matching_withRank_dt$matching_region_ID)]
stopifnot(!is.na(matching_withRank_dt$matching_region_cptmt_full))

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
matching_withRank_dt$norm2tumor_cptmtChange_binary <- paste0(matching_withRank_dt$normCptmt_binary, "->", matching_withRank_dt$tumorCptmt_binary)

################################################################## RANK DIFF AND COMPARTMENT CHANGE
matching_withRank_dt$sameCptmt_binary <- matching_withRank_dt$ref_region_cptmt == matching_withRank_dt$matching_region_cptmt
matching_withRank_dt$sameCptmt_eight <- matching_withRank_dt$ref_region_cptmt_eight == matching_withRank_dt$matching_region_cptmt_eight

save(matching_withRank_dt, file="tmp_dt.Rdata", version=2)


######################## ######################## ######################## 

# filter here notPF TADs

matching_withRank_dt <- matching_withRank_dt[matching_withRank_dt$ref_region_ID %in% tokeep_tads &
                                               matching_withRank_dt$matching_region_ID %in% tokeep_tads,]

stopifnot(nrow(matching_withRank_dt) > 0)

######################## 


noGap_dt <- matching_withRank_dt[!grepl("g", matching_withRank_dt$norm2tumor_cptmtChange),]
stopifnot(nrow(noGap_dt) == nrow(matching_withRank_dt))

matching_withRank_dt$abs_rankDiff <- abs(matching_withRank_dt$rankDiff)
summary(matching_withRank_dt$abs_rankDiff[])

vartype="eight"
for(vartype in c("binary", "eight")) {
  
  legTitle <- paste0("sameCptmt_", vartype, "\n(same=",sum(matching_withRank_dt[,paste0("sameCptmt_", vartype)]) , ")")
  my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], sort(unique(matching_withRank_dt[,paste0("sameCptmt_", vartype)])))
  
  plotTit <- paste0("TAD rank diff. and cptmt changes (notPF)")
  
  all_cmps <- unique(file.path(matching_withRank_dt$matching_hicds, matching_withRank_dt$matching_exprds,
                               matching_withRank_dt$ref_hicds, matching_withRank_dt$ref_exprds))
  
  mySub <- paste0("# DS comparisons = ", length(all_cmps), "; # TADs = ", nrow(matching_withRank_dt), 
                  " (signif.: ", sum(matching_withRank_dt$adjPval <= tadSignifThresh), ")")
  
  p3 <- gghistogram(matching_withRank_dt,
                  x = paste0("rankDiff"),
                  # y = "..count..",
                  y = "..density..",
                  bins=100,
                  # combine = TRUE,                  # Combine the 3 plots
                  xlab = paste0("Rank diff. with matched TAD"),
                  # ylab = "# of TADs",
                  ylab="Density",
                  # add = "median",                  # Add median line.
                  rug = FALSE,                      # Add marginal rug
                  color = paste0("sameCptmt_", vartype),
                  fill = paste0("sameCptmt_", vartype),
                  palette = "jco"
  ) +
    ggtitle(plotTit, subtitle = mySub)+
    scale_color_manual(values=my_cols)+
    scale_fill_manual(values=my_cols)  +
    labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
    guides(color=FALSE)+
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
       mytheme
  
  outFile <- file.path(outFolder, paste0("tad_rankDiff_same_diff_cptmt_", vartype, "_density.", plotType))
  ggsave(p3, file=outFile, height=myHeightGG, width=myWidthGG)
  cat(paste0("... written: ", outFile, "\n"))
  
  print(paste0(vartype, " - abs. diff. - diff. cptmt"))
  print(summary(matching_withRank_dt$abs_rankDiff[!matching_withRank_dt[,paste0("sameCptmt_", vartype)]]))
  
  
}

################################################################## PLOT THE B->B
plot_types <- c("withGap", "noGap")

plot_types="noGap"

for(p_type in plot_types) {
  if(p_type == "withGap") {
    plot_dt <-   matching_withRank_dt
  } else if(p_type == "noGap"){
    plot_dt <- noGap_dt
  } else {
    stop("error\n")
  }
  save(plot_dt, file=file.path(outFolder, "plot_dt.Rdata"), version=2)
  cat(paste0("... written: ", file.path(outFolder, "plot_dt.Rdata"), "\n"))
  
  ### >> for this version, plot only the B->B changes
  # plot_dt <- plot_dt[plot_dt$norm2tumor_cptmtChange == "B->B",] 
  plot_dt <- plot_dt ## plot here all for PF !!!
  stopifnot(nrow(plot_dt) > 0)
  init_plot_dt <- plot_dt
  
  all_suffices <- c("_binary", "_eight")
    all_suffices <- c("", "_eight")
   all_suffices <- c( "_eight")
  
  for(suffix in all_suffices) {
    
  
    for(bplot in c("all", "grouped")) {
      
      if(bplot == "grouped") {
        plot_dt <- init_plot_dt
        plot_dt[,paste0("norm2tumor_cptmtChange", suffix)][plot_dt[,paste0("tumorCptmt", suffix)] == plot_dt[,paste0("normCptmt", suffix)]] <- "same"
      } else {
        plot_dt <- init_plot_dt
      }
      
      my_cols <- setNames(pal_jama()(5)[c(3, 2,4)], sort(unique(plot_dt$ref_tadSignif)))
      legTitle <- ""
      
      all_cmps <- unique(file.path(plot_dt$matching_hicds, plot_dt$matching_exprds,
                                   plot_dt$ref_hicds, plot_dt$ref_exprds))
      
      mySub <- paste0("# DS comparisons = ", length(all_cmps), "; # TADs = ", nrow(plot_dt), 
                      " (signif.: ", sum(plot_dt$adjPval <= tadSignifThresh), ")")
      
      ################# BARPLOT DISTRIBUTION CPTMT CHANGE CATEGORY SIGNIF/NOT SIGNIF - BINARY
      
      # barplot: ratio of b->a ou a->b changes in signif. and not signif
      agg_dt <- aggregate(as.formula(paste0("ref_region_ID~norm2tumor_cptmtChange", suffix, " + ref_tadSignif")),
                          data=plot_dt, FUN=length)
      colnames(agg_dt)[colnames(agg_dt)=="ref_region_ID"] <- "nTADs"
      tmp_dt <- aggregate(ref_region_ID~ref_tadSignif, data=plot_dt, FUN=length)
      nTot <- setNames(tmp_dt$ref_region_ID, tmp_dt$ref_tadSignif)
      agg_dt$ratioTADs <- agg_dt$nTADs/nTot[paste0(agg_dt$ref_tadSignif)]
      plotTit <- paste0("norm vs. tumor ratio TADs by cptmt", suffix, " change (notPF; ", p_type, ")")
      agg_dt$norm2tumor_cptmtChange <- factor(agg_dt[,paste0("norm2tumor_cptmtChange", suffix)], levels=rev(c("A->A", "B->B", "A->B", "B->A")))
      stopifnot(!is.na(agg_dt[,paste0("norm2tumor_cptmtChange", suffix)]))
      ggbar_p <-  ggbarplot(agg_dt, 
                            y="ratioTADs",
                            x="ref_tadSignif", 
                            fill=paste0("norm2tumor_cptmtChange", suffix)) +
        ggtitle(plotTit, subtitle=mySub)+
        mytheme +
        labs(x="" , y ="ratio of TADs", color=paste0(legTitle),fill=paste0(legTitle)) + 
        theme(
          axis.text.x = element_text(hjust=1, vjust=0.5,size=10,angle=90)
        )
      outFile <- file.path(outFolder,paste0("norm_vs_tumor_cptmt", suffix, "_change_signif_notSignif_", p_type, "_", bplot, "_barplot.", plotType))
      ggsave(ggbar_p, filename = outFile, height=myHeightGG*1.2, width=myWidthGG)
      cat(paste0("... written: ", outFile, "\n"))
      
      ################# BOXPLOT FC BY CPT CHANGE - BINARY 
      plotTit <- paste0("norm vs. tumor FC change and cptmt", suffix, " change (notPF; ", p_type, ")")
      
      ggbox_p <- ggboxplot(
        data=plot_dt,
        xlab="",
        color = "ref_tadSignif",
        x=paste0("norm2tumor_cptmtChange", suffix), 
        y="tumorMinusNorm_meanLog2FC"
      ) + mytheme +
        ggtitle(plotTit, subtitle = mySub)+
        scale_color_manual(values=my_cols)+
        scale_fill_manual(values=my_cols)  +
        labs(color=paste0(legTitle),fill=paste0(legTitle), x="", y="(tumor-norm) meanLog2FC") +
        # guides(color=FALSE)+
        scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
      
      outFile <- file.path(outFolder,paste0("norm_vs_tumor_meanLog2FC_cptmt", suffix, "_change_", p_type, "_", bplot, "_boxplot.", plotType))
      ggsave(ggbox_p, filename = outFile, height=myHeightGG, width=myWidthGG)
      cat(paste0("... written: ", outFile, "\n"))
      
      
      
    }
  }
}


# ggboxplot(
#   data=plot_dt,
#   x="norm2tumor_cptmtChange", 
#   y="tumorMinusNorm_meanLog2FC"
# ) +
#   mytheme +
#   ggtitle(plotTit, subtitle = mySub)+
#   scale_color_manual(values=my_cols)+
#   scale_fill_manual(values=my_cols)  +
#   labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
#   # guides(color=FALSE)+
#   scale_y_continuous(breaks = scales::pretty_breaks(n = 10))




# subA_dt <- noGap_dt[grepl("A->", noGap_dt$norm2tumor_cptmtChange),]
# 
# subB_dt <- noGap_dt[grepl("B->", noGap_dt$norm2tumor_cptmtChange),]
# 
# ggdensity(subA_dt,
#           x = paste0("tumorMinusNorm_meanLog2FC"),
#           y = "..density..",
#           # combine = TRUE,                  # Combine the 3 plots
#           xlab = paste0("Rank diff. with matched TAD"),
#           # add = "median",                  # Add median line.
#           rug = FALSE,                      # Add marginal rug
#           color = "norm2tumor_cptmtChange",
#           fill = "norm2tumor_cptmtChange",
#           palette = "jco"
# ) 
# 
# ggdensity(subB_dt,
#           x = paste0("tumorMinusNorm_meanLog2FC"),
#           y = "..density..",
#           # combine = TRUE,                  # Combine the 3 plots
#           xlab = paste0("Rank diff. with matched TAD"),
#           # add = "median",                  # Add median line.
#           rug = FALSE,                      # Add marginal rug
#           color = "norm2tumor_cptmtChange",
#           fill = "norm2tumor_cptmtChange",
#           palette = "jco"
# ) 
# 
# 
# ggdensity(noGap_dt,
#           x = paste0("tumorMinusNorm_meanLog2FC"),
#           y = "..density..",
#           # combine = TRUE,                  # Combine the 3 plots
#           xlab = paste0("Rank diff. with matched TAD"),
#           # add = "median",                  # Add median line.
#           rug = FALSE,                      # Add marginal rug
#           color = "norm2tumor_cptmtChange",
#           fill = "norm2tumor_cptmtChange",
#           palette = "jco"
# ) 

# same with difference signif. not signif.





# matching_withRank_dt[
#   matching_withRank_dt$matching_region_ID == "LI_40kb/TCGAlihc_norm_lihc/chr1_TAD12" | 
#     matching_withRank_dt$ref_region_ID == "LI_40kb/TCGAlihc_norm_lihc/chr1_TAD12",]
# 
# matchingID_maxOverlapBp      refID            ref_hicds         ref_exprds
# 1                  chr1_TAD10 chr1_TAD12              LI_40kb TCGAlihc_norm_lihc
# 17046              chr1_TAD12 chr1_TAD10 GSE105381_HepG2_40kb TCGAlihc_norm_lihc
# matching_hicds    matching_exprds   adjPval chromo.x start.x   end.x refID_rank
# 1     GSE105381_HepG2_40kb TCGAlihc_norm_lihc 0.4143104     chr1 3400001 3760000  0.6879433
# 17046              LI_40kb TCGAlihc_norm_lihc 0.3482057     chr1 3520001 3800000  0.8377029
# normMeanFC chromo.y start.y   end.y matchingID_rank tumorMeanFC   rankDiff
# 1     -0.2967625     chr1 3520001 3800000       0.8377029  -0.2442317 -0.1497596
# 17046 -0.2967625     chr1 3400001 3760000       0.6879433  -0.2442317  0.1497596
# ref_region_ID
# 1                  LI_40kb/TCGAlihc_norm_lihc/chr1_TAD12
# 17046 GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/chr1_TAD10
# matching_region_ID ref_region_pval matching_region_pval
# 1     GSE105381_HepG2_40kb/TCGAlihc_norm_lihc/chr1_TAD10       0.4143104            0.3482057
# 17046              LI_40kb/TCGAlihc_norm_lihc/chr1_TAD12       0.3482057            0.4143104
# ref_tadSignif
# 1       not signif.
# 17046   not signif.
# # normMeanFC always the same -> correct ! 




