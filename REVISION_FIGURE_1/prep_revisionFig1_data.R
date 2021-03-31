### CREATE TABLE WITH THE FOLLOWING COLUMNS
# - TAD ID
# - compartment
# - meanCorr
# - meanAggExpr
# - purityCorr

# Rscript prep_revisionFig1_data.R

outFolder <- file.path("PREP_REVISIONFIG1_DATA")
dir.create(outFolder, recursive = TRUE)

source("../full_dataset_names.R")

runFolder <- "../../v2_Yuanlong_Cancer_HiC_data_TAD_DA"

merging_cols <- c("hicds", "exprds", "region")

# outcols_order <- c(merging_cols, "binaryCptmtLab", "eightCptmtLab", "aggLog10Expr", "meanCorr",  "meanLogFC" ,"adjPvalComb", "purityCorr")
outcols_order <- c(merging_cols, "binaryCptmtLab", "eightCptmtLab", "aggLog10Expr", "meanCorr" ,"adjPvalComb", "purityCorr")

#### #### #### #### #### #### #### #### #### #### #### #### 
#### prep signif data
#### #### #### #### #### #### #### #### #### #### #### #### 
result_file <- file.path(runFolder,"CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData <- resultData[,c(merging_cols, "meanCorr", "meanLogFC", "adjPvalComb")]

#### #### #### #### #### #### #### #### #### #### #### #### 
#### prep purity data
#### #### #### #### #### #### #### #### #### #### #### #### 

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
agg_purity$hicds <- dirname(agg_purity$dataset) 
agg_purity$exprds <- basename(agg_purity$dataset) 
agg_purity$dataset <- NULL

merge_dt1 <- merge(agg_purity, resultData, by=c(merging_cols), all.x=TRUE, all.y=FALSE)  # !!! WILL DISCARD DATA WITHOUT PURITY SCORE !!!

stopifnot(nrow(merge_dt1) == nrow(agg_purity))

#### #### #### #### #### #### #### #### #### #### #### #### 
#### prep compartment data
#### #### #### #### #### #### #### #### #### #### #### #### 

inFolder <- file.path(runFolder, "REVISION_CHANGES_CPTMTLABELS_ALLDS")
outFile <- file.path(inFolder, "tad2cptmt_dt.Rdata")
tad2cptmt_dt <- get(load(outFile))
stopifnot(tad2cptmt_dt$startCptmtLabel == tad2cptmt_dt$tadCptmtLabel)
stopifnot(tad2cptmt_dt$startCptmtLabel == tad2cptmt_dt$endCptmtLabel)
stopifnot(substr(tad2cptmt_dt$tadCptmtLabel, start=1, stop=1) == tad2cptmt_dt$tad_binaryCptmtLab)

merge_dt2 <- merge(merge_dt1, tad2cptmt_dt[,c(merging_cols, c("tad_eightCptmtLab", "tad_binaryCptmtLab"))], 
                   by=c(merging_cols), all.x=TRUE, all.y=FALSE)  # !!! WILL DISCARD DATA WITHOUT PURITY SCORE !!!

stopifnot(nrow(merge_dt1) == nrow(merge_dt2))


#### #### #### #### #### #### #### #### #### #### #### #### 
#### prep expression level data
#### #### #### #### #### #### #### #### #### #### #### #### 

expr_var <- "aggLog10Expr"
all_exprLevel_dt <- get(load(file.path(runFolder, "REVISION_EXPRESSION_LEVEL", paste0(expr_var, "_aggByTAD_mean.Rdata"))))
all_exprLevel_dt$hicds <- dirname(dirname(all_exprLevel_dt$regionID))
all_exprLevel_dt$exprds <- basename(dirname(all_exprLevel_dt$regionID))
all_exprLevel_dt$region <- basename(all_exprLevel_dt$regionID)

merge_dt3 <- merge(merge_dt2, all_exprLevel_dt[,c(merging_cols, c("aggLog10Expr"))], 
                   by=c(merging_cols), all.x=TRUE, all.y=FALSE)  # !!! WILL DISCARD DATA WITHOUT PURITY SCORE !!!

stopifnot(nrow(merge_dt3) == nrow(merge_dt2))

#### #### #### #### #### #### #### #### #### #### #### #### 
#### change labels for output
#### #### #### #### #### #### #### #### #### #### #### #### 

colnames(merge_dt3) <- gsub("tad_", "", colnames(merge_dt3))

stopifnot(outcols_order %in% colnames(merge_dt3))

out_dt <- merge_dt3[,outcols_order]
outFile <- file.path(outFolder, "revision_fig1_cptmtAnnot_with_corr_purity.Rdata")
save(out_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile,"\n"))


stopifnot(merge_dt3$hicds %in% names(hicds_names))
merge_dt3$hicds <- hicds_names[paste0(merge_dt3$hicds)]

stopifnot(merge_dt3$exprds %in% names(exprds_names))
merge_dt3$exprds <- exprds_names[paste0(merge_dt3$exprds)]

outFile <- file.path(outFolder, "revision_fig1_cptmtAnnot_with_corr_purity.txt")
write.table(merge_dt3[,outcols_order], file=outFile, col.names=T, row.names=F, sep="\t", quote=F)
cat(paste0("... written: ", outFile,"\n"))


# out_dt <- merge_dt3
# out_dt$purityCorr <- round(out_dt$purityCorr, 4)
# out_dt$aggLog10Expr <- round(out_dt$aggLog10Expr, 4)
# out_dt$meanCorr <- round(out_dt$meanCorr, 4)
# out_dt$meanLogFC <- round(out_dt$meanLogFC, 4)
# out_dt$adjPvalComb <- round(out_dt$adjPvalComb, 4)
