outFolder <- file.path("PREP_PRAD_DATA")
dir.create(outFolder,recursive = TRUE)

# Rscript prep_prad_data.R
source("../full_dataset_names.R")

keepCols <- c("hicds", "exprds", "region", "start", "end", "region_genes", "meanLogFC", "meanCorr", "adjPvalComb",
              "tad_eightCptmtLab", "tadCptmtNormRank",  "tadCptmtScore"  )

pvalthresh <- 0.01

runFolder <- "../../v2_Yuanlong_Cancer_HiC_data_TAD_DA"
result_file <- file.path(runFolder,"CREATE_FINAL_TABLE", "all_result_dt.Rdata")
resultData <- get(load(result_file))
resultData$dataset <- file.path(resultData$hicds, resultData$exprds)
resultData$region_id <- file.path(resultData$dataset, resultData$region)
resultData$signif_lab <- ifelse(resultData$adjPvalComb <= pvalthresh, "signif", "notsignif")
resultData$direction_lab <- ifelse(resultData$meanLogFC <0, "down", 
                                   ifelse(resultData$meanLogFC >0, "up",NA))
stopifnot(!is.na(resultData$direction_lab))
resultData$signif_lab <- paste0(resultData$signif_lab, ".", resultData$direction_lab)
signif_labs <- setNames(resultData$signif_lab, resultData$region_id)

hicds_of_interest <- c( "GSE118514_RWPE1_40kb",  "GSE118514_22Rv1_40kb")
# hicds_of_interest <- c( "GSE118514_RWPE1_40kb")
exprds_of_interest <- c("TCGAprad_norm_prad")

# sum(resultData$hicds=="GSE118514_RWPE1_40kb")
# [1] 1758
# sum(resultData$hicds=="GSE118514_22Rv1_40kb")
# [1] 1743

dt <- resultData[resultData$hicds %in% hicds_of_interest & resultData$exprds %in% exprds_of_interest ,]
dt <- dt[order(dt$hicds, dt$exprds, dt$adjPvalComb),]


cptmt_dt <- get(load(file.path(runFolder, "REVISION_EXPRESSIONLEVEL_CPTMTS/tad2cptmt_dt.Rdata")))
stopifnot(cptmt_dt$startCptmtLabel==cptmt_dt$tadCptmtLabel)
head(cptmt_dt)

out_dt <- merge(dt, cptmt_dt[,c("hicds", "exprds", "region", "tad_eightCptmtLab", "tadCptmtNormRank","tadCptmtScore")],
                by=c("hicds", "exprds", "region"),all.x=T,all.y=F)
out_dt <- out_dt[order(out_dt$hicds, out_dt$exprds, out_dt$adjPvalComb),]
rownames(out_dt) <- NULL
stopifnot(!is.na(out_dt))

change_cptmt_dt <- get(load(file.path(runFolder,"REVISION_CHANGES_CPTMTLABELS_V2/plot_dt.Rdata")))
change_cptmt_dt <- change_cptmt_dt[change_cptmt_dt$matching_hicds %in% hicds_of_interest & 
                                   change_cptmt_dt$matching_exprds %in% exprds_of_interest &
                                     change_cptmt_dt$ref_hicds %in% hicds_of_interest & 
                                     change_cptmt_dt$ref_exprds %in% exprds_of_interest,]

change_cptmt_dt$hicds <- change_cptmt_dt$ref_hicds
change_cptmt_dt$exprds <- change_cptmt_dt$ref_exprds
change_cptmt_dt$region <- change_cptmt_dt$refID

out_dt2 <- merge(dt, change_cptmt_dt[,c("hicds", "exprds", "region",
                                 "ref_region_cptmt_eight", "ref_region_pval",
                                 "matching_hicds", "matching_exprds",
                                 "matching_region_cptmt_eight", "matching_region_pval", "rankDiff")],
                by=c("hicds", "exprds", "region"),all.x=T,all.y=F)
# stopifnot(!is.na(out_dt2))
which(apply(out_dt2, 1, function(x)any(is.na(x))))
out_dt2 <- out_dt2[order(out_dt2$hicds, out_dt2$exprds, out_dt2$ref_region_pval),]
rownames(out_dt2) <- NULL

# out_dt$meanLogFC <- round(out_dt$meanLogFC,4)
# out_dt$meanCorr <- round(out_dt$meanCorr,4)
# # out_dt$adjPvalComb <- round(out_dt$adjPvalComb,4)
# # out_dt$ratioDown <- round(out_dt$ratioDown,4)
# out_dt$tadCptmtNormRank <- round(out_dt$tadCptmtNormRank,4)
out_dt$tadCptmtScore <- as.numeric(as.character(out_dt$tadCptmtScore))
# stopifnot(!is.na(out_dt$tadCptmtScore))
# out_dt$tadCptmtScore <- round(out_dt$tadCptmtScore,4)


out_dt$dataset <- file.path(out_dt$hicds, out_dt$exprds)
head_out_dt = do.call(rbind, by(out_dt, out_dt$dataset, function(x) {tmp=x[1:10,];tmp$tad_rank <- rank(tmp$adjPvalComb); tmp}))

out_dt <- out_dt[,keepCols]
colnames(out_dt) <- gsub("tad_", "", colnames(out_dt))
colnames(out_dt) <- gsub("tadCpt", "cpt", colnames(out_dt))

outFile <- file.path(outFolder, paste0("RWPE1_22Rv1_prad_results.Rdata"))
save(out_dt, version=2,file=outFile)
cat(paste0("... written: ", outFile, "\n"))

stopifnot(out_dt$hicds %in% names(hicds_names))
out_dt$hicds <- hicds_names[paste0(out_dt$hicds)]

stopifnot(out_dt$exprds %in% names(exprds_names))
out_dt$exprds <- exprds_names[paste0(out_dt$exprds)]

outFile <- file.path(outFolder, paste0("RWPE1_22Rv1_prad_results.txt"))
write.table(out_dt, file=outFile,sep="\t", quote=F, col.names=T, row.names=F)
cat(paste0("... written: ", outFile, "\n"))

# check
dt1=get(load("PREP_PRAD_DATA/RWPE1_22Rv1_prad_results.Rdata"))
dt2=get(load("PREP_REVISIONFIG1_DATA/revision_fig1_cptmtAnnot_with_corr_purity.Rdata"))
m_dt <- merge(dt1, dt2, by=c("hicds", "exprds", "region"))
m_dt <- na.omit(m_dt)
stopifnot(m_dt$eightCptmtLab.x == m_dt$eightCptmtLab.y)