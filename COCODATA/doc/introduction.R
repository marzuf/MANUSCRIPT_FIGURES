## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------

rm(list=ls())

if(!require(COCODATA))
  devtools::install_github("marzuf/MANUSCRIPT_FIGURES", subdir="COCODATA")
  # alternatively: 
  # install.packages("COCODATA_0.0.0.1.tar.gz", repos = NULL, type ="source")
 # data("norm_ID")
library(COCODATA)


nCpu <- 2
library(doMC)
library(foreach)
registerDoMC(nCpu)

## ----load_data-----------------------------------------------------------
# ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_DE_topTable.Rdata
# ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_meanLogFC_permDT.Rdata
# ENCSR489OCU_NCI-H460_all_assigned_regions.txt
# ENCSR489OCU_NCI-H460_all_genes_positions.txt
# luad_ID.Rdata
# norm_ID.Rdata
# TCGAluad_norm_luad_fpkmDT.Rdata
# TCGAluad_norm_luad_rnaseqDT_v2.Rdata



## ----TO_FUNCTION_TO_DOCUMENT---------------------------------------------
# auc.Rd
# get_downByRegion_v2.Rd
# get_fcc.Rd
# get_FCdownByRegion_v2.Rd
# get_multiShuffledPositions_vFunct.Rd
# get_ratioDown.Rd
# get_ratioFC.Rd
# get_ShuffledPositions_vFunct.Rd
# get_statFromShuffle_para.Rd
# madNorm.Rd # not used, this was for microarray data
# printAndLog.Rd
# quantNorm.Rd
# stouffer.Rd
# get_auc_ratio

## ----prep_obs------------------------------------------------------------
# table from gene-level DE analysis:
data("ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_DE_topTable") # this loads DE_topTable
DE_topTable$genes <- as.character(DE_topTable$genes)

# list of genes used in the pipeline
data("ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_pipeline_geneList") # this loads pipeline_geneList.Rdata

# for those genes, I have logFC data
stopifnot(names(pipeline_geneList) %in% DE_topTable$genes)
DE_topTable <- DE_topTable[DE_topTable$genes %in% names(pipeline_geneList),]
DE_topTable$entrezID <- pipeline_geneList[DE_topTable$genes]
stopifnot(!is.na(DE_topTable$entrezID))


# table with gene-to-TAD assignment
gene2tad_dt <- read.delim(system.file("extdata", "ENCSR489OCU_NCI-H460_all_genes_positions.txt", package = "COCODATA"),
                          stringsAsFactors = FALSE,
                          header=FALSE, 
                          col.names=c("entrezID", "chromosome", "start", "end", "region"))
gene2tad_dt$entrezID <- as.character(gene2tad_dt$entrezID)


# take only the genes used in the pipeline
pip_g2t_dt <- gene2tad_dt[gene2tad_dt$entrezID %in% pipeline_geneList,]

# merge to match gene-to-TAD and logFC
merged_dt <- merge(pip_g2t_dt[,c("entrezID", "region")], DE_topTable[,c("logFC", "entrezID")], by="entrezID", all.x=TRUE, all.y=FALSE)
stopifnot(!is.na(merged_dt))

## ----ratioDown_obs_example-----------------------------------------------
# compute the FCC for each TAD
FCC_dt <- aggregate(logFC ~ region, FUN=get_ratioDown, data=merged_dt)
colnames(FCC_dt)[colnames(FCC_dt) == "logFC"] <- "ratioDown"
obs_ratioDown <- setNames(FCC_dt$ratioDown, FCC_dt$region)
save(obs_ratioDown, file="package_obs_ratioDown.Rdata")

## ----ratioFC_obs_example-------------------------------------------------
# compute the FCC for each TAD
FCC_dt <- aggregate(logFC ~ region, FUN=get_ratioFC, data=merged_dt)
colnames(FCC_dt)[colnames(FCC_dt) == "logFC"] <- "ratioFC"
obs_ratioNegFC <- setNames(FCC_dt$ratioFC, FCC_dt$region)
save(obs_ratioNegFC, file="package_obs_ratioNegFC.Rdata")

## ----FCC_obs_example-----------------------------------------------------
# compute the FCC for each TAD
FCC_dt <- aggregate(logFC ~ region, FUN=get_fcc, data=merged_dt)
colnames(FCC_dt)[colnames(FCC_dt) == "logFC"] <- "FCC"
obs_FCC <- setNames(FCC_dt$FCC, FCC_dt$region)
save(obs_FCC, file="package_obs_FCC.Rdata")

## ----FCC_permut_example--------------------------------------------------

# compute the FCC for the permutation data
data("cut1000_ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_permutationsDT") # this loads permutationsDT
head_sq(permutationsDT)

stopifnot(setequal(rownames(permutationsDT), pipeline_geneList))

tad_levels <- as.character(FCC_dt$region)

all_permut_FCC_dt <- foreach(i = 1:ncol(permutationsDT), .combine='cbind') %dopar% {
  perm_g2t <- data.frame(entrezID=as.character(rownames(permutationsDT)),
                        region = as.character(permutationsDT[,i]),
                        stringsAsFactors = FALSE)
  perm_merged_dt <- merge(perm_g2t, DE_topTable[,c("logFC", "entrezID")], by="entrezID", all.x=TRUE, all.y=FALSE)
  stopifnot(!is.na(perm_merged_dt))
  
  # compute the FCC for each TAD
  permut_FCC_dt <- aggregate(logFC ~ region, FUN=get_fcc, data=perm_merged_dt)
  colnames(permut_FCC_dt)[colnames(permut_FCC_dt) == "logFC"] <- "FCC"
  rownames(permut_FCC_dt) <- as.character(permut_FCC_dt$region)
  stopifnot(setequal(rownames(permut_FCC_dt), tad_levels))
  permut_FCC_dt[tad_levels, "FCC", drop=FALSE]
}
head_sq(all_permut_FCC_dt)

## ----FCC_AUC_ratio_example-----------------------------------------------
stopifnot(setequal(names(obs_FCC), rownames(all_permut_FCC_dt)))
get_auc_ratio(fcc_vect=obs_FCC, 
              fcc_permDT=all_permut_FCC_dt, 
              doPlot=T, 
              main="FCC AUC cumsum curves")
mtext(text = paste0("ENCSR489OCU_NCI-H460_40kb - TCGAluad_norm_luad"), side=3)

