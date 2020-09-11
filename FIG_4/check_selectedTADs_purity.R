purity_dt <- get(load("../../v2_Yuanlong_Cancer_HiC_data_TAD_DA/ALLTADS_AND_PURITY_FINAL/aran/CPE/log10/all_ds_corrPurity_dt.Rdata"))
purity_dt$region_id <- file.path(purity_dt$dataset, purity_dt$region)


 selectedTADs <- c("ENCSR489OCU_NCI-H460_40kb/TCGAluad_mutKRAS_mutEGFR/chr10_TAD16",
 "ENCSR489OCU_NCI-H460_40kb/TCGAluad_mutKRAS_mutEGFR/chr17_TAD162"      ,
"ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/chr10_TAD268",
"ENCSR489OCU_NCI-H460_40kb/TCGAluad_norm_luad/chr11_TAD390" )
 
 out_dt <- purity_dt[purity_dt$region_id %in% selectedTADs,c("region_id", "purityCorr")]
 rownames(out_dt) <- NULL
 agg_out_dt <- aggregate(purityCorr~region_id, FUN=mean, data=out_dt)
 colnames(agg_out_dt)[2] <- "meanCorr"
 agg_out_dt[order(agg_out_dt$meanCorr, decreasing=F),]

 