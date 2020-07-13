# cmp_go_all_partial_filter.R
runFolder <- file.path("../../v2_Yuanlong_Cancer_HiC_data_TAD_DA")

outFolder <- "CMP_GO_ALL_PARTIAL_FILTER"
dir.create(outFolder, recursive = T)

final_dt <- get(load(file.path(runFolder, "CREATE_FINAL_TABLE/all_result_dt.Rdata")))
final_dt$regID <- file.path(final_dt$hicds, final_dt$exprds, final_dt$region)

############################################################################################################## COMPARE THE CONSERVED REGIONS -> RATIO DISC. IN CONSERV. REGIONS

conservRegions_all_dt <- get(load(file.path(runFolder,
                                            "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2", 
                                            "conserved_regions_with_genes_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))
conservRegions_all_dt$nDS <- as.numeric(sapply(conservRegions_all_dt$corresp_tads, function(x) length(unlist(strsplit(x, ",")))))
conservRegions_all_dt <- conservRegions_all_dt[order(conservRegions_all_dt$nDS, decreasing = T),]
head(conservRegions_all_dt[,c("conserved_region", "intersect_genes_symbol", "nDS")],20)

conservRegions_tads_dt <- do.call(rbind, apply(conservRegions_all_dt, 1, function(x) data.frame(conserved_region=unique(x["conserved_region"]), 
                                                       region=unlist(strsplit(x["corresp_tads"], ",")),
                                                        stringsAsFactors = FALSE)))
rownames(conservRegions_tads_dt) <- NULL
all_conservTADs <- unlist(strsplit(conservRegions_all_dt$corresp_tads, split = ","))


discardTADs_aran <- get(load(file.path(runFolder, 
                              "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER", "log10", 
                              "signifTADs_to_discard.Rdata")))


discardTADs_EPIC <- get(load(file.path(runFolder, 
                                       "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER", "EPIC", "log10", 
                                       "signifTADs_to_discard.Rdata")))

conservRegions_tads_dt$discardAran <- conservRegions_tads_dt$region %in% discardTADs_aran
conservRegions_tads_dt$discardEPIC <- conservRegions_tads_dt$region %in% discardTADs_EPIC
  

aggAran_conservRegions_dt <- aggregate(discardAran~conserved_region, FUN=mean, data=conservRegions_tads_dt)
stopifnot(aggAran_conservRegions_dt$discardAran <=1 & aggAran_conservRegions_dt$discardAran >=0)
aggAran_conservRegions_dt <- aggAran_conservRegions_dt[order(aggAran_conservRegions_dt$discardAran, decreasing = TRUE),]
colnames(aggAran_conservRegions_dt)[colnames(aggAran_conservRegions_dt) == "discardAran"] <- "ratioDiscAran"

aggEPIC_conservRegions_dt <- aggregate(discardEPIC~conserved_region, FUN=mean, data=conservRegions_tads_dt)
stopifnot(aggEPIC_conservRegions_dt$discardEPIC <=1 & aggEPIC_conservRegions_dt$discardEPIC >=0)
aggEPIC_conservRegions_dt <- aggEPIC_conservRegions_dt[order(aggEPIC_conservRegions_dt$discardEPIC, decreasing = TRUE),]
colnames(aggEPIC_conservRegions_dt)[colnames(aggEPIC_conservRegions_dt) == "discardEPIC"] <- "ratioDiscEPIC"

out_disc_dt <- merge( merge(conservRegions_all_dt[,c("conserved_region","nDS", "intersect_genes_symbol")],aggAran_conservRegions_dt, by="conserved_region", all=TRUE ),
                      aggEPIC_conservRegions_dt, by="conserved_region", all=TRUE )
stopifnot(!is.na(out_disc_dt))
out_disc_dt <- out_disc_dt[order(out_disc_dt$ratioDiscAran, out_disc_dt$ratioDiscEPIC,decreasing = TRUE),]
head(out_disc_dt, 10)
out_disc_dt <- out_disc_dt[order(out_disc_dt$nDS, out_disc_dt$ratioDiscAran, out_disc_dt$ratioDiscEPIC,decreasing = TRUE),]
head(out_disc_dt, 10)



############################################################################################################## COMPARE THE CONSERVED REGIONS -> CORRESP. TADs - epic
final_dt <- get(load(file.path(runFolder, "CREATE_FINAL_TABLE/all_result_dt.Rdata")))
final_dt$regID <- file.path(final_dt$hicds, final_dt$exprds, final_dt$region)
final_dt$region <- final_dt$regID

conservRegions_all_dt <- get(load(file.path(runFolder,
                                            "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2", 
                                            "conserved_regions_with_genes_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))
conservRegions_all_dt$nDS <- as.numeric(sapply(conservRegions_all_dt$corresp_tads, function(x) length(unlist(strsplit(x, ",")))))
conservRegions_all_dt <- conservRegions_all_dt[order(conservRegions_all_dt$nDS, decreasing = T),]
head(conservRegions_all_dt[,c("conserved_region", "intersect_genes_symbol", "nDS")],20)
conservRegions_tads_dt <- do.call(rbind, apply(conservRegions_all_dt, 1, function(x) data.frame(conserved_region=unique(x["conserved_region"]), 
                                                                                                region=unlist(strsplit(x["corresp_tads"], ",")),
                                                                                                stringsAsFactors = FALSE)))
rownames(conservRegions_tads_dt) <- NULL

conserv_thresh <- 8
conservedRegions_aboveThresh <- as.character(conservRegions_all_dt$conserved_region[conservRegions_all_dt$nDS >= conserv_thresh])
all_conservedTADs_aboveThresh <- conservRegions_tads_dt$region[conservRegions_all_dt$conserved_region %in% conservedRegions_aboveThresh]

conservRegions_filterEPIC_dt <- get(load(file.path(runFolder,
                                                   "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER",
                                                   "EPIC", "log10", 
                                                   "conserved_regions_with_genes_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))
conservRegions_filterEPIC_dt$nDS <- as.numeric(sapply(conservRegions_filterEPIC_dt$corresp_tads, function(x) length(unlist(strsplit(x, ",")))))
conservRegions_filterEPIC_dt <- conservRegions_filterEPIC_dt[order(conservRegions_filterEPIC_dt$nDS, decreasing = T),]
head(conservRegions_filterEPIC_dt[,c("conserved_region", "intersect_genes_symbol", "nDS")],20)
conservRegions_filterEPIC_tads_dt <- do.call(rbind, apply(conservRegions_filterEPIC_dt, 1, function(x) data.frame(conserved_region=unique(x["conserved_region"]), 
                                                                                                                  region=unlist(strsplit(x["corresp_tads"], ",")),
                                                                                                                  stringsAsFactors = FALSE)))
missedEPIC_all_conservedTADs_aboveThresh <- all_conservedTADs_aboveThresh[!all_conservedTADs_aboveThresh %in% conservRegions_filterEPIC_tads_dt$region ]
cat(paste0(length(missedEPIC_all_conservedTADs_aboveThresh) , "/", length(all_conservedTADs_aboveThresh), "\n"))

out_dt <- merge(conservRegions_tads_dt, conservRegions_all_dt[,c("conserved_region", "intersect_genes_symbol", "nDS")], all=T, by="conserved_region")
stopifnot(!is.na(out_dt))
out_dt <- merge(out_dt, final_dt[,c("region", "region_genes")], by="region", all.x=T, all.y=F)
stopifnot(!is.na(out_dt))
missed_out_dt <- out_dt[out_dt$region %in% missedEPIC_all_conservedTADs_aboveThresh & out_dt$conserved_region %in% conservedRegions_aboveThresh,]
stopifnot(missed_out_dt$conserved_region %in% conservedRegions_aboveThresh)
stopifnot(missed_out_dt$region %in% all_conservedTADs_aboveThresh)
stopifnot(missed_out_dt$nDS >= conserv_thresh)
aboveThresh_missedEPIC_dt <- missed_out_dt[order(missed_out_dt$nDS, decreasing = T),]
rownames(aboveThresh_missedEPIC_dt) <- NULL
head(aboveThresh_missedEPIC_dt[,c("region", "conserved_region", "region_genes")], 20)
unique(aboveThresh_missedEPIC_dt[,c("conserved_region", "intersect_genes_symbol", "nDS")])


############################################################################################################## COMPARE THE CONSERVED REGIONS -> CORRESP. TADs - aran
final_dt <- get(load(file.path(runFolder, "CREATE_FINAL_TABLE/all_result_dt.Rdata")))
final_dt$regID <- file.path(final_dt$hicds, final_dt$exprds, final_dt$region)
final_dt$region <- final_dt$regID

conservRegions_all_dt <- get(load(file.path(runFolder,
                                            "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2", 
                                            "conserved_regions_with_genes_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))
conservRegions_all_dt$nDS <- as.numeric(sapply(conservRegions_all_dt$corresp_tads, function(x) length(unlist(strsplit(x, ",")))))
conservRegions_all_dt <- conservRegions_all_dt[order(conservRegions_all_dt$nDS, decreasing = T),]
head(conservRegions_all_dt[,c("conserved_region", "intersect_genes_symbol", "nDS")],20)
conservRegions_tads_dt <- do.call(rbind, apply(conservRegions_all_dt, 1, function(x) data.frame(conserved_region=unique(x["conserved_region"]), 
                                                                                                region=unlist(strsplit(x["corresp_tads"], ",")),
                                                                                                stringsAsFactors = FALSE)))
rownames(conservRegions_tads_dt) <- NULL

conserv_thresh <- 8
conservedRegions_aboveThresh <- as.character(conservRegions_all_dt$conserved_region[conservRegions_all_dt$nDS >= conserv_thresh])
all_conservedTADs_aboveThresh <- conservRegions_tads_dt$region[conservRegions_all_dt$conserved_region %in% conservedRegions_aboveThresh]


conservRegions_filterAran_dt <- get(load(file.path(runFolder,
                                                    "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER", "log10", 
                                                    "conserved_regions_with_genes_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))
conservRegions_filterAran_dt$nDS <- as.numeric(sapply(conservRegions_filterAran_dt$corresp_tads, function(x) length(unlist(strsplit(x, ",")))))
conservRegions_filterAran_dt <- conservRegions_filterAran_dt[order(conservRegions_filterAran_dt$nDS, decreasing = T),]
head(conservRegions_filterAran_dt[,c("conserved_region", "intersect_genes_symbol", "nDS")],20)
conservRegions_filterAran_tads_dt <- do.call(rbind, apply(conservRegions_filterAran_dt, 1, function(x) data.frame(conserved_region=unique(x["conserved_region"]), 
                                                                                                                  region=unlist(strsplit(x["corresp_tads"], ",")),
                                                                                                                  stringsAsFactors = FALSE)))
missedAran_all_conservedTADs_aboveThresh <- all_conservedTADs_aboveThresh[!all_conservedTADs_aboveThresh %in% conservRegions_filterAran_tads_dt$region]
cat(paste0(length(missedAran_all_conservedTADs_aboveThresh) , "/", length(all_conservedTADs_aboveThresh), "\n"))

out_dt <- merge(conservRegions_tads_dt, conservRegions_all_dt[,c("conserved_region", "intersect_genes_symbol", "nDS")], all=T, by="conserved_region")
stopifnot(!is.na(out_dt))
out_dt <- merge(out_dt, final_dt[,c("region", "region_genes")], by="region", all.x=T, all.y=F)
stopifnot(!is.na(out_dt))
missed_out_dt <- out_dt[out_dt$region %in% missedAran_all_conservedTADs_aboveThresh & out_dt$conserved_region %in% conservedRegions_aboveThresh,]
stopifnot(missed_out_dt$conserved_region %in% conservedRegions_aboveThresh)
stopifnot(missed_out_dt$region %in% all_conservedTADs_aboveThresh)
stopifnot(missed_out_dt$nDS >= conserv_thresh)
aboveThresh_missedAran_dt <- missed_out_dt[order(missed_out_dt$nDS, decreasing = T),]
rownames(aboveThresh_missedAran_dt) <- NULL
head(aboveThresh_missedAran_dt, 20)
unique(aboveThresh_missedEPIC_dt[,c("conserved_region", "intersect_genes_symbol", "nDS")])





############################################################################################################## 
############################################################################################################## genes behind GOs
############################################################################################################## 

setDir <- "/media/electron"
# setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
entrez2symb <- setNames(gff_dt$symbol,gff_dt$entrezID)

signifGO_thresh <- 0.05

go_dt <- get(load("GO_SIGNIF_ACROSS_HICDS_v2/conserved_signif_enrich_resultDT.Rdata"))
stopifnot(go_dt$log10_pval == -log10(go_dt$p.adjust))
stopifnot(rownames(go_dt) == go_dt$ID)
rownames(go_dt) <- NULL
signif_go_dt <- go_dt[go_dt$p.adjust <= signifGO_thresh,c("ID", "geneID", "p.adjust")]

signif_go_dt$symbols <- sapply(signif_go_dt$geneID, function(x){
 all_entrez <- unlist(strsplit(x, split="/")) 
 all_symbols <- entrez2symb[paste0(all_entrez)]
 stopifnot(!is.na(all_symbols))
 paste0(all_symbols, collapse=",")
})
head(signif_go_dt[,c("ID", "p.adjust", "symbols")],10)

outFile <- file.path(outFolder, "signif_go.txt")
write.table(signif_go_dt, file=outFile, sep="\t", quote=F, col.names=TRUE,row.names=F, append=F)

partial_go_dt <- get(load("GO_SIGNIF_ACROSS_HICDS_v2_PARTIAL//conserved_signif_enrich_resultDT.Rdata"))
stopifnot(partial_go_dt$log10_pval == -log10(partial_go_dt$p.adjust))
stopifnot(rownames(partial_go_dt) == partial_go_dt$ID)
rownames(partial_go_dt) <- NULL
signif_partial_go_dt <- partial_go_dt[partial_go_dt$p.adjust <= signifGO_thresh,c("ID", "geneID", "p.adjust")]

signif_partial_go_dt$symbols <- sapply(signif_partial_go_dt$geneID, function(x){
  all_entrez <- unlist(strsplit(x, split="/")) 
  all_symbols <- entrez2symb[paste0(all_entrez)]
  stopifnot(!is.na(all_symbols))
  paste0(all_symbols, collapse=",")
})
head(signif_partial_go_dt[,c("ID", "p.adjust", "symbols")],10)
outFile <- file.path(outFolder, "signif_partial_go.txt")
write.table(signif_partial_go_dt, file=outFile, sep="\t", quote=F, col.names=TRUE,row.names=F, append=F)

aran_go_dt <- get(load("GO_SIGNIF_ACROSS_HICDS_v2_PURITYFILTER/log10//conserved_signif_enrich_resultDT.Rdata"))
stopifnot(aran_go_dt$log10_pval == -log10(aran_go_dt$p.adjust))
stopifnot(rownames(aran_go_dt) == aran_go_dt$ID)
rownames(aran_go_dt) <- NULL
signif_aran_go_dt <- aran_go_dt[aran_go_dt$p.adjust <= signifGO_thresh,c("ID", "geneID", "p.adjust")]

signif_aran_go_dt$symbols <- sapply(signif_aran_go_dt$geneID, function(x){
  all_entrez <- unlist(strsplit(x, split="/")) 
  all_symbols <- entrez2symb[paste0(all_entrez)]
  stopifnot(!is.na(all_symbols))
  paste0(all_symbols, collapse=",")
})
head(signif_aran_go_dt[,c("ID", "p.adjust", "symbols")],10)
outFile <- file.path(outFolder, "signif_aran_go.txt")
write.table(signif_aran_go_dt, file=outFile, sep="\t", quote=F, col.names=TRUE,row.names=F, append=F)


epic_go_dt <- get(load("GO_SIGNIF_ACROSS_HICDS_v2_PURITYFILTER/EPIC/log10//conserved_signif_enrich_resultDT.Rdata"))
stopifnot(epic_go_dt$log10_pval == -log10(epic_go_dt$p.adjust))
stopifnot(rownames(epic_go_dt) == epic_go_dt$ID)
rownames(epic_go_dt) <- NULL
signif_epic_go_dt <- epic_go_dt[epic_go_dt$p.adjust <= signifGO_thresh,c("ID", "geneID", "p.adjust")]

signif_epic_go_dt$symbols <- sapply(signif_epic_go_dt$geneID, function(x){
  all_entrez <- unlist(strsplit(x, split="/")) 
  all_symbols <- entrez2symb[paste0(all_entrez)]
  stopifnot(!is.na(all_symbols))
  paste0(all_symbols, collapse=",")
})
head(signif_epic_go_dt[,c("ID", "p.adjust", "symbols")],10)

outFile <- file.path(outFolder, "signif_epic_go.txt")
write.table(signif_epic_go_dt, file=outFile, sep="\t", quote=F, col.names=TRUE,row.names=F, append=F)






# 
# 
# 
# missed_Aran_dt <- conservRegions_tads_dt[conservRegions_tads_dt$region %in% missedAran_all_conservedTADs_aboveThresh,]
# missed_Aran_dt <- merge(missed_Aran_dt, conservRegions_all_dt, by="conserved_region", all.x=T, all.y=F)
# missed_Aran_dt <- missed_Aran_dt[missed_Aran_dt$nDS >= conserv_thresh,]
# stopifnot(missedAran_all_conservedTADs_aboveThresh %in% missed_Aran_dt$region)
# 
# 
# 
# 
# conservRegions_filterAran_tads_dt <- do.call(rbind, apply(conservRegions_filterAran_dt, 1, function(x) data.frame(conserved_region=unique(x["conserved_region"]), 
#                                                                                                 region=unlist(strsplit(x["corresp_tads"], ",")),
#                                                                                                 stringsAsFactors = FALSE)))
# aran_conservedRegions_aboveThresh <- conservRegions_filterAran_dt$conserved_region[conservRegions_filterAran_dt$nDS >= conserv_thresh]
# aran_conservedTADs_aboveThresh <- conservRegions_filterAran_tads_dt$region[conservRegions_filterAran_dt$conserved_region %in% aran_conservedRegions_aboveThresh]
# 
# 
# #### regions missed by Aran filter
# missed_regs <- setdiff(all_conservedTADs_aboveThresh, aran_conservedTADs_aboveThresh)
# stopifnot(missed_regs %in% all_conservedTADs_aboveThresh)
# stopifnot(!missed_regs %in% aran_conservedTADs_aboveThresh)
# 
# all_missedByAran <- final_dt[final_dt$regID %in% missed_regs, c("regID", "region_genes")]
# colnames(all_missedByAran)[1] <- "region"
# all_missedByAran <- merge(all_missedByAran, conservRegions_tads_dt[,c("region", "conserved_region")], all.x=T, all.y=F, by="region")
# stopifnot(all_missedByAran$conserved_region %in% conservedRegions_aboveThresh)
# 
# 
# all_missedByAran <- all_missedByAran[all_missedByAran$conserved_region %in% conservedRegions_aboveThresh,]
# stopifnot(setequal(all_missedByAran$region, missed_regs))
# 
# 
# all_missedByAran <- merge(all_missedByAran, conservRegions_all_dt[,c("conserved_region", "nDS")], all.x=T, all.y=F, by="conserved_region")
# stopifnot(setequal(all_missedByAran$region, setdiff(all_conservedTADs_aboveThresh, aran_conservedTADs_aboveThresh)))
# all_missedByAran <- all_missedByAran[all_missedByAran$nDS >= conserv_thresh,]
# stopifnot(all_missedByAran$region %in% all_conservedTADs_aboveThresh)
# stopifnot(all_missedByAran$conserved_region %in% conservedRegions_aboveThresh)
# stopifnot(setequal(all_missedByAran$region, setdiff(all_conservedTADs_aboveThresh, aran_conservedTADs_aboveThresh)))
# 
# 
# 
# 
# 
# 
# 
# conservRegions_filterEPIC_dt <- get(load(file.path(runFolder,
#                                                    "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER", "EPIC", "log10", 
#                                                    "conserved_regions_with_genes_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))
# 
# conservRegions_filterEPIC_dt$nDS <- as.numeric(sapply(conservRegions_filterEPIC_dt$corresp_tads, function(x) length(unlist(strsplit(x, ",")))))
# conservRegions_filterEPIC_dt <- conservRegions_filterEPIC_dt[order(conservRegions_filterEPIC_dt$nDS, decreasing = T),]
# head(conservRegions_filterEPIC_dt[,c("conserved_region", "intersect_genes_symbol", "nDS")],20)
# 
# 
# conservRegions_filterEPIC_tads_dt <- do.call(rbind, apply(conservRegions_filterEPIC_dt, 1, function(x) data.frame(conserved_region=unique(x["conserved_region"]), 
#                                                                                                                   region=unlist(strsplit(x["corresp_tads"], ",")),
#                                                                                                                   stringsAsFactors = FALSE)))
# EPIC_conservedRegions_aboveThresh <- conservRegions_filterEPIC_dt$conserved_region[conservRegions_filterEPIC_dt$nDS >= conserv_thresh]
# EPIC_conservedTADs_aboveThresh <- conservRegions_filterEPIC_tads_dt$region[conservRegions_filterEPIC_dt$conserved_region %in% EPIC_conservedRegions_aboveThresh]
# 
# 
# 
# 
# 
# 
# 
# 
# 
#   
# keepTADs_aran <- get(load(file.path(runFolder, 
#                               "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER", "log10", 
#                               "signifTADs_to_keep.Rdata")))
# 
# 
# 
# 
# keepTADs_EPIC <- get(load(file.path(runFolder, 
#                               "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER", "EPIC", "log10", 
#                               "signifTADs_to_keep.Rdata")))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# conservRegions_filter_aran_dt <- get(load(file.path(runFolder,
#                                                     "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER", "log10", 
#                                                     "conserved_regions_with_genes_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))
# filterAran_conservTADs <- unlist(strsplit(conservRegions_filter_aran_dt$corresp_tads, split = ","))
# 
# conservRegions_filter_epic_dt <- get(load(file.path(runFolder,
#                                                     "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER", "EPIC", "log10", 
#                                                     "conserved_regions_with_genes_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))
# filterEPIC_conservTADs <- unlist(strsplit(conservRegions_filter_epic_dt$corresp_tads, split = ","))
# 
# conservRegions_partial_dt <- get(load(file.path(runFolder,
#                                                 "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PARTIAL", 
#                                                 "conserved_signif_tadsadjPvalComb0.01_minBpRatio0.8_minInterGenes3.Rdata")))
# partial_conservTADs <- unlist(strsplit(conservRegions_partial_dt$corresp_tads, split = ","))
# 
# 
# # retrieve TADs that are missed
# notAran_tads <- all_conservTADs[all_conservTADs %in% filterAran_conservTADs]
# 
# 
# 
# 
# 
# go_all_dt <- get(load("GO_SIGNIF_ACROSS_HICDS_v2/conserved_signif_enrich_resultDT.Rdata"))
# go_filter_epic_dt <- get(load("GO_SIGNIF_ACROSS_HICDS_v2_PURITYFILTER//conserved_signif_enrich_resultDT.Rdata"))
# go_filter_aran_dt <- get(load("GO_SIGNIF_ACROSS_HICDS_v2_PURITYFILTER//conserved_signif_enrich_resultDT.Rdata"))
# go_partial_dt <- get(load("GO_SIGNIF_ACROSS_HICDS_v2_PARTIAL/conserved_signif_enrich_resultDT.Rdata"))