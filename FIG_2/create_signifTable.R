# Rscript create_signifTable.R

plotType <- "svg"
source("../settings.R")

outFolder <- "CREATE_SIGNIF_TABLE"
dir.create(outFolder, recursive = TRUE)

geneSignifThresh <- 0.01
tadSignifThresh <- 0.01

final_dt <- get(load(file.path(runFolder, "CREATE_FINAL_TABLE/all_result_dt.Rdata")))

agg_final_dt <- aggregate(adjPvalComb~hicds+exprds, function(x)sum(x<=tadSignifThresh), data = final_dt)
agg_final_dt <- agg_final_dt[order(agg_final_dt$hicds, agg_final_dt$exprds),]

signif_final_dt <- final_dt[final_dt$adjPvalComb <= tadSignifThresh, c("hicds", "exprds", "region", "adjPvalComb")]
signif_final_dt <- signif_final_dt[order(signif_final_dt$hicds, signif_final_dt$exprds, signif_final_dt$region), ]


all_ranks_dt <- get(load(file.path(runFolder, "GENE_RANK_TAD_RANK/all_gene_tad_signif_dt.Rdata")))

agg_all_dt <- aggregate(adj.P.Val~hicds+exprds, data=all_ranks_dt, FUN=function(x)sum(x<=geneSignifThresh))
agg_all_dt <- agg_all_dt[order(agg_all_dt$hicds, agg_all_dt$exprds),]


signif_tad_dt <- all_ranks_dt[all_ranks_dt$tad_adjCombPval <= tadSignifThresh,c("hicds", "exprds", "region", "tad_adjCombPval")]
signif_tad_dt <- signif_tad_dt[order(signif_tad_dt$hicds, signif_tad_dt$exprds, signif_tad_dt$region), ]
signif_tad_dt <- unique(signif_tad_dt)

stopifnot(nrow(signif_tad_dt) == nrow(signif_final_dt))
stopifnot(all.equal(signif_tad_dt, signif_final_dt, check.attributes=F))

signif_gene_dt <- all_ranks_dt[all_ranks_dt$adj.P.Val <= geneSignifThresh,c("hicds", "exprds", "entrezID", "adj.P.Val")]
signif_gene_dt <- signif_gene_dt[order(signif_gene_dt$hicds, signif_gene_dt$exprds, signif_gene_dt$entrezID), ]
signif_gene_dt <- unique(signif_gene_dt)

nSignif_tad_dt <- aggregate(region~hicds+exprds, FUN=length, data=signif_tad_dt)
nSignif_tad_dt <- nSignif_tad_dt[order(nSignif_tad_dt$hicds, nSignif_tad_dt$exprds),]

stopifnot(all.equal(nSignif_tad_dt, agg_final_dt, check.attributes=F))

nSignif_gene_dt <- aggregate(entrezID~hicds+exprds, FUN=length, data=signif_gene_dt)
nSignif_gene_dt <- nSignif_gene_dt[order(nSignif_gene_dt$hicds, nSignif_gene_dt$exprds),]

stopifnot(all.equal(nSignif_gene_dt, agg_all_dt, check.attributes=F))

colnames(nSignif_gene_dt)[3] <- paste0("nSignifGenes_", geneSignifThresh)
colnames(nSignif_tad_dt)[3] <- paste0("nSignifTADs_", tadSignifThresh)

final_dt <- merge(nSignif_gene_dt, nSignif_tad_dt, by=c("hicds", "exprds"), all.x=T, all.y=T)
stopifnot(!is.na(final_dt))
final_dt <- final_dt[order(final_dt$hicds, final_dt$exprds),]

outFile <- file.path(outFolder, paste0("nbrGeneSignif", geneSignifThresh, "_nbrTadSignif", tadSignifThresh,"_byDS.txt"))
write.table(final_dt, file=outFile, sep="\t", quote=F, col.names=T, row.names=F)
cat(paste0("... written:" , outFile, "\n"))



