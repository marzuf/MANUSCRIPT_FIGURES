options(scipen=100)

SSHFS=F

buildData <- F
buildTable <- F
plotOnly <- T

# Rscript go_signif_across_hicds_v2_withThresh_allBars.R

script_name <- "go_signif_across_hicds_v2_withThresh_allBars.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")
require(clusterProfiler)
require(foreach)
require(doMC)
require(GenomicRanges)
require(ggplot2)
require(RColorBrewer)
require(reshape2)
require(ggpubr)
# require(gplots)
registerDoMC(4)

#source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
#source("my_heatmap.2.R")


args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) {
  data_cmpType <- ""
  file_prefix <- ""
} else if(length(args) == 1) {
  data_cmpType <- args[1]  
  stopifnot(data_cmpType %in% c("norm_vs_tumor", "subtypes", "wt_vs_mut"))
  file_prefix <- paste0(data_cmpType, "_")
} else {
  stop("error\n") 
}



plotType <- "svg"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
axisCex <- 1.4

myHeightHeat <- myHeight * 1.8
myWidthHeat <- myWidth * 1.8

myWidthGG <- 12
myHeightGG <- 12


#source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")
source(file.path(runFolder, "enricher_settings.R"))
# from enricher_settings, load:
# enricher_ontologyType <- "BP"
# enricher_pvalueCutoff <- 1
# enricher_pAdjustMethod <- "BH"
# enricher_minGSSize <- 1
# enricher_maxGSSize <- 500
# enricher_qvalueCutoff <- 1
# enricher_results_sortGOby <- "p.adjust"

conservThresh <- 8


script0_name <- "0_prepGeneData"

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)

mainFolder <- file.path(runFolder)
stopifnot(dir.exists(mainFolder))
pipFolder <- file.path(mainFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipFolder))
all_hicds <- list.files(pipFolder)

all_hicds <- all_hicds[!grepl("RANDOM", all_hicds) & !grepl("PERMUT", all_hicds)]

file.path(mainFolder, all_hicds)[!dir.exists(file.path(mainFolder, all_hicds))]
stopifnot(dir.exists(file.path(mainFolder, all_hicds)))

all_exprds <- lapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))
names(all_exprds) <- all_hicds

outFolder <- file.path("GO_SIGNIF_ACROSS_HICDS_v2_WITHTHRESH_ALLBARS", conservThresh, data_cmpType)
dir.create(outFolder, recursive = TRUE)

all_datasets <- unlist(lapply(1:length(all_exprds), function(x) file.path(names(all_exprds)[x], all_exprds[[x]])))
nDS <- length(all_datasets)

cat(paste0("n allDS = ", length(all_datasets), "\n"))

# => best TAD matching
# in # of genes
# in bp

final_dt_file <- file.path(runFolder, "CREATE_FINAL_TABLE", "all_result_dt.Rdata")
stopifnot(file.exists(final_dt_file))
final_dt <- get(load(final_dt_file))

signif_column <- "adjPvalComb"
signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", signifThresh)

final_dt[, paste0(signifcol)] <- final_dt[, paste0(signif_column)] <= signifThresh

minOverlapBpRatio <- 0.8
minIntersectGenes <- 3

nRegionLolli <- 10

strwdth <- 30


cat(paste0("> signif_column\t=\t", signif_column, "\n"))
cat(paste0("> signifThresh\t=\t", signifThresh, "\n"))
cat(paste0("> minOverlapBpRatio\t=\t", minOverlapBpRatio, "\n"))
cat(paste0("> minIntersectGenes\t=\t", minIntersectGenes, "\n"))
cat(paste0("> nRegionLolli\t=\t", nRegionLolli, "\n"))


inFolder <- file.path(runFolder, "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2", data_cmpType)
stopifnot(dir.exists(inFolder))

inFile <- file.path(inFolder, paste0(file_prefix, "conserved_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
stopifnot(file.exists(inFile))
conserved_signif_tads <- get(load(inFile))

## update here to filter: take only regions with conservThresh
conserved_signif_tads <- conserved_signif_tads[lengths(conserved_signif_tads) >= conservThresh]

stopifnot(lengths(conserved_signif_tads) >= conservThresh)

kept_conserved_signif_tads <- conserved_signif_tads
outFile <- file.path(outFolder, paste0(file_prefix, "kept_conserved_signif_tads.Rdata"))
save(kept_conserved_signif_tads, file = outFile, version=2)

#checked: there is no duplicate in the conserved datasets
ds_signif_tads <- lapply(conserved_signif_tads, function(x) unique(dirname(x)))
stopifnot(lengths(conserved_signif_tads) == lengths(ds_signif_tads))


inFile <- file.path(inFolder, paste0(file_prefix, "signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
stopifnot(file.exists(inFile))
signif_tads <- get(load(inFile))

all_conserved_signif_tads <- unique(unlist(conserved_signif_tads))
stopifnot(all_conserved_signif_tads %in% signif_tads)

all_not_conserved_signif_tads <- signif_tads[! signif_tads %in% all_conserved_signif_tads]

####################################################################################################################################### >>> prepare the data
if(buildTable) {
  cat("... start preparing data before matching \n")
  
  hicds = all_hicds[1]
  all_gene_list <- foreach(hicds = all_hicds) %dopar% {
    
    hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_assigned_regions.txt")
    stopifnot(file.exists(hicds_file))
    tadpos_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))
    tadpos_dt <- tadpos_dt[grepl("_TAD", tadpos_dt$region),]
    stopifnot(nrow(tadpos_dt) > 0 )
    
    hicds_file <- file.path(mainFolder, hicds, "genes2tad", "all_genes_positions.txt")
    stopifnot(file.exists(hicds_file))
    g2t_dt <- read.delim(hicds_file, header=F, stringsAsFactors = FALSE, col.names=c("entrezID", "chromo", "start", "end", "region"))
    g2t_dt$entrezID <- as.character(g2t_dt$entrezID)
    stopifnot(nrow(g2t_dt) > 0 )
    
    stopifnot(g2t_dt$entrezID %in% gff_dt$entrezID)
    
    exprds = all_exprds[[paste0(hicds)]][1]
    exprds_list <- foreach(exprds = all_exprds[[paste0(hicds)]]) %do% {
      
      
      signif_tads <- final_dt$region[final_dt$hicds == hicds & final_dt$exprds == exprds & final_dt[, paste0(signifcol)] ]
      
      gene_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_geneList.Rdata")
      stopifnot(file.exists(gene_file))  
      
      region_file <- file.path(pipFolder, hicds, exprds, script0_name, "pipeline_regionList.Rdata")
      stopifnot(file.exists(region_file))  
      
      geneList <- get(load(gene_file))
      regionList <- get(load(region_file))
      
      stopifnot(geneList %in% g2t_dt$entrezID)
      ref_g2t_dt <- g2t_dt[g2t_dt$entrezID %in% geneList,]
      stopifnot(setequal(ref_g2t_dt$region, regionList))
      
      stopifnot(regionList %in% tadpos_dt$region)
      ref_tadpos_dt <- tadpos_dt[tadpos_dt$region %in% regionList,]
      stopifnot(setequal(ref_g2t_dt$region, ref_tadpos_dt$region))
      
      # take only the data from signif TADs
      stopifnot(setequal(ref_g2t_dt$entrezID, geneList))
      
      
      all_gene_tads <- setNames(file.path(hicds, exprds, ref_g2t_dt$region), ref_g2t_dt$entrezID)

      conserved_signif_tads_genes <- names(all_gene_tads[all_gene_tads %in% all_conserved_signif_tads]) 
      
      not_conserved_signif_tads_genes <- names(all_gene_tads[all_gene_tads %in% all_not_conserved_signif_tads]) 
      
      stopifnot(length(intersect(conserved_signif_tads_genes, not_conserved_signif_tads_genes)) == 0)
      
      not_signif_tads_genes <- geneList[!geneList %in% c(conserved_signif_tads_genes, not_conserved_signif_tads_genes) ]

      length(conserved_signif_tads_genes) + length(not_conserved_signif_tads_genes) + length(not_signif_tads_genes) == length(geneList)

      list(
        universe_genes = as.character(geneList),
        conserved_signif_tads_genes = as.character(conserved_signif_tads_genes),
        not_conserved_signif_tads_genes = as.character(not_conserved_signif_tads_genes) #,
      )
    }
    names(exprds_list) <- file.path(hicds, all_exprds[[paste0(hicds)]])
    exprds_list
  }
  all_gene_list <- unlist(all_gene_list, recursive=FALSE)
  stopifnot(length(all_gene_list) == length(all_datasets))
  outFile <- file.path(outFolder, paste0(file_prefix, "all_gene_list.Rdata"))
  save(all_gene_list, file = outFile, version=2)
} else {
  outFile <- file.path(outFolder, paste0(file_prefix, "all_gene_list.Rdata"))
  cat("... load data\n")
  all_gene_list <- get(load(outFile))
}

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

padjVarGO <- "p.adjust" # p.adjust or qvalue ???

logFile <- file.path(outFolder, "go_signif_across_hicds_logFile.txt")
if(buildData) file.remove(logFile)

barplot_vars <- c("foldEnrichment", "geneRatio", "log10_pval")
barplot_vars_tit <- setNames(c("Fold enrichment", "Gene ratio", paste0("-log10(", padjVarGO,  ")")), barplot_vars)
barplot_vars <- c( "log10_pval")
barplot_vars_tit <- setNames(c(paste0("-log10(", padjVarGO,  ")")), barplot_vars)

plotMaxBars <- 10

par_bot <- 12

padjVarGO_plotThresh <- 0.05

# GO for BP nad MF [do not take c5_CC]
if(enricher_ontologyType == "BP" | enricher_ontologyType == "MF" | enricher_ontologyType == "BP_MF"){
  gmtFile <- file.path(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2", paste0("c5.", tolower(enricher_ontologyType), ".v6.1.entrez.gmt"))
} else {
  stop(paste0(enricher_ontologyType, " is not a valid ontologyType\n"))
}
stopifnot(file.exists(gmtFile))
c5_msigdb <- read.gmt(gmtFile)


printAndLog <- function(txt, logFile=NULL) {
  cat(txt)
  cat(txt, file=logFile, append=T)
}


txt <- paste0("!!! HARD CODED SETTINGS !!! \n")
printAndLog(txt, logFile)
txt <- paste0("... gmtFile:\t", gmtFile, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher background genes:\t", "universe", "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher pvalueCutoff:\t", enricher_pvalueCutoff, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher pAdjustMethod:\t", enricher_pAdjustMethod, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher minGSSize:\t", enricher_minGSSize, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher maxGSSize:\t", enricher_maxGSSize, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher qvalueCutoff:\t", enricher_qvalueCutoff, "\n")
printAndLog(txt, logFile)
txt <- paste0("... enricher ontologyType:\t", enricher_ontologyType, "\n")
printAndLog(txt, logFile)
# txt <- paste0("... pval thresh select GOs:\t", pvalSelectGO, "\n")
# printAndLog(txt, logFile)
txt <- paste0("... enricher_results_sortGOby:\t", enricher_results_sortGOby, "\n")
printAndLog(txt, logFile)

txt <- paste0("... padjVarGO_plotThresh:\t", padjVarGO_plotThresh, "\n")
printAndLog(txt, logFile)
txt <- paste0("... plotMaxBars:\t", plotMaxBars, "\n")
printAndLog(txt, logFile)


all_universe_genes <- unique(unlist(lapply(all_gene_list, function(x) x[["universe_genes"]])))
all_conserved_signif_tads_genes <- unique(unlist(lapply(all_gene_list, function(x) x[["conserved_signif_tads_genes"]]))) 
all_not_conserved_signif_tads_genes <- unique(unlist(lapply(all_gene_list, function(x) x[["not_conserved_signif_tads_genes"]])))

go_all_universe_genes <- as.character(all_universe_genes)[as.character(all_universe_genes) %in% as.character(c5_msigdb$gene)]
go_all_conserved_signif_tads_genes <- as.character(all_conserved_signif_tads_genes)[as.character(all_conserved_signif_tads_genes) %in% as.character(c5_msigdb$gene)]
go_all_not_conserved_signif_tads_genes <- as.character(all_not_conserved_signif_tads_genes)[as.character(all_not_conserved_signif_tads_genes) %in% as.character(c5_msigdb$gene)]

cat("... available annot. for all_universe_genes:\t", length(go_all_universe_genes) , "/", length(all_universe_genes), "\n")
cat("... available annot. for all_conserved_signif_tads_genes:\t", length(go_all_conserved_signif_tads_genes) , "/", length(all_conserved_signif_tads_genes), "\n")
cat("... available annot. for all_not_conserved_signif_tads_genes:\t", length(go_all_not_conserved_signif_tads_genes) , "/", length(all_not_conserved_signif_tads_genes), "\n")

stopifnot(go_all_conserved_signif_tads_genes %in% go_all_universe_genes)
stopifnot(go_all_not_conserved_signif_tads_genes %in% go_all_universe_genes)

if(! plotOnly) {
  
  
  #***** 1) conserved_signif
  cat(paste0(">  start enricher for conserved_signif \n"))
  
  if(length(go_all_conserved_signif_tads_genes) > 0) {
    
    conserved_signif_enrich <- enricher(gene = go_all_conserved_signif_tads_genes, 
                                        TERM2GENE=c5_msigdb,
                                        universe = go_all_universe_genes,
                                        pvalueCutoff = enricher_pvalueCutoff, 
                                        pAdjustMethod = enricher_pAdjustMethod, 
                                        minGSSize = enricher_minGSSize, 
                                        maxGSSize = enricher_maxGSSize, 
                                        qvalueCutoff =enricher_qvalueCutoff)
    
    conserved_signif_enrich_resultDT <- conserved_signif_enrich@result
    conserved_signif_enrich_resultDT <- conserved_signif_enrich_resultDT[order(conserved_signif_enrich_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
    conserved_signif_enrich_resultDT$log10_pval <- -log10(conserved_signif_enrich_resultDT[,paste0(padjVarGO)])
    conserved_signif_enrich_resultDT$geneRatio <- as.numeric(sapply(conserved_signif_enrich_resultDT$GeneRatio, function(x) {
      gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
    conserved_signif_enrich_resultDT$bgRatio <- as.numeric(sapply(conserved_signif_enrich_resultDT$BgRatio, function(x) {
      gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
    conserved_signif_enrich_resultDT$foldEnrichment <- conserved_signif_enrich_resultDT$geneRatio/conserved_signif_enrich_resultDT$bgRatio
    
    genes_signif_plotMax <- min(c(plotMaxBars, nrow(conserved_signif_enrich_resultDT)))
    if(genes_signif_plotMax > 0) {
      for(var_plot in barplot_vars) {
        myTit <- paste0(barplot_vars_tit[var_plot], " ", data_cmpType, " (conserved_signif)")
        conserved_signif_enrich_resultDT <- conserved_signif_enrich_resultDT[order(conserved_signif_enrich_resultDT[,var_plot], decreasing=TRUE),]
        
        my_x <- gsub("GO_", "", conserved_signif_enrich_resultDT$Description[1:genes_signif_plotMax])
        
        #### UPDATE HERE !!!!
        conserved_signif_dt <- conserved_signif_enrich_resultDT[,c(var_plot, "Description"), drop=FALSE]
        # conserved_signif_dt <- conserved_signif_enrich_resultDT[1:genes_signif_plotMax,c(var_plot, "Description"), drop=FALSE]
        conserved_signif_dt$labs <-  gsub("GO_", "", conserved_signif_dt$Description)
        conserved_signif_dt$plot_labs <-  unlist(lapply(strwrap(gsub("_", " ", conserved_signif_dt$labs), 
                                                                width = strwdth, simplify=FALSE), 
                                                        function(x) paste0(x, collapse="\n")))
        
        outFile <- file.path(outFolder,paste0(file_prefix, "conserved_signif", var_plot, "_GO_",
                                              enricher_ontologyType, "_", var_plot, "_barplot_dt.Rdata"))
        save(conserved_signif_dt, file=outFile, version=2)
        
        
        outFile <- file.path(outFolder,paste0(file_prefix, "conserved_signif", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
        do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.5))    
        par(oma=c(par_bot,1,1,1))
        
        
        
        barplot(conserved_signif_enrich_resultDT[1:genes_signif_plotMax,var_plot],
                main = myTit, 
                # ylab = paste0("-log10 (", padjVarGO, ")"),
                ylab = paste0(barplot_vars_tit[var_plot]),
                names.arg =unlist(lapply(strwrap(gsub("_", " ", my_x), width = strwdth, 
                                                 simplify=FALSE), function(x) paste0(x, collapse="\n"))), 
                las=2,
                cex.lab = plotCex,
                cex.main = plotCex,
                cex.axis=plotCex,
                cex.names = 0.9
        )
        mtext(side=3, text=paste0("(# datasets = ", nDS, ")"), font=3)
        foo <- dev.off()
        cat(paste0("... written: ", outFile, "\n"))
      }
    }
    conserved_signif_enrich_resultDT <- data.frame(conserved_signif_enrich_resultDT)
    
  } else {
    conserved_signif_enrich_resultDT <- NULL
  }
  
  
  #***** 2) not_conserved_signif
  
  cat(paste0("> start enricher for not_conserved_signif \n"))
  
  if(length(go_all_not_conserved_signif_tads_genes) > 0) {
    
    not_conserved_signif_enrich <- enricher(gene = go_all_not_conserved_signif_tads_genes, 
                                            TERM2GENE=c5_msigdb,
                                            universe = go_all_universe_genes,
                                            pvalueCutoff = enricher_pvalueCutoff, 
                                            pAdjustMethod = enricher_pAdjustMethod, 
                                            minGSSize = enricher_minGSSize, 
                                            maxGSSize = enricher_maxGSSize, 
                                            qvalueCutoff =enricher_qvalueCutoff)
    
    not_conserved_signif_enrich_resultDT <- not_conserved_signif_enrich@result
    not_conserved_signif_enrich_resultDT <- not_conserved_signif_enrich_resultDT[order(not_conserved_signif_enrich_resultDT[,enricher_results_sortGOby], decreasing=FALSE), ]
    not_conserved_signif_enrich_resultDT$log10_pval <- -log10(not_conserved_signif_enrich_resultDT[,paste0(padjVarGO)])
    not_conserved_signif_enrich_resultDT$geneRatio <- as.numeric(sapply(not_conserved_signif_enrich_resultDT$GeneRatio, function(x) {
      gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
    not_conserved_signif_enrich_resultDT$bgRatio <- as.numeric(sapply(not_conserved_signif_enrich_resultDT$BgRatio, function(x) {
      gr <- as.numeric(eval(parse(text = x ))); stopifnot(!is.na(gr)); gr})) 
    not_conserved_signif_enrich_resultDT$foldEnrichment <- not_conserved_signif_enrich_resultDT$geneRatio/not_conserved_signif_enrich_resultDT$bgRatio
    
    genes_signif_plotMax <- min(c(plotMaxBars, nrow(not_conserved_signif_enrich_resultDT)))
    if(genes_signif_plotMax > 0) {
      for(var_plot in barplot_vars) {
        myTit <- paste0(barplot_vars_tit[var_plot], " ", data_cmpType, " (not_conserved_signif)")
        not_conserved_signif_enrich_resultDT <- not_conserved_signif_enrich_resultDT[order(not_conserved_signif_enrich_resultDT[,var_plot], decreasing=TRUE),]
        
        my_x <- gsub("GO_", "", not_conserved_signif_enrich_resultDT$Description[1:genes_signif_plotMax])
        
        
        save(not_conserved_signif_enrich_resultDT, file="not_conserved_signif_enrich_resultDT.Rdata", version=2)
        save(barplot_vars_tit, file="barplot_vars_tit.Rdata", version=2)
        save(my_x, file="my_x.Rdata", version=2)
        
        #### changed here !!!
        not_conserved_signif_dt <- not_conserved_signif_enrich_resultDT[,c(var_plot, "Description"),drop=FALSE]
        # not_conserved_signif_dt <- not_conserved_signif_enrich_resultDT[1:genes_signif_plotMax,c(var_plot, "Description"),drop=FALSE]
        not_conserved_signif_dt$labs <-  gsub("GO_", "", not_conserved_signif_dt$Description)
        not_conserved_signif_dt$plot_labs <-  unlist(lapply(strwrap(gsub("_", " ", not_conserved_signif_dt$labs), 
                                                                    width = strwdth, simplify=FALSE), 
                                                            function(x) paste0(x, collapse="\n")))
        
        outFile <- file.path(outFolder,paste0(file_prefix, "not_conserved_signif", var_plot, "_GO_",
                                              enricher_ontologyType, "_", var_plot, "_barplot_dt.Rdata"))
        save(not_conserved_signif_dt, file=outFile, version=2)
        
        
        
        outFile <- file.path(outFolder,paste0(file_prefix, "not_conserved_signif", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplot", ".", plotType))
        do.call(plotType, list(outFile, height = myHeight*1.2, width = myWidth*1.5))    
        
        
        par(oma=c(par_bot,1,1,1))
        barplot(not_conserved_signif_enrich_resultDT[1:genes_signif_plotMax,var_plot],
                main = myTit, 
                # ylab = paste0("-log10 (", padjVarGO, ")"),
                ylab = paste0(barplot_vars_tit[var_plot]),
                names.arg =unlist(lapply(strwrap(gsub("_", " ", my_x), width = strwdth, simplify=FALSE), 
                                         function(x) paste0(x, collapse="\n"))),
                # names.arg = "",
                las=2,
                cex.main=plotCex,
                cex.lab = plotCex,
                cex.axis=plotCex,
                cex.names = 0.9
        )
        # 
        # axis(1, at = 1:genes_signif_plotMax,
        #      labels = unlist(lapply(strwrap(gsub("_", " ", my_x), width = 20, simplify=FALSE), function(x) paste0(x, collapse="\n"))),)
        # 
        mtext(side=3, text=paste0("(# datasets = ", nDS, ")"), font=3)
        foo <- dev.off()
        cat(paste0("... written: ", outFile, "\n"))
      }
    }
    not_conserved_signif_enrich_resultDT <- data.frame(not_conserved_signif_enrich_resultDT)
    # not_conserved_signif_enrich_resultDT$hicds <- hicds
    # not_conserved_signif_enrich_resultDT$exprds <- exprds
    
    
  } else {
    not_conserved_signif_enrich_resultDT <- NULL
  }
  
  
  outFile <- file.path(outFolder,paste0(file_prefix, "not_conserved_signif_enrich_resultDT.Rdata"))
  save(not_conserved_signif_enrich_resultDT, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  outFile <- file.path(outFolder,paste0(file_prefix, "conserved_signif_enrich_resultDT.Rdata"))
  save(conserved_signif_enrich_resultDT, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
  
  
} else {
  
  for(var_plot in barplot_vars) {
      
    # inFile <- file.path(outFolder,paste0(file_prefix, "not_conserved_signif", var_plot, "_GO_",
    #                                       enricher_ontologyType, "_", var_plot, "_barplot_dt.Rdata"))
    # not_conserved_signif_dt <- get(load(inFile))
    
    inFile <- file.path(outFolder,paste0(file_prefix, "conserved_signif", var_plot, "_GO_",
                                         enricher_ontologyType, "_", var_plot, "_barplot_dt.Rdata"))
    conserved_signif_dt <- get(load(inFile))
    conserved_signif_dt <- conserved_signif_dt[order(conserved_signif_dt$log10_pval, decreasing=TRUE),]
    conserved_signif_dt$plot_labs <- factor(conserved_signif_dt$plot_labs, levels=as.character(conserved_signif_dt$plot_labs))
    
    myTit <- paste0(barplot_vars_tit[var_plot], " ", data_cmpType, " (conserved_signif)")
    myTit <- "Signif. enriched GO from conserved region"
    
    subTit <- paste0("(conserved in >= ", conservThresh, "/", nDS,  " datasets)")
    
    if(var_plot == "log10_pval") {
      my_ylab <- "adj. p-val [-log10]"
    } else {
      my_ylab <- var_plot
    }
    
    
    strwdth <- 25
    
    conserved_signif_dt$plot_labs <-  unlist(lapply(strwrap(gsub("_", " ", conserved_signif_dt$labs), 
                                                                width = strwdth, simplify=FALSE), 
                                                        function(x) paste0(x, collapse="\n")))
    
    
   ggbar_p <-  ggbarplot(conserved_signif_dt, 
              x="plot_labs", y=var_plot, fill="darkgrey") +
      ggtitle(myTit, subtitle=subTit)+
      my_box_theme +
      labs(x="" , y = my_ylab) + 
      theme(
        axis.text.x = element_text(hjust=1, vjust=0.5,size=10,angle=90)
    ) +
     scale_y_continuous( breaks= scales::pretty_breaks(n = 8))
     
   outFile <- file.path(outFolder,paste0(file_prefix, "conserved_signif", var_plot, "_GO_", enricher_ontologyType, "_", var_plot, "_barplotGG", ".", plotType))
   ggsave(ggbar_p, filename = outFile, height=myHeightGG, width=myWidthGG*1.2)
   cat(paste0("... written: ", outFile, "\n"))


topSignifGO_conservedRegions_dt <- conserved_signif_dt
saveFile <- file.path(outFolder, paste0("fig4B_topSignifGO_conservedRegions_dt.Rdata"))
save(topSignifGO_conservedRegions_dt, file=saveFile, version=2)
cat(paste0("... written:" , saveFile, "\n"))

   
   
  }
  
  
}


######################################################################################
######################################################################################
######################################################################################
cat(paste0("... written: ", logFile, "\n"))
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


  
  
