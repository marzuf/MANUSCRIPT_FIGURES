# Rscript check_missing.R

outFolder <-  "CHECK_MISSING"
dir.create(outFolder, recursive = TRUE)

require(foreach)


mainFolder <- "../..//v2_Yuanlong_Cancer_HiC_data_TAD_DA"

minCount <- 5

hicds <- "ENCSR489OCU_NCI-H460_40kb"
exprds1 <- "TCGAluad_norm_luad"

all_conds <- c("lusc", "luad")

setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)
stopifnot(!duplicated(gff_dt$entrezID))
entrez2symb <- setNames(gff_dt$symbol, gff_dt$entrezID)



plot_dt <- foreach(cond = all_conds, .combine='rbind')  %dopar% {
  
  exprds1 <- paste0("TCGA", cond, "_norm_", cond)
  
  gl1 <- get(load(file.path(mainFolder, "PIPELINE/OUTPUT_FOLDER", hicds, exprds1,"0_prepGeneDataCheckRm", "pipeline_geneList_noCpmFilter.Rdata")))
  
  gl1_kept <- get(load(file.path(mainFolder, "PIPELINE/OUTPUT_FOLDER", hicds, exprds1,"0_prepGeneDataCheckRm", "pipeline_geneList.Rdata")))
  
  
  stopifnot(gl1 %in% names(entrez2symb))
  
  matchNames_dt1 <- data.frame(
    entrezRef = gl1,
    entrezRNA = names(gl1),
    symbol = entrez2symb[gl1],
    kept = gl1 %in% gl1_kept,
    stringsAsFactors = FALSE
  ) 
  
  
  stopifnot(!is.na(matchNames_dt1 ))
  
  
  matchNames_dt1_mmp <- matchNames_dt1[grepl("^MMP", matchNames_dt1$symbol),]
  
  exprds1_e2s <- setNames(matchNames_dt1_mmp$symbol, matchNames_dt1_mmp$entrezRNA)
  exprds1_kept <- setNames(matchNames_dt1_mmp$kept, matchNames_dt1_mmp$entrezRNA)
  
  exprds1_countdt <- get(load(file.path(mainFolder, "PIPELINE/OUTPUT_FOLDER", hicds, exprds1,"0_prepGeneDataCheckRm", "countFilter_rnaseqDT.Rdata")))
  exprds1_countdt[1:5,1:5]
  
  stopifnot(rownames(exprds1_countdt) %in% matchNames_dt1$entrezRNA)
  
  
  exprds1_countdt_mmp <- exprds1_countdt[rownames(exprds1_countdt) %in% matchNames_dt1_mmp$entrezRNA,]
  
  s1 <- get(load(file.path(setDir, paste0("//mnt/ed4/marie/other_datasets/TCGA", cond, "_norm_", cond, "/", "norm_ID.Rdata"))))
  s2 <- get(load(file.path(setDir, paste0("//mnt/ed4/marie/other_datasets/TCGA", cond, "_norm_", cond, "/", cond, "_ID.Rdata"))))
  
  
  dt1_sums <- apply(exprds1_countdt_mmp, 1, function(x) sum(x >= minCount)) 
  
  totSamp <- length(s1) + length(s2)
  
  stopifnot(dt1_sums <= totSamp)
  
  exprds1_sum_dt <- data.frame(
    exprds = exprds1,
    entrezRNA = names(dt1_sums),
    nSampMinCount = dt1_sums,
    nSampMinRatio = dt1_sums/totSamp,
    stringsAsFactors = FALSE
  )
  stopifnot(exprds1_sum_dt$entrezRNA %in% names(exprds1_e2s))
  exprds1_sum_dt$symbol <- exprds1_e2s[exprds1_sum_dt$entrezRNA]
  
  stopifnot(exprds1_sum_dt$entrezRNA %in% names(exprds1_kept))
  exprds1_sum_dt$kept <- exprds1_kept[exprds1_sum_dt$entrezRNA]
  
  
  stopifnot(!is.na(exprds1_sum_dt))
  
  exprds1_sum_dt
  
}

outFile <- file.path("CHECK_MISSING", "plot_dt.Rdata")
save(plot_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

# stop("_ok")

# plot_dt <- get(load("CHECK_MISSING/plot_dt.Rdata"))


# boxplot(nSampMinCount ~ exprds + symbol, data=plot_dt)

require(ggplot2)
# p1 <- ggplot(data= plot_dt, aes( x= symbol, y = nSampMinCount, color = exprds)) + 
#   geom_bar(stat="identity", position="dodge")

anyKept <- unique(plot_dt$symbol[plot_dt$kept])

# p2 <- ggplot(data= plot_dt[plot_dt$symbol %in% anyKept,], aes( x= symbol, y = nSampMinCount, color = kept, fill  = exprds)) + 
#   geom_bar(stat="identity", position="dodge")+
#   scale_fill_manual(values=c("red", "blue")) +
#   scale_color_manual(values=c("green", "yellow"))

p3 <- ggplot(data= plot_dt[plot_dt$symbol %in% anyKept,], aes( x= symbol, y = nSampMinRatio, color = kept, fill  = exprds)) + 
  geom_bar(stat="identity", position="dodge")+
  scale_fill_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("green", "yellow")) +
  geom_hline(yintercept=0.8)+
  theme(
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
    axis.title.x=element_blank()
  )
outFile <- file.path(outFolder, paste0("nSampMinRatio_", hicds, "_luad_lusc.svg"))
ggsave(p3, filename = outFile , height=5, width=7)
cat(paste0("... written: ", outFile, "\n"))

