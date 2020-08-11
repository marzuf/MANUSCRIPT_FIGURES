options(scipen=100)

SSHFS=F

setDir <- "/media/electron"
setDir <- ""

# Rscript ds_conservation_acrossDS_vSameNbr_purityFilter_final.R

purity_ds <- "aran"
pm <- "CPE"
purity_plot_name <- "Aran - CPE"
transfExpr="log10"

hicds="K562_40kb"
exprds="TCGAlaml_wt_mutFLT3"

script_name <- "ds_conservation_acrossDS_vSameNbr_purityFilter_final.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(foreach)
require(doMC)
require(RColorBrewer)
require(reshape2)
require(ggplot2)
require(ggpubr)
# require(gplots)
registerDoMC(4)
source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")

args <- commandArgs(trailingOnly = T)

if(length(args) == 0){
  cmpType <- ""
  filePrefix <- ""
  cmpTit <- paste0("all")
} else {
  cmpType <- args[1]  
  filePrefix <- paste0(cmpType, "_")
  cmpTit <- cmpType
  stopifnot(cmpType %in% c("subtypes", "wt_vs_mut", "norm_vs_tumor"))
}

outFolder <- file.path("DS_CONSERVATION_ACROSSDS_VSAMEBR_PURITYFILTER_FINAL", cmpType, purity_ds, pm, transfExpr)
dir.create(outFolder, recursive=TRUE)

entrez2symb_dt <- read.delim(file.path(setDir,
                                       "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"),
                             header=T, stringsAsFactors = FALSE)
entrez2symb_dt$entrezID <- as.character(entrez2symb_dt$entrezID)
entrez2symb <- setNames(entrez2symb_dt$symbol, entrez2symb_dt$entrezID)



buildTable <- TRUE

plotType <- "svg"

source("../settings.R")

result_dt <- get(load(file.path(runFolder, "CREATE_FINAL_TABLE/all_result_dt.Rdata")))
nDS <- length(unique(file.path(result_dt$hicds, result_dt$exprds)))

myHeightGG <- 5
myWidthGG <- 7

tad_pval <- 0.01
minMatchBp_ratio <- 0.8
minMatch_genes <- 3

atLeast <- 10
atMin <- 2

all_datasets <- sort(unique(file.path(result_dt$hicds, result_dt$exprds)))

# retrieve the max from permut
all_perm_dt <- get(load(file.path(runFolder,
                      "CONSERV_SIGNIF_OBSSHUFFLE_PURITYFILTER_VSAMENBR_FINAL", purity_ds, pm, transfExpr, "all_perm_dt.Rdata")))
max_perm <- max(all_perm_dt$conserved)

permline_col <- "black"

add_maxPerm_line <- function(p) {
   return(p + 
    geom_hline(yintercept=max_perm, linetype=2, color=permline_col) + 
    annotate("text", x = 1, y= max_perm+0.5, 
              label=paste0("max permut. (vSameNbr) conserved = ", max_perm),
              hjust=0, vjust=0, 
              fontface ="italic", color=permline_col))
}

inFile <- file.path(runFolder,
                    "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER_FINAL",
                    cmpType,purity_ds, pm, transfExpr,
                    paste0("conserved_signif_tadsadjPvalComb", tad_pval, "_minBpRatio", minMatchBp_ratio, "_minInterGenes", minMatch_genes, ".Rdata"))

obs_data <- get(load(inFile))

all_nCons <- lengths(obs_data)
dsByReg_dt <- do.call(rbind, lapply(obs_data, function(x) as.numeric(all_datasets %in% unique(dirname(x)))))
colnames(dsByReg_dt) <- all_datasets
stopifnot(sum(dsByReg_dt) == sum(all_nCons))
stopifnot(max(rowSums(dsByReg_dt)) == max(all_nCons))
obsCons_dt <- dsByReg_dt

m_dt <- melt(obsCons_dt)

colnames(m_dt) <- c("region", "dataset", "conserved")
stopifnot(grepl("region", m_dt$region))
m_dt$dataset <- as.character(m_dt$dataset)
m_dt$hicds <- dirname(m_dt$dataset)
m_dt$exprds <- basename(m_dt$dataset)
m_dt$cmpType <- all_cmps[paste0(m_dt$exprds)]
m_dt$cmpCol <- all_cols[m_dt$cmpType]

agg_dt <- aggregate(conserved ~ region, data=m_dt, FUN=sum)
agg_dt <- agg_dt[order(agg_dt$conserved),]
reg_levels <- as.character(agg_dt$region)
stopifnot(agg_dt$conserved >= atMin)
# agg_dt$region_rank <-

# plot(
#   x = 1:nrow(agg_dt),
#   y = agg_dt$conserved,
#   type="b"
# )

inFile <- file.path(runFolder,
                    "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER_FINAL",
                    purity_ds, pm, transfExpr,
                    paste0("conserved_regions_with_genes_signif_tadsadjPvalComb", tad_pval, "_minBpRatio", minMatchBp_ratio, "_minInterGenes", minMatch_genes, ".Rdata"))
content_dt <- get(load(inFile))
agg_dt2 <- agg_dt
colnames(agg_dt2)[colnames(agg_dt2) == "region"] <- "conserved_region"
print_dt <- merge(content_dt, agg_dt2, by=c("conserved_region"), all.x=F, all.y=T)
colnames(print_dt)[colnames(print_dt) == "conserved"] <- "nConserved"
print_dt <- print_dt[order(print_dt$nConserved, decreasing = T),]

outFile <- file.path(outFolder, paste0(filePrefix, "conserved_regions_table.txt"))
write.table(print_dt, file = outFile, sep="\t", col.names =T, row.names = F, quote=F)
cat(paste0("... written: ", outFile, "\n"))


agg_cmp_dt <- aggregate(conserved ~ region + cmpType, data=m_dt, FUN=sum)
agg_cmp_dt$cmpCol <- all_cols[agg_cmp_dt$cmpType]


conserv_atLeast <- agg_dt$region[agg_dt$conserved >= atLeast]
reg_levels2 <- reg_levels[reg_levels %in% conserv_atLeast]

#######################################################
# plot all
#######################################################

plot_dt <- agg_cmp_dt

plot_dt$region <- factor(plot_dt$region, levels = reg_levels)
stopifnot(!is.na(plot_dt$region))

plot_dt$cmpType_lab <- cmp_names[paste0(plot_dt$cmpType)]
stopifnot(!is.na(plot_dt$cmpType_lab))

all_cols_lab <- all_cols
names(all_cols_lab) <- cmp_names[names(all_cols_lab)]

cmp_levels <- sort(cmp_names)

stopifnot(cmp_levels %in% names(all_cols_lab))

plot_dt$region_rank <- as.numeric(plot_dt$region)

nCons <- length(unique(plot_dt$region))

my_ylab <- "# datasets conserved"
# my_xlab <- "Ranked conserved regions"
my_xlab <- paste0("Ranked conserved regions (", nCons, " with >= ", atMin, " DS cons.)")
myTit <- "Signif. conserved regions"
subTit <- paste0(purity_plot_name , " data - TAD adj. p-val <= ", tad_pval, "\nconserv. match ratio bp >= ", minMatchBp_ratio, " and match genes >= ", minMatch_genes)


plot_dt$cmpType_lab <- factor(plot_dt$cmpType_lab, levels=cmp_levels)
stopifnot(!is.na(plot_dt$cmpType_lab))

save(plot_dt, file=file.path(outFolder, "plot_dt.Rdata"), version=2)

conservedRegions_dt <- plot_dt
saveFile <- file.path(outFolder, paste0("fig4A_conservedRegions_dt.Rdata"))
save(conservedRegions_dt, file=saveFile, version=2)
cat(paste0("... written:" , saveFile, "\n"))

bar_p <- ggplot(plot_dt, aes(x=region_rank, y=conserved, fill = cmpType_lab, color=cmpType_lab)) +
  ggtitle(myTit, subtitle = subTit)+
  geom_bar(stat="identity") +
  scale_y_continuous(name=my_ylab, breaks= scales::pretty_breaks(n = 5))+
  scale_x_continuous(name=my_xlab, breaks= noZero_breaks(n=10), expand=c(0,0))+
  scale_fill_manual(values=all_cols_lab[cmp_levels], labels=cmp_levels)+
  scale_color_manual(values=all_cols_lab[cmp_levels], labels=cmp_levels)+
  labs(x=my_xlab, y=my_ylab, color="", fill="")+
  my_box_theme +
  theme(
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x  =  element_blank(),
    panel.grid.minor.y  =  element_blank(),
    legend.text = element_text(size=12),
    # axis.line = element_line()
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    axis.line = element_line()
  )



bar_p <- add_maxPerm_line(bar_p)

outFile <- file.path(outFolder, paste0(filePrefix, "nConserved_byRegion_allCond.", plotType))
ggsave(bar_p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))

# stop("-ok\n")

#######################################################
# plot only from at least top 10
#######################################################
plot_dt <- agg_cmp_dt[agg_cmp_dt$region %in% conserv_atLeast,]

plot_dt$region <- factor(plot_dt$region, levels = reg_levels2)
stopifnot(!is.na(plot_dt$region))

plot_dt$cmpType_lab <- cmp_names[paste0(plot_dt$cmpType)]
stopifnot(!is.na(plot_dt$cmpType_lab))

all_cols_lab <- all_cols
names(all_cols_lab) <- cmp_names[names(all_cols_lab)]

cmp_levels <- sort(cmp_names)

stopifnot(cmp_levels %in% names(all_cols_lab))

plot_dt$region_rank <- as.numeric(plot_dt$region)

nCons <- length(unique(plot_dt$region))

my_ylab <- "# datasets conserved"
# my_xlab <- "Ranked conserved regions"
my_xlab <- paste0("Ranked conserved regions (", nCons, " with >= ", atLeast, " DS cons.)")
myTit <- "Signif. conserved regions"
subTit <- paste0(purity_plot_name, " data - TAD adj. p-val <= ", tad_pval, "\nconserv. match ratio bp >= ", minMatchBp_ratio, " and match genes >= ", minMatch_genes)


plot_dt$cmpType_lab <- factor(plot_dt$cmpType_lab, levels=cmp_levels)
stopifnot(!is.na(plot_dt$cmpType_lab))

save(plot_dt, file=file.path(outFolder, "plot_dt.Rdata"), version=2)

bar_p <- ggplot(plot_dt, aes(x=region_rank, y=conserved, fill = cmpType_lab, color=cmpType_lab)) +
  ggtitle(myTit, subtitle = subTit)+
  geom_bar(stat="identity") +
  scale_y_continuous(name=my_ylab, breaks= scales::pretty_breaks(n = 5))+
  scale_x_continuous(name=my_xlab, breaks= noZero_breaks(n=5), expand=c(0,0))+
  scale_fill_manual(values=all_cols_lab[cmp_levels], labels=cmp_levels)+
  scale_color_manual(values=all_cols_lab[cmp_levels], labels=cmp_levels)+
  labs(x=my_xlab, y=my_ylab, color="", fill="")+
  my_box_theme +
  theme(
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x  =  element_blank(),
    legend.text = element_text(size=12),
    # axis.line = element_line()
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    axis.line.x = element_line()
  )

outFile <- file.path(outFolder, paste0(filePrefix, "nConserved_byRegion_allCond_atLeast", atLeast, ".", plotType))
ggsave(bar_p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



#####

allCmpTypes <- unique(all_cmps)
cmpT=allCmpTypes[1]
byCmp_dt <- foreach(cmpT = allCmpTypes, .combine='rbind') %dopar% {
  inFile <- file.path(runFolder,
                      "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2_PURITYFILTER_FINAL",
                      cmpT, purity_ds, pm, transfExpr,
                      paste0("plot_matching_dt_signif_tadsadjPvalComb", tad_pval, "_minBpRatio", minMatchBp_ratio, "_minInterGenes", minMatch_genes, ".Rdata"))
  
  matching_dt <- get(load(inFile))
  
  dt <- t(matching_dt)
  m_dt <- melt(dt)
  colnames(m_dt) <- c("region", "dataset", "conserved")
  m_dt$dataset <- as.character(m_dt$dataset)
  m_dt$hicds <- dirname(m_dt$dataset)
  m_dt$exprds <- basename(m_dt$dataset)
  m_dt$cmpType <- all_cmps[paste0(m_dt$exprds)]
  stopifnot(as.character(m_dt$cmpType) == cmpT)
  agg_dt <- aggregate(conserved ~ region, data=m_dt, FUN=sum)
  agg_dt <- agg_dt[order(agg_dt$conserved),]
  reg_levels <- as.character(agg_dt$region)
  stopifnot(agg_dt$conserved >= 2)
  agg_dt$cmpType <- cmpT
  agg_dt$reg_rank <- 1:nrow(agg_dt)
  agg_dt
}
# byCmp_dt$reg_rank <- 1:nrow(byCmp_dt)
# byCmp_dt$reg_rank <- factor(as.character(byCmp_dt$reg_rank))

my_ylab <- "# datasets conserved"
my_xlab <- "Ranked conserved regions"
myTit <- "Signif. conserved regions"
subTit <- paste0(purity_plot_name, " data - TAD adj. p-val <= ", tad_pval, "\nconserv. match ratio bp >= ", minMatchBp_ratio, " and match genes >= ", minMatch_genes)

byCmp_dt$cmpType_lab <- cmp_names[paste0(byCmp_dt$cmpType)]
stopifnot(!is.na(byCmp_dt$cmpType_lab))

byCmp_dt$cmpType_lab <- factor(byCmp_dt$cmpType_lab, levels=cmp_levels)
stopifnot(!is.na(byCmp_dt$cmpType_lab))

save(byCmp_dt, file=file.path(outFolder, "byCmp_dt.Rdata"), version=2)


p_wrap <- ggplot(byCmp_dt, aes(x=reg_rank, y=conserved, group=1))+
  ggtitle(myTit, subtitle = subTit)+
  geom_line()+
  facet_wrap(. ~ cmpType_lab, scales="free") +
  # scale_y_continuous(name=my_ylab)+
  scale_y_continuous(name=my_ylab, breaks= scales::pretty_breaks(n = 5))+
  scale_x_continuous(name=my_xlab, breaks= noZero_breaks(n=5))+
  my_box_theme +
  theme(
    strip.text = element_text(size=10, face="bold"),# color=all_cols_lab),
    # strip.background = element_rect(color=all_cols_lab),
    panel.grid.major.x =  element_blank(),
    panel.grid.minor.x  =  element_blank(),
    # axis.line = element_line()
    # axis.ticks.x = element_blank(),
    # axis.text.x = element_blank(),
    axis.line.x = element_line()
  )


outFile <- file.path(outFolder, paste0(filePrefix, "nConserved_byRegion_wrapCond1.", plotType))
ggsave(p_wrap, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



g <- ggplot_gtable(ggplot_build(p_wrap))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- all_cols_lab[cmp_levels]
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  # curr_sub <- g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill
  # g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- all_cols_lab[names(curr_sub)]
  k <- k+1
}

outFile <- file.path(outFolder, paste0(filePrefix, "nConserved_byRegion_wrapCond.", plotType))
do.call(plotType, list(outFile, height=myHeightGG, width=myWidthGG))
grid::grid.draw(g)
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

######################################################################################
######################################################################################
######################################################################################

# all_dt <- get(load("DS_CONSERVATION_ACROSSDS/fig4A_conservedRegions_dt.Rdata"))
# agg_all_dt <- aggregate(conserved~region, FUN=sum, data=all_dt)
# 
# cpe_dt <- get(load("DS_CONSERVATION_ACROSSDS_PURITYFILTER_FINAL/CPE/log10/fig4A_conservedRegions_dt.Rdata"))
# cpe_agg_all_dt <- aggregate(conserved~region, FUN=sum, data=cpe_dt)
# 
# 
# epic_dt <- get(load("DS_CONSERVATION_ACROSSDS_PURITYFILTER/EPIC/log10/fig4A_conservedRegions_dt.Rdata"))
# epic_agg_all_dt <- aggregate(conserved~region, FUN=sum, data=epic_dt)
# 
# aran_dt <- get(load("DS_CONSERVATION_ACROSSDS_PURITYFILTER/log10/fig4A_conservedRegions_dt.Rdata"))
# aran_agg_all_dt <- aggregate(conserved~region, FUN=sum, data=aran_dt)
# 
# source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")
# 
# png("cmp_nDS_conserved_purity_filter.png", height=400, width=600)
# plot_multiDens(
#   list(noFilt = agg_all_dt$conserved,
#        aranCpeFilt = cpe_agg_all_dt$conserved,
#        aranFilt = aran_agg_all_dt$conserved,
#        epicFilt = epic_agg_all_dt$conserved
#        ), my_xlab = "# datasets conserved", plotTit = "Recurrent DA TADs - # datasets conservation"
# )
# foo <- dev.off()
                                                
######################################################################################
######################################################################################
######################################################################################
cat("*** DONE\n")
cat(paste0(startTime, "\n", Sys.time(), "\n"))


