# Rscript revision_fig1_allPlots.R

outFolder <- file.path("REVISION_FIG1_ALLPLOTS_CPTMTS")
dir.create(outFolder, recursive = TRUE)

runFolder <- "../../v2_Yuanlong_Cancer_HiC_data_TAD_DA"

tadSignifThresh <- 0.01
corrPurityQtThresh <- 0.05

revFig_dt <- get(load(file.path( "PREP_REVISIONFIG1_DATA/revision_fig1_cptmtAnnot_with_corr_purity.Rdata")))
revFig_dt$regionID <- file.path(revFig_dt$hicds, revFig_dt$exprds, revFig_dt$region)

## go ahead with only the notPF data !!!
revFig_dt$signif <- revFig_dt$adjPvalComb <= tadSignifThresh
purityCorrThresh <- as.numeric(quantile(revFig_dt$purityCorr[!revFig_dt$signif], probs = corrPurityQtThresh ))
tokeep_tads <- revFig_dt$regionID[revFig_dt$purityCorr > purityCorrThresh]
revFig_dt <- revFig_dt[revFig_dt$purityCorr > purityCorrThresh,]

require(ggplot2)
require(ggrepel)
require(ggsci)
require(ggpubr)
source(file.path(runFolder, "revision_settings.R"))
source(file.path("..", "revision_settings.R"))

plotType <- "svg"
myHeightGG_pie <- 5
myWidthGG_pie <- 5.5

myHeightGG_bp8 <- 5
myWidthGG_bp8 <- 6

myWidthGG_boxCptmt <- 7
myHeightGG_boxCptmt <- 5.5

#*********************************************************************************************************
################### PLOT 1: PIE CHART SIGNIF/NOT SIGNIF - pie chart binary
#*********************************************************************************************************

cptmt_var <- "binaryCptmtLab"
cptmt_var_lab <- "CD 2-compartment"

tmp_dt <- aggregate(as.formula(paste0("regionID ~ ", cptmt_var)), data=revFig_dt, FUN=length)
colnames(tmp_dt)[colnames(tmp_dt) == "regionID"] <- "nTADs"
tmp_dt$ratioTADs <- tmp_dt$nTADs/nrow(revFig_dt)
tmp_dt$ratioTADs_lab <- paste0(round(tmp_dt$ratioTADs*100, 2), "%")

myTit <- paste0("Dist. of all CDs across ", cptmt_var_lab)
mysub <- paste0("purity-filtered: # DS = ", length(unique(file.path(revFig_dt$hicds, revFig_dt$exprds))) , "; # CDs = ", nrow(revFig_dt))

p_dist <- ggplot(tmp_dt, aes_string(x="1", y="ratioTADs", fill=cptmt_var)) +
  geom_col() +
  scale_fill_manual(values=mycolsBinary)+
  geom_text_repel(aes(label = ratioTADs_lab), size=6, position = position_stack(vjust = 0.5))+
  coord_polar(theta = "y") + 
  ggtitle(myTit, subtitle=mysub) +
  theme_void() +    
  labs(fill="") +
  blank_theme +
  title_theme


outFile <- file.path(outFolder, paste0( "distAllTADs_by_", cptmt_var, "_pie.", plotType))
ggsave(p_dist, file=outFile, height=myHeightGG_pie, width=myWidthGG_pie)
cat(paste0("... written: ", outFile, "\n"))

signif_dt <- revFig_dt[revFig_dt$adjPvalComb <= tadSignifThresh,]
agg_dt <- aggregate(as.formula(paste0("regionID ~ ", cptmt_var)), data=signif_dt, FUN=length)
colnames(agg_dt)[colnames(agg_dt) == "regionID"] <- "nTADs"
agg_dt$ratioTADs <- agg_dt$nTADs/nrow(signif_dt)

agg_dt$ratioTADs_lab <- paste0(round(agg_dt$ratioTADs*100, 2), "%")

myTit <- paste0("Dist. signif. CDs across ", cptmt_var_lab)
mysub <- paste0("purity-filtered: # DS = ", length(unique(file.path(signif_dt$hicds, signif_dt$exprds))) , "; # CDs = ", nrow(signif_dt))

p_dist <- ggplot(agg_dt, aes_string(x="1", y="ratioTADs", fill=cptmt_var)) +
  geom_col() +
  scale_fill_manual(values=mycolsBinary)+
  geom_text_repel(aes(label = ratioTADs_lab), size=6, position = position_stack(vjust = 0.5))+
  coord_polar(theta = "y") + 
  ggtitle(myTit, subtitle=mysub) +
  theme_void() +    
  labs(fill="") +
  blank_theme +
  theme(
    plot.title = element_text(size=14, face="bold", hjust=0.5),
    plot.subtitle = element_text(size=12, face="italic", hjust=0.5)
  )

outFile <- file.path(outFolder, paste0( "distSignifTADs_by_", cptmt_var, "_pie.", plotType))
ggsave(p_dist, file=outFile, height=myHeightGG_pie, width=myWidthGG_pie)
cat(paste0("... written: ", outFile, "\n"))


#*********************************************************************************************************
################### PLOT 2: BARPLOT CHART SIGNIF/NOT SIGNIF - barplot enrichment eight lab
#*********************************************************************************************************

cptmt_var <- "eightCptmtLab"
cptmt_var_lab <- "CD 8-compartment"

all_agg_dt <- aggregate(as.formula(paste0("regionID ~ ", cptmt_var)), data=revFig_dt, FUN=length)
colnames(all_agg_dt)[colnames(all_agg_dt) == "regionID"] <- "nTADs"
all_agg_dt$ratioTADs <- all_agg_dt$nTADs/nrow(revFig_dt)
all_agg_dt$ratioTADs_lab <- paste0(round(all_agg_dt$ratioTADs*100, 2), "%")


signif_dt <- revFig_dt[revFig_dt$adjPvalComb <= tadSignifThresh,]
signif_agg_dt <- aggregate(as.formula(paste0("regionID ~ ", cptmt_var)), data=signif_dt, FUN=length)
colnames(signif_agg_dt)[colnames(signif_agg_dt) == "regionID"] <- "nTADs"
signif_agg_dt$ratioTADs <- signif_agg_dt$nTADs/nrow(signif_dt)

signif_agg_dt$ratioTADs_lab <- paste0(round(signif_agg_dt$ratioTADs*100, 2), "%")

bp_all_dt <- all_agg_dt
bp_signif_dt <- signif_agg_dt
both_dt <- merge(bp_all_dt, bp_signif_dt, by=c(cptmt_var), suffixes = c("_all", "_signif"), all=TRUE)

both_dt$signif_enrich <-  both_dt$ratioTADs_signif/both_dt$ratioTADs_all
stopifnot(!is.na(both_dt$signif_enrich))

both_dt$signif_enrich_hk <-  both_dt$signif_enrich-1

shift_trans = function(d = 0) {
  scales::trans_new("shift", transform = function(x) x - d, inverse = function(x) x + d)
}

myTit <- paste0("Signif./all CDs dist. across ", cptmt_var_lab)


enrich_p <- ggbarplot(data=both_dt,x=cptmt_var, y = "signif_enrich", fill = cptmt_var,
                      xlab="", ylab ="signif./all ratio of CDs") + 
  geom_hline(yintercept = 1) +
  scale_y_continuous(trans = shift_trans(1)) +
  ggtitle(myTit, subtitle=mysub)+
  eval(parse(text = paste0("scale_fill_", eightCptmtPalette, "()")))+ #
  labs(fill="")+
  title_theme
# blank_theme

outFile <- file.path(outFolder, paste0( "ratio_distSignifOverAllTADs_by_", cptmt_var, "_barplot.", plotType))
ggsave(enrich_p, file=outFile, height=myHeightGG_bp8, width=myWidthGG_bp8)
cat(paste0("... written: ", outFile, "\n"))

#*********************************************************************************************************
################### PLOT 3: distribution signif not signif distribution across compartments
#*********************************************************************************************************

# !!! need to use aes and not aes_string -> see: https://github.com/tidyverse/ggplot2/issues/2896
# geom_split_violin from https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
# see also interesting with quantile  https://stackoverflow.com/questions/47651868/split-violin-plot-with-ggplot2-with-quantiles/47652563#47652563

cptmt_var <- "eightCptmtLab"
cptmt_var_lab <- "CD 8-compartment"
plot_var <- "meanCorr"

plot_var_labs <- c(meanCorr = "CD mean corr.",
                   aggLog10Expr = "CD mean expression [log10]",
                   purityCorr = "CD mean corr. with purity"
                   )

legTitle <- ""

myTit <- paste0("")
mysub <- paste0("purity-filtered: # DS = ", length(unique(file.path(revFig_dt$hicds, revFig_dt$exprds))) , 
                "; # CDs = ", nrow(revFig_dt), " (signif.: ", sum(revFig_dt$signif), ")")


for(plot_var in names(plot_var_labs)) {
  revFig_dt$plot_var <- revFig_dt[,plot_var]
  revFig_dt$cptmt_var <- revFig_dt[,cptmt_var]
  
  revFig_dt$signif_lab <- ifelse(revFig_dt$signif, "signif.", "not signif.")
  
  p_uniqueBox1 <- ggplot(revFig_dt, aes(x = cptmt_var, y = plot_var)) + 
    # geom_split_violin() +
    geom_violin(fill="azure3") +
    scale_fill_manual(values=mycolsSignif)+
    geom_boxplot(fill="white", color="dimgray", notch=TRUE, width=0.1) + 
    ggtitle(myTit, subtitle=mysub)+
    labs(x="", y=plot_var_labs[plot_var], fill=legTitle, col=legTitle)+ 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    mytheme +
    title_theme + 
    leg_theme +
    theme(axis.line=element_line(), 
          axis.title=element_text(color="black", size=14),
          axis.text=element_text(color="black", size=12)
    )  
  outFile <- file.path(outFolder,paste0(plot_var, "_byCptmt_", cptmt_var, "_1dist_with1box_violinplot.", plotType))
  ggsave(p_uniqueBox1, filename = outFile, height=myHeightGG_boxCptmt, width=myWidthGG_boxCptmt)
  cat(paste0("... written: ", outFile, "\n"))
  next
  
  p_uniqueBox <- ggplot(revFig_dt, aes(x = cptmt_var, y = plot_var, fill=signif_lab)) + 
    geom_split_violin() +
    scale_fill_manual(values=mycolsSignif)+
    geom_boxplot(fill="white", color="dimgray", notch=TRUE, width=0.1) + 
    ggtitle(myTit, subtitle=mysub)+
    labs(x="", y=plot_var_labs[plot_var], fill=legTitle, col=legTitle)+ 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    mytheme +
    title_theme + 
    leg_theme +
    theme(axis.line=element_line(), 
          axis.title=element_text(color="black", size=14),
          axis.text=element_text(color="black", size=12)
    )  
  outFile <- file.path(outFolder,paste0(plot_var, "_byCptmt_", cptmt_var, "_2dist_with1box_violinplot.", plotType))
  ggsave(p_uniqueBox, filename = outFile, height=myHeightGG_boxCptmt, width=myWidthGG_boxCptmt)
  cat(paste0("... written: ", outFile, "\n"))
  
  p_twoBox <- ggplot(revFig_dt, aes(x = cptmt_var, y = plot_var, fill=signif_lab)) + 
    geom_split_violin() +
    scale_fill_manual(values=mycolsSignif)+
    geom_boxplot(width=0.1, color="dimgray") +
    ggtitle(myTit, subtitle=mysub)+
    labs(x="", y=plot_var_labs[plot_var], fill=legTitle, col=legTitle)+ 
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    mytheme +
    title_theme + 
    leg_theme + 
    theme(axis.line=element_line(), 
          axis.title=element_text(color="black", size=14),
          axis.text=element_text(color="black", size=12)
          )
  
  outFile <- file.path(outFolder,paste0(plot_var, "_byCptmt_", cptmt_var, "_2dist_with2box_violinplot.", plotType))
  ggsave(p_twoBox, filename = outFile, height=myHeightGG_boxCptmt, width=myWidthGG_boxCptmt)
  cat(paste0("... written: ", outFile, "\n"))
  
}


















################### ################### ################### ################### ################### 
################### ################### ################### ################### ################### TRASH
################### ################### ################### ################### ################### 

inFolder <- file.path(runFolder, "REVISION_CHANGES_CPTMTLABELS_ALLDS")
outFile <- file.path(inFolder, "tad2cptmt_dt.Rdata")
all_tad2cptmt_dt <- get(load(outFile))
## go ahead with only the notPF data !!!
stopifnot(tokeep_tads %in% all_tad2cptmt_dt$region_ID)
all_tad2cptmt_dt <- all_tad2cptmt_dt[all_tad2cptmt_dt$region_ID %in% tokeep_tads,]

stopifnot(all_tad2cptmt_dt$tad_binaryCptmtLab == all_tad2cptmt_dt$start_binaryCptmtLab)

# TO CHECK THE DATA INPUT
source("../full_dataset_names.R")
# hicds_names <- setNames(names(hicds_names), hicds_names)
# exprds_names <- setNames(names(exprds_names), exprds_names)
# stopifnot(revFig_dt$hicds %in% names(hicds_names))
# revFig_dt$hicds_lab <- hicds_names[paste0(revFig_dt$hicds)]
# stopifnot(revFig_dt$exprds %in% names(exprds_names))
# revFig_dt$exprds_lab <- exprds_names[paste0(revFig_dt$exprds)]
# revFig_dt$regionID <- file.path(revFig_dt$hicds_lab, revFig_dt$exprds_lab, revFig_dt$region)
stopifnot(all_tad2cptmt_dt$hicds %in% names(hicds_names))
all_tad2cptmt_dt$hicds_lab <- hicds_names[paste0(all_tad2cptmt_dt$hicds)]
stopifnot(all_tad2cptmt_dt$exprds %in% names(exprds_names))
all_tad2cptmt_dt$exprds_lab <- exprds_names[paste0(all_tad2cptmt_dt$exprds)]
all_tad2cptmt_dt$regionID <- file.path(all_tad2cptmt_dt$hicds_lab, all_tad2cptmt_dt$exprds_lab, all_tad2cptmt_dt$region)
all_tad2cptmt_dt$regionID <- file.path(all_tad2cptmt_dt$hicds, all_tad2cptmt_dt$exprds, all_tad2cptmt_dt$region)

tmp_labs <- setNames(revFig_dt$eightCptmtLab, revFig_dt$regionID)
tmp_pvals <- setNames(revFig_dt$adjPvalComb, revFig_dt$regionID)

tmp_labs2 <- setNames(all_tad2cptmt_dt$tad_eightCptmtLab, all_tad2cptmt_dt$regionID)
stopifnot(names(tmp_labs) %in% names(tmp_labs2))
tmp_labs2 <- tmp_labs2[paste0(names(tmp_labs))]
stopifnot(!is.na(tmp_labs2))
stopifnot(tmp_labs == tmp_labs2)
stopifnot(names(tmp_labs) == names(tmp_labs2))

tmp_pvals2 <- setNames(all_tad2cptmt_dt$adjPvalComb, all_tad2cptmt_dt$regionID)
stopifnot(names(tmp_pvals) %in% names(tmp_pvals2))
tmp_pvals2 <- tmp_pvals2[paste0(names(tmp_pvals))]
stopifnot(!is.na(tmp_pvals2))
stopifnot(tmp_pvals == tmp_pvals2)
stopifnot(names(tmp_pvals) == names(tmp_pvals2))


