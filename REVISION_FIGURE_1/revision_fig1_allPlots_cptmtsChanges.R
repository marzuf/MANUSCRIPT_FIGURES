# Rscript revision_fig1_allPlots_cptmtsChanges.R

outFolder <- file.path("REVISION_FIG1_ALLPLOTS_CPTMTSCHANGES")
dir.create(outFolder, recursive = TRUE)

runFolder <- "../../v2_Yuanlong_Cancer_HiC_data_TAD_DA"

tadSignifThresh <- 0.01
corrPurityQtThresh <- 0.05

cptmtQtThresh <- 0.75 

nHistBins <- 100

qtcol <- "red"

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


myHeightGG_bp <- 6
myWidthGG_bp <- 5.5

myHeightGG_dist <- 5.5
myWidthGG_dist <- 6.5


matching_dt <- get(load(file.path("PREP_REVISIONMATCHING_DATA/matching_withRank_dt.Rdata")))
### go ahead with notPF data !!!
matching_dt <- matching_dt[matching_dt$ref_region_ID %in% tokeep_tads & 
                             matching_dt$matching_region_ID %in% tokeep_tads 
                             ,] 
stopifnot(nrow(matching_dt)>0)

init_matching_dt <- matching_dt

#*********************************************************************************************************
################### PLOT 1: barplot cptmt changes (2 cptmt - all changes)
#*********************************************************************************************************
matching_dt <- init_matching_dt
suffix <- ""
cptmt_var_lab <- "CD 2-compartment"
cptmt_levels <- c("A->A", "B->B", "A->B", "B->A") 
bplot <- "all"

my_cols_palette <- "nejm"
legTitle <- ""

all_cmps <- unique(file.path(matching_dt$matching_hicds, matching_dt$matching_exprds,
                             matching_dt$ref_hicds, matching_dt$ref_exprds))

mySub <- paste0("# DS comparisons = ", length(all_cmps), "; # CDs = ", nrow(matching_dt), 
                " (signif.: ", sum(matching_dt$adjPval <= tadSignifThresh), ")")

agg_dt <- aggregate(as.formula(paste0("ref_region_ID~norm2tumor_cptmtChange", suffix, " + ref_tadSignif")),
                    data=matching_dt, FUN=length)
colnames(agg_dt)[colnames(agg_dt)=="ref_region_ID"] <- "nTADs"
tmp_dt <- aggregate(ref_region_ID~ref_tadSignif, data=matching_dt, FUN=length)
nTot <- setNames(tmp_dt$ref_region_ID, tmp_dt$ref_tadSignif)
agg_dt$ratioTADs <- agg_dt$nTADs/nTot[paste0(agg_dt$ref_tadSignif)]
plotTit <- paste0("norm vs. tumor ratio CDs by 2-cptmt change")
agg_dt$norm2tumor_cptmtChange <- factor(agg_dt[,paste0("norm2tumor_cptmtChange", suffix)], levels=cptmt_levels)
stopifnot(!is.na(agg_dt[,paste0("norm2tumor_cptmtChange", suffix)]))
ggbar_p <-  ggplot(agg_dt, aes_string(y="ratioTADs",
                      x="ref_tadSignif", 
                      fill=paste0("norm2tumor_cptmtChange", suffix))) +
  ggtitle(plotTit, subtitle=mySub)+
  geom_bar(position = position_stack(), stat="identity") +
  labs(x="" , y ="ratio of CDs", color=paste0(legTitle),fill=paste0(legTitle)) +
  eval(parse(text = paste0("scale_fill_", my_cols_palette, "()")))+ #
  mytheme + title_theme+
    theme(
      plot.title = element_text(size=14, face="bold", hjust=0.5),
      plot.subtitle = element_text(size=12, face="italic", hjust=0.5),
      axis.line=element_line(),
    axis.text.y = element_text(size=16, color="black"),
    axis.text.x = element_text(size=16, color="black"),
    axis.title.y = element_text(size=20, color="black")
  )

  #   axis.text.x = element_text(hjust=1, vjust=0.5,size=10,angle=90)
  # )
outFile <- file.path(outFolder,paste0("norm_vs_tumor_cptmt", suffix, "_change_signif_notSignif", "_", bplot, "_barplot.", plotType))
ggsave(ggbar_p, filename = outFile, height=myHeightGG_bp, width=myWidthGG_bp)
cat(paste0("... written: ", outFile, "\n"))



#*********************************************************************************************************
################### PLOT 2: barplot cptmt changes (8 cptmt - B->B)
#*********************************************************************************************************

matching_dt <- init_matching_dt[init_matching_dt$norm2tumor_cptmtChange == "B->B",]
stopifnot(nrow(matching_dt) > 0)
suffix <- "_eight"
my_cols_palette <- "simpsons"


matching_dt[,paste0("norm2tumor_cptmtChange", suffix)][matching_dt[,paste0("tumorCptmt", suffix)] == matching_dt[,paste0("normCptmt", suffix)]] <- "same"
bplot <- "grouped"

legTitle <- ""

all_cmps <- unique(file.path(matching_dt$matching_hicds, matching_dt$matching_exprds,
                             matching_dt$ref_hicds, matching_dt$ref_exprds))

mySub <- paste0("# DS comparisons = ", length(all_cmps), "; # TADs = ", nrow(matching_dt), 
                " (signif.: ", sum(matching_dt$adjPval <= tadSignifThresh), ")")

agg_dt <- aggregate(as.formula(paste0("ref_region_ID~norm2tumor_cptmtChange", suffix, " + ref_tadSignif")),
                    data=matching_dt, FUN=length)
colnames(agg_dt)[colnames(agg_dt)=="ref_region_ID"] <- "nTADs"
tmp_dt <- aggregate(ref_region_ID~ref_tadSignif, data=matching_dt, FUN=length)
nTot <- setNames(tmp_dt$ref_region_ID, tmp_dt$ref_tadSignif)
agg_dt$ratioTADs <- agg_dt$nTADs/nTot[paste0(agg_dt$ref_tadSignif)]
plotTit <- paste0("norm vs. tumor ratio B->B CDs by 8-cptmt change")
agg_dt$norm2tumor_cptmtChange <- factor(agg_dt[,paste0("norm2tumor_cptmtChange", suffix)], levels=rev(cptmt_levels))
stopifnot(!is.na(agg_dt[,paste0("norm2tumor_cptmtChange", suffix)]))
ggbar_p <-  ggplot(agg_dt, aes_string(y="ratioTADs",
                      x="ref_tadSignif",
                      fill=paste0("norm2tumor_cptmtChange", suffix))) +
  geom_bar(position = position_stack(), stat="identity") +
  ggtitle(plotTit, subtitle=mySub)+
  labs(x="" , y ="ratio of CDs", color=paste0(legTitle),fill=paste0(legTitle)) + 
  eval(parse(text = paste0("scale_fill_", my_cols_palette, "()")))+ #
  mytheme +
   title_theme+
  theme(
    axis.line=element_line(),
    plot.title = element_text(size=14, face="bold", hjust=0.5),
    plot.subtitle = element_text(size=12, face="italic", hjust=0.5),
    axis.text.y = element_text(size=16, color="black"),
    axis.text.x = element_text(size=16, color="black"),
    axis.title.y = element_text(size=20, color="black")
  )

outFile <- file.path(outFolder,paste0("norm_vs_tumor_cptmt", suffix, "_change_signif_notSignif", "_", bplot, "_barplot.", plotType))
ggsave(ggbar_p, filename = outFile, height=myHeightGG_bp, width=myWidthGG_bp)
cat(paste0("... written: ", outFile, "\n"))


#*********************************************************************************************************
################### PLOT 3: distribution of  abs rank diff
#*********************************************************************************************************

matching_dt <- init_matching_dt

cptmt_type <- "eight"
legTitle <- paste0("sameCptmt_", cptmt_type, "\n(same=",sum(matching_dt[,paste0("sameCptmt_", cptmt_type)]) , ")")
legTitle <- ""

plotTit <- paste0("CD rank diff. and 8-cptmt change")

all_cmps <- unique(file.path(matching_dt$matching_hicds, matching_dt$matching_exprds,
                             matching_dt$ref_hicds, matching_dt$ref_exprds))

mySub <- paste0("# DS comparisons = ", length(all_cmps), "; # TADs = ", nrow(matching_dt), 
                " (signif.: ", sum(matching_dt$adjPval <= tadSignifThresh), ")")

matching_dt$absRankDiff <- abs(matching_dt$rankDiff)

matching_dt$plot_lab <- ifelse(matching_dt[,paste0("sameCptmt_", cptmt_type)], paste0("same cptmt\n(", sum(matching_dt[,paste0("sameCptmt_", cptmt_type)]), ")"), 
                               paste0("diff. cptmt\n(", sum(!matching_dt[,paste0("sameCptmt_", cptmt_type)]), ")"))

my_cols <- setNames(pal_jama()(3)[c(2, 3)], sort(unique(matching_dt$plot_lab)))


p3 <- gghistogram(matching_dt,
                  x = paste0("absRankDiff"),
                  # y = "..count..",
                  y = "..density..",
                  bins=100,
                  # combine = TRUE,                  # Combine the 3 plots
                  xlab = paste0("Rank diff. with matched TAD"),
                  # ylab = "# of TADs",
                  ylab="Density",
                  # add = "median",                  # Add median line.
                  rug = FALSE,                      # Add marginal rug
                  color = "plot_lab", 
                  fill = "plot_lab", #paste0("sameCptmt_", cptmt_type),
                  palette = "jco"
) +
  ggtitle(plotTit, subtitle = mySub)+
  scale_color_manual(values=my_cols)+
  scale_fill_manual(values=my_cols)  +
  labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
  guides(color=FALSE)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  mytheme

outFile <- file.path(outFolder, paste0("tad_absRankDiff_same_diff_cptmt_", cptmt_type, "_density.", plotType))
ggsave(p3, file=outFile, height=myHeightGG_dist, width=myWidthGG_dist)
cat(paste0("... written: ", outFile, "\n"))

#*********************************************************************************************************
################### PLOT 4: enrichment signif in rank diff 
#*********************************************************************************************************

thresh1 <- quantile(abs(matching_dt$rankDiff)[!matching_dt$sameCptmt_eight], probs=cptmtQtThresh)
  
plot_dt=matching_dt

nsignif <- sum(plot_dt$ref_tadSignif == "signif.")
nnonsignif <- sum(plot_dt$ref_tadSignif == "not signif.")

col_var <- "absRankDiff"

curr_thresh <- thresh1
curr_thresh_name <- cptmtQtThresh
threshLab <-   paste0(cptmtQtThresh, "-qt for cptmt-changing TADs (", round(curr_thresh, 4), ")") 


print(paste0("curr_thresh=", round(curr_thresh, 4), "\n"))
print(paste0("curr_thresh_name=", round(curr_thresh_name, 4), "\n"))

tmp_ctg_dt <- matching_dt
tmp_ctg_dt$aboveQt <- tmp_ctg_dt$absRankDiff >= curr_thresh
ctg_mat <- table(tmp_ctg_dt$aboveQt, tmp_ctg_dt$ref_tadSignif)
names(dimnames(ctg_mat)) <- c("GeQt", "TAD signif." )
print(ctg_mat)


obs_geQt_signif <- ctg_mat["TRUE", "signif."]
obs_geQt_nonsignif <- ctg_mat["TRUE", "not signif."]


all_cmps <- unique(file.path(matching_dt$matching_hicds, matching_dt$matching_exprds,
                             matching_dt$ref_hicds, matching_dt$ref_exprds))

mySub <- paste0("# DS comparisons = ", length(all_cmps), "; # TADs = ", nrow(matching_dt), 
                " (signif.: ", sum(matching_dt$adjPval <= tadSignifThresh), ")")

mylineTxt <- paste(
  threshLab,
  paste0(round(100*obs_geQt_signif/sum(ctg_mat[, "signif."]) , 2), "% of signif. TADs"),
  paste0(round(100*obs_geQt_nonsignif/sum(ctg_mat[, "not signif."]) , 2), "% of not signif. TADs"),
  paste0("Fisher test p-val = ", formatC(fisher.test(ctg_mat)$p.value, format = "e", digits = 4)),
  # paste0("emp. p-val signif. (n=", nRandom_stats,") = ", formatC(empPval_signif, format = "e", digits = 4)),
  # paste0("emp. p-val nonsignif. (n=", nRandom_stats,") = ", formatC(empPval_nonsignif, format = "e", digits = 4)),
  sep="\n")

plotTit <- paste0("CD abs. rank diff. and signif. (not PF)")

p3 <- ggdensity(matching_dt,
                  x = paste0("absRankDiff"),
                  y = "..density..",
                  # y = "..count..",
                  bins=nHistBins,
                  # ylab ="# TADs",
                  ylab ="Density",
                  # combine = TRUE,                  # Combine the 3 plots
                  xlab = paste0("Abs. rank diff. with matched CD"),
                  # add = "median",                  # Add median line.
                  rug = FALSE,                      # Add marginal rug
                  color = "ref_tadSignif",
                  fill = "ref_tadSignif",
                  palette = "jco"
) +
  ggtitle(plotTit, subtitle = mySub)+
  scale_color_manual(values=mycolsSignif)+
  scale_fill_manual(values=mycolsSignif)  +
  guides(color=FALSE)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(color=paste0(legTitle),fill=paste0(legTitle), y="Density") +
  mytheme +
  geom_vline(xintercept=curr_thresh, linetype=2, col=qtcol)


txtY <- ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range[2]

p3 <- p3 +   annotate("text", 
                      label=mylineTxt, 
                      x=curr_thresh, 
                      y =txtY, color = qtcol, hjust = 0, vjust=1)


outFile <- file.path(outFolder, paste0("tad_absRankDiff_signif_notsignif_thresh", 
                                       gsub("\\.", "", as.character(curr_thresh_name)), "_density.", plotType))
ggsave(p3, file=outFile, height=myHeightGG_dist, width=myWidthGG_dist)
cat(paste0("... written: ", outFile, "\n"))



