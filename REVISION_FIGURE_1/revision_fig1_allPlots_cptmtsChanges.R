# Rscript revision_fig1_allPlots_cptmtsChanges.R

outFolder <- file.path("REVISION_FIG1_ALLPLOTS_CPTMTSCHANGES")
dir.create(outFolder, recursive = TRUE)

runFolder <- "../../v2_Yuanlong_Cancer_HiC_data_TAD_DA"

tadSignifThresh <- 0.01
corrPurityQtThresh <- 0.05

cptmtQtThresh <- 0.9

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
require(scatterpie)
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

cptmt_var_lab <- "2-cptmt change"

tmp_dt <- matching_dt
tmp_dt$agg_col <- tmp_dt[,paste0("norm2tumor_cptmtChange", suffix)]

nonsignif_dt <- tmp_dt[tmp_dt$ref_region_pval > tadSignifThresh,]
nonsignif_agg_dt <- aggregate(as.formula(paste0("ref_region_ID ~ ", "agg_col")), data=nonsignif_dt, FUN=length)
colnames(nonsignif_agg_dt)[colnames(nonsignif_agg_dt) == "ref_region_ID"] <- "nTADs"
nonsignif_agg_dt$ratioTADs <- nonsignif_agg_dt$nTADs/nrow(nonsignif_dt)
nonsignif_agg_dt$ratioTADs_lab <- paste0(round(nonsignif_agg_dt$ratioTADs*100, 2), "%")

signif_dt <- tmp_dt[tmp_dt$ref_region_pval <= tadSignifThresh,]
signif_agg_dt <- aggregate(as.formula(paste0("ref_region_ID ~ ", "agg_col")), data=signif_dt, FUN=length)
colnames(signif_agg_dt)[colnames(signif_agg_dt) == "ref_region_ID"] <- "nTADs"
signif_agg_dt$ratioTADs <- signif_agg_dt$nTADs/nrow(signif_dt)
signif_agg_dt$ratioTADs_lab <- paste0(round(signif_agg_dt$ratioTADs*100, 2), "%")

bp_nonsignif_dt <- nonsignif_agg_dt
bp_signif_dt <- signif_agg_dt
both_dt <- merge(bp_nonsignif_dt, bp_signif_dt, by=c("agg_col"), suffixes = c("_nonsignif", "_signif"), all=TRUE)

both_dt$signif_enrich <-  both_dt$ratioTADs_signif/both_dt$ratioTADs_nonsignif
stopifnot(!is.na(both_dt$signif_enrich))

both_dt$signif_enrich_hk <-  both_dt$signif_enrich-1

shift_trans = function(d = 0) {
  scales::trans_new("shift", transform = function(x) x - d, inverse = function(x) x + d)
}

myTit <- paste0("Signif./non signif. CDs dist. across ", cptmt_var_lab)

mysub <- paste0("# DS comparisons = ", length(all_cmps), "; # CDs = ", nrow(matching_dt), 
                " (signif.: ", sum(matching_dt$adjPval <= tadSignifThresh), ")")

both_dt$agg_col <- factor(both_dt$agg_col, levels=cptmt_levels)
stopifnot(!is.na(both_dt$agg_col))

enrich_p <- ggbarplot(data=both_dt,x="agg_col", y = "signif_enrich", fill = "agg_col",
                      xlab="", ylab ="signif./non signif. ratio of CDs") + 
  geom_hline(yintercept = 1) +
  scale_y_continuous(trans = shift_trans(1)) +
  ggtitle(myTit, subtitle=mysub)+
  eval(parse(text = paste0("scale_fill_", eightCptmtPalette, "()")))+ #
  labs(fill="")+
  title_theme +
  theme(
    # plot.title = element_text(size=14, face="bold", hjust=0.5),
    # plot.subtitle = element_text(size=12, face="italic", hjust=0.5),
    axis.line=element_line(),
    axis.text.y = element_text(size=16, color="black"),
    axis.text.x = element_text(size=16, color="black"),
    axis.title.y = element_text(size=20, color="black")
  )
# blank_theme

outFile <- file.path(outFolder, paste0( "ratio_distSignifOverAllTADs_by_", "cptmt_change", "_barplot.", plotType))
ggsave(enrich_p, file=outFile, height=myHeightGG_bp, width=myWidthGG_bp)
cat(paste0("... written: ", outFile, "\n"))

save(both_dt, file=file.path(outFolder, "both_dt.Rdata"), version=2)
save(nonsignif_agg_dt, file=file.path(outFolder, "nonsignif_agg_dt.Rdata"), version=2)
save(signif_agg_dt, file=file.path(outFolder, "signif_agg_dt.Rdata"), version=2)


#*********************************************************************************************************
################### PLOT 2: barplot cptmt changes with pie chart on the top same/diff 8 cptmt
#*********************************************************************************************************

enrich_p <- ggbarplot(data=both_dt,x="agg_col", y = "signif_enrich", fill = "agg_col",
                      xlab="", ylab ="signif./non signif. ratio of CDs") + 
  geom_hline(yintercept = 1) +
  scale_y_continuous(trans = shift_trans(1)) +
  ggtitle(myTit, subtitle=mysub)+
  eval(parse(text = paste0("scale_fill_", eightCptmtPalette, "()")))+ #
  labs(fill="")+
  title_theme +
  theme(
    # plot.title = element_text(size=14, face="bold", hjust=0.5),
    # plot.subtitle = element_text(size=12, face="italic", hjust=0.5),
    axis.line=element_line(),
    axis.text.y = element_text(size=16, color="black"),
    axis.text.x = element_text(size=16, color="black"),
    axis.title.y = element_text(size=20, color="black")
  )
# blank_theme

sub_matching_dt <- init_matching_dt[init_matching_dt$norm2tumor_cptmtChange == "B->B" |
                                      init_matching_dt$norm2tumor_cptmtChange == "A->A",]
stopifnot(nrow(sub_matching_dt) > 0)
suffix <- "_eight"
mycol <- paste0("norm2tumor_cptmtChange", suffix) 
sub_matching_dt[,mycol][sub_matching_dt[,paste0("tumorCptmt", suffix)] == sub_matching_dt[,paste0("normCptmt", suffix)]] <- "same 8-cptmt label"
sub_matching_dt[,mycol][sub_matching_dt[,paste0("tumorCptmt", suffix)] != sub_matching_dt[,paste0("normCptmt", suffix)]] <- "diff. 8-cptmt label"
# stopifnot(unique(sub_matching_dt[,mycol]) %in% c("same8", "diff8"))
stopifnot(unique(sub_matching_dt[,mycol]) %in% c("same 8-cptmt label", "diff. 8-cptmt label"))

agg_dt <- aggregate(as.formula(paste0("ref_region_ID", "~", mycol , " + norm2tumor_cptmtChange")), data=sub_matching_dt, FUN=length)

agg_dt$x_pos <- as.numeric(sapply(agg_dt$norm2tumor_cptmtChange, function(x)which(levels(both_dt$agg_col) == x)))
stopifnot(!is.na(agg_dt$x_pos))

both_dt$nTADs <- both_dt$nTADs_nonsignif+both_dt$nTADs_signif
agg_dt$totTADs <- as.numeric(sapply(agg_dt$norm2tumor_cptmtChange, function(x)both_dt$nTADs[which(levels(both_dt$agg_col) == x)]))
agg_dt$value <- agg_dt$ref_region_ID/agg_dt$totTADs
stopifnot(!is.na(agg_dt$x_pos))
agg_dt$y_pos_tmp <- as.numeric(sapply(agg_dt$norm2tumor_cptmtChange, function(x)both_dt$signif_enrich_hk[which(levels(both_dt$agg_col) == x)]))
stopifnot(!is.na(agg_dt$y_pos_tmp))

yoffset <- 0.15
pieradius <- 0.2
textoffset <- 0.25

stopifnot(sum(agg_dt$ref_region_ID) == nrow(init_matching_dt[init_matching_dt$norm2tumor_cptmtChange=="B->B" | 
                                                     init_matching_dt$norm2tumor_cptmtChange=="A->A",]))

agg_dt$y_pos <- ceiling(agg_dt$y_pos_tmp) 
agg_dt$y_pos[agg_dt$y_pos < 1] <- 1 + yoffset#0
agg_dt$y_pos <- agg_dt$y_pos + yoffset
# agg_dt$y_pos <- 2
agg_dt$radius <- pieradius

enrich_p1 <- enrich_p +
  geom_scatterpie(aes(x=x_pos, y=y_pos, r=radius),
                  # r=radius,
                  # pie_scale=15,
                  data=agg_dt, cols="norm2tumor_cptmtChange_eight", long_format=TRUE) 
  
annot_vect <- c(by(agg_dt, agg_dt$norm2tumor_cptmtChange, function(x) {
  tmp_dt <- x[order(x$norm2tumor_cptmtChange_eight),]
  paste0(paste0(round(tmp_dt$value*100, 2), "%"), collapse="-")
}))
annot_dt <- agg_dt[,c("norm2tumor_cptmtChange", "x_pos", "y_pos")]
annot_dt$lab <- annot_vect[agg_dt$norm2tumor_cptmtChange]
stopifnot(!is.na(annot_dt$lab))

enrich_p2 <- enrich_p1 + 
  annotate("text", x = annot_dt$x_pos, y = annot_dt$y_pos + textoffset, label=annot_dt$lab )#+
# annotate("text", x =1, y = 1, label="test" )


# save(agg_dt, file=file.path(outFolder, "agg_dt_pie.Rdata"), version=2)


outFile <- file.path(outFolder, paste0( "ratio_distSignifOverAllTADs_by_", "cptmt_change", "_barplot_withPies.", plotType))
ggsave(enrich_p2, file=outFile, height=myHeightGG_bp, width=myWidthGG_bp*1.5)
cat(paste0("... written: ", outFile, "\n"))

#*********************************************************************************************************
################### PLOT XXX  - pie chart 8 cptmt change 1) all CDs, 2) only signif CDs
#*********************************************************************************************************
save(init_matching_dt, file=file.path(outFolder, "init_matching_dt.Rdata"), version=2)
myHeightGG_pie <- 5
myWidthGG_pie <- 5.5

all_cmps <- unique(file.path(init_matching_dt$matching_hicds, init_matching_dt$matching_exprds,
                             init_matching_dt$ref_hicds, init_matching_dt$ref_exprds))
mysub <- paste0("# DS comparisons = ", length(all_cmps), "; # CDs = ", nrow(init_matching_dt), 
                " (signif.: ", sum(init_matching_dt$adjPval <= tadSignifThresh), ")")


cptmt_var <- "sameCptmt_eight"
cptmt_var_lab <- "CD 8-compartment"

mycolsSame <- c("chocolate1", "darkslateblue")

tmp_dt <- aggregate(as.formula(paste0("ref_region_ID ~ ", cptmt_var)), data=init_matching_dt, FUN=length)
colnames(tmp_dt)[colnames(tmp_dt) == "ref_region_ID"] <- "nTADs"
tmp_dt$ratioTADs <- tmp_dt$nTADs/nrow(init_matching_dt)
tmp_dt$ratioTADs_lab <- paste0(round(tmp_dt$ratioTADs*100, 2), "%")

tmp_dt$cptmt_var2 <- tmp_dt[,cptmt_var]
tmp_dt$cptmt_var2[tmp_dt[,cptmt_var]] <- "same"
tmp_dt$cptmt_var2[!tmp_dt[,cptmt_var]] <- "diff."

myTit <- paste0("Dist. of all CDs across ", cptmt_var_lab)

p_dist <- ggplot(tmp_dt, aes_string(x="1", y="ratioTADs", fill="cptmt_var2")) +
  geom_col() +
  scale_fill_manual(values=mycolsSame)+
  geom_text_repel(aes(label = ratioTADs_lab), size=6, position = position_stack(vjust = 0.5))+
  coord_polar(theta = "y") + 
  ggtitle(myTit, subtitle=mysub) +
  theme_void() +    
  labs(fill="8-cptmt label") +
  blank_theme +
  title_theme + 
  theme(legend.text=element_text(size=10),
        legend.title=element_text(face="bold"))


outFile <- file.path(outFolder, paste0( "distAllTADs_by_", cptmt_var, "_pie.", plotType))
ggsave(p_dist, file=outFile, height=myHeightGG_pie, width=myWidthGG_pie)
cat(paste0("... written: ", outFile, "\n"))


signif_dt <- init_matching_dt[init_matching_dt$ref_region_pval <= tadSignifThresh,]

mysub <- paste0("# DS comparisons = ", length(all_cmps), "; # CDs = ", nrow(signif_dt))

agg_dt <- aggregate(as.formula(paste0("ref_region_ID ~ ", cptmt_var)), data=signif_dt, FUN=length)
colnames(agg_dt)[colnames(agg_dt) == "ref_region_ID"] <- "nTADs"
agg_dt$ratioTADs <- agg_dt$nTADs/nrow(signif_dt)

agg_dt$ratioTADs_lab <- paste0(round(agg_dt$ratioTADs*100, 2), "%")
agg_dt$cptmt_var2 <- agg_dt[,cptmt_var]
agg_dt$cptmt_var2[agg_dt[,cptmt_var]] <- "same"
agg_dt$cptmt_var2[!agg_dt[,cptmt_var]] <- "diff."

myTit <- paste0("Dist. signif. CDs across ", cptmt_var_lab)

p_dist <- ggplot(agg_dt, aes_string(x="1", y="ratioTADs", fill="cptmt_var2")) +
  geom_col() +
  scale_fill_manual(values=mycolsSame)+
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

all_by8cptmt_pie <- tmp_dt
save(tmp_dt, file=file.path(outFolder, "all_by8cptmt_pie.Rdata"), version=2)
signif_by8cptmt_pie <- agg_dt
save(signif_by8cptmt_pie, file=file.path(outFolder, "signif_by8cptmt_pie.Rdata"), version=2)

 # stop("--ok\n")




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
matching_dt$absRankDiff <- abs(matching_dt$rankDiff)


# for plot 3 and 4:
col_var <- "absRankDiff"
thresh1 <- quantile(matching_dt$absRankDiff[!matching_dt$sameCptmt_eight], probs=cptmtQtThresh)
curr_thresh <- thresh1
curr_thresh_name <- cptmtQtThresh
threshLab <-   paste0(cptmtQtThresh, "-qt for cptmt-changing CDs (", round(curr_thresh, 4), ")") 

cptmt_type <- "eight"
legTitle <- paste0("sameCptmt_", cptmt_type, "\n(same=",sum(matching_dt[,paste0("sameCptmt_", cptmt_type)]) , ")")
legTitle <- ""

plotTit <- paste0("CD rank diff. and 8-cptmt change")

all_cmps <- unique(file.path(matching_dt$matching_hicds, matching_dt$matching_exprds,
                             matching_dt$ref_hicds, matching_dt$ref_exprds))

mySub <- paste0("# DS comparisons = ", length(all_cmps), "; # CDs = ", nrow(matching_dt), 
                " (signif.: ", sum(matching_dt$adjPval <= tadSignifThresh), ")")


matching_dt$plot_lab <- ifelse(matching_dt[,paste0("sameCptmt_", cptmt_type)], paste0("same cptmt\n(", sum(matching_dt[,paste0("sameCptmt_", cptmt_type)]), ")"), 
                               paste0("diff. cptmt\n(", sum(!matching_dt[,paste0("sameCptmt_", cptmt_type)]), ")"))

my_cols <- setNames(pal_jama()(3)[c(2, 3)], sort(unique(matching_dt$plot_lab)))

mylineTxt <- paste(
  threshLab)

save(matching_dt, file=file.path(outFolder, "density_matching_dt.Rdata"), version=2)

p3 <- gghistogram(matching_dt,
                  x = paste0("absRankDiff"),
                  # y = "..count..",
                  y = "..density..",
                  bins=100,
                  # combine = TRUE,                  # Combine the 3 plots
                  xlab = paste0("Abs. rank diff. with matched CD"),
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
  mytheme + 
theme(legend.text = element_text(size=12))+
  geom_vline(xintercept=curr_thresh, linetype=2, col=qtcol)


txtY <- ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range[2]

p3 <- p3 +   annotate("text", 
                      label=mylineTxt, 
                      x=curr_thresh, 
                      y =txtY, color = qtcol, hjust = 0, vjust=1)


outFile <- file.path(outFolder, paste0("tad_absRankDiff_same_diff_cptmt_", cptmt_type, "_", gsub("\\.", "",cptmtQtThresh), "qt_density.", plotType))
ggsave(p3, file=outFile, height=myHeightGG_dist, width=myWidthGG_dist)
cat(paste0("... written: ", outFile, "\n"))

#*********************************************************************************************************
################### PLOT 4: enrichment signif in rank diff 
#*********************************************************************************************************

  
plot_dt=matching_dt

nsignif <- sum(plot_dt$ref_tadSignif == "signif.")
nnonsignif <- sum(plot_dt$ref_tadSignif == "not signif.")


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

mySub <- paste0("# DS comparisons = ", length(all_cmps), "; # CDs = ", nrow(matching_dt), 
                " (signif.: ", sum(matching_dt$adjPval <= tadSignifThresh), ")")

mylineTxt <- paste(
  threshLab,
  paste0(round(100*obs_geQt_signif/sum(ctg_mat[, "signif."]) , 2), "% of signif. CDs"),
  paste0(round(100*obs_geQt_nonsignif/sum(ctg_mat[, "not signif."]) , 2), "% of not signif. CDs"),
  paste0("Fisher test p-val = ", formatC(fisher.test(ctg_mat)$p.value, format = "e", digits = 4)),
  # paste0("emp. p-val signif. (n=", nRandom_stats,") = ", formatC(empPval_signif, format = "e", digits = 4)),
  # paste0("emp. p-val nonsignif. (n=", nRandom_stats,") = ", formatC(empPval_nonsignif, format = "e", digits = 4)),
  sep="\n")

plotTit <- paste0("CD abs. rank diff. and signif.")

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



