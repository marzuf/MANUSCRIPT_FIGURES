tadSignifThresh <- 0.01
corrPurityQtThresh <- 0.05

revFig_dt <- get(load(file.path( "../REVISION_FIGURE_1/PREP_REVISIONFIG1_DATA/revision_fig1_cptmtAnnot_with_corr_purity.Rdata")))
revFig_dt$regionID <- file.path(revFig_dt$hicds, revFig_dt$exprds, revFig_dt$region)

## go ahead with only the notPF data !!!
revFig_dt$signif <- revFig_dt$adjPvalComb <= tadSignifThresh
purityCorrThresh <- as.numeric(quantile(revFig_dt$purityCorr[!revFig_dt$signif], probs = corrPurityQtThresh ))
tokeep_tads <- revFig_dt$regionID[revFig_dt$purityCorr > purityCorrThresh]
filt_revFig_dt <- revFig_dt[revFig_dt$purityCorr > purityCorrThresh,]


stopifnot("GSE118514_RWPE1_40kb/TCGAprad_norm_prad/chr12_TAD194" %in% filt_revFig_dt$regionID)
stopifnot("GSE118514_RWPE1_40kb/TCGAprad_norm_prad/chr7_TAD424"%in% filt_revFig_dt$regionID)
stopifnot("GSE118514_RWPE1_40kb/TCGAprad_norm_prad/chr17_TAD174"%in% filt_revFig_dt$regionID)
stopifnot("GSE118514_RWPE1_40kb/TCGAprad_norm_prad/chr17_TAD147"%in% filt_revFig_dt$regionID)


stopifnot("GSE118514_22Rv1_40kb/TCGAprad_norm_prad/chr12_TAD196"%in% filt_revFig_dt$regionID)
stopifnot("GSE118514_22Rv1_40kb/TCGAprad_norm_prad/chr1_TAD460"%in% filt_revFig_dt$regionID)
stopifnot("GSE118514_22Rv1_40kb/TCGAprad_norm_prad/chr17_TAD135"%in% filt_revFig_dt$regionID)
stopifnot("GSE118514_22Rv1_40kb/TCGAprad_norm_prad/chr17_TAD269"%in% filt_revFig_dt$regionID)
