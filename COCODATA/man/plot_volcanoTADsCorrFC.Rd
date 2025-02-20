% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_func.R
\name{plot_volcanoTADsCorrFC}
\alias{plot_volcanoTADsCorrFC}
\title{Volcano plot TAD meanLogFC and meanCorr}
\usage{
plot_volcanoTADsCorrFC(meanCorr, meanFC, comb_pval, padjusted = FALSE,
  tads_to_annot = NULL, x_lab = "TAD mean LogFC",
  y_lab = "TAD mean intraCorr", plotTit = NULL, subTit = NULL,
  fcUpDownThresh = 0, signifThresh = 0.05, verySignifThresh = 0.01,
  strongUp = "#EE0011FF", lightUp = "#EE001132",
  strongDown = "#0C5BB0FF", lightDown = "#0C5BB032")
}
\arguments{
\item{meanCorr}{Vector of intra-TAD mean correlations (vector names should be TAD IDs).}

\item{meanFC}{Vector of TAD meanLogFC (vector names should be TAD IDs).}

\item{comb_pval}{Vector of TAD p-values (vector names should be TAD IDs).}

\item{padjusted}{If the p-values are already adjusted (BH method).}

\item{tads_to_annot}{(optional) TAD that should be annotated with a label.}

\item{...}{Other parameters for fine-tuning the plot.}
}
\value{
The plot (ggplot).
}
\description{
Function that returns a volcano-like plot of TAD meanLogFC and intra-TAD meanCorr.
}
