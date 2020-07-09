## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
rm(list=ls())

if(!require(COCODATA))
  devtools::install_github("marzuf/MANUSCRIPT_FIGURES", subdir="COCODATA")
  # alternatively: 
  # install.packages("COCODATA_0.0.0.1.tar.gz", repos = NULL, type ="source")
 # data("norm_ID")
library(COCODATA)

## ----plot_conserved, fig.height=8, fig.width=14--------------------------
data("conserved_region_130_genes_plot_dt") # this loads genes_plot_dt
head(genes_plot_dt)

data("conserved_region_130_tads_plot_dt") # this loads tads_plot_dt
head(tads_plot_dt)

plot_conservedRegions(genes_dt=genes_plot_dt, 
                      tads_dt=tads_plot_dt,
                      dsCat_cols = setNames(c("firebrick3", "navy", "gray50"), c("wt_vs_mut", "norm_vs_tumor", "subtypes")))

## ----plot_volcano, fig.height=6, fig.width=8-----------------------------
data("ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_all_meanCorr_TAD.RData") # this loads all_meanCorr_TAD
data("ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_all_meanLogFC_TAD.RData") # this loads all_meanLogFC_TAD
data("ENCSR489OCU_NCI-H460_40kb_TCGAluad_norm_luad_emp_pval_combined.RData") # this loads emp_pval_combined

plot_volcanoTADsCorrFC(meanCorr=all_meanCorr_TAD, 
                       meanFC=all_meanLogFC_TAD, 
                       comb_pval=emp_pval_combined,
                       tads_to_annot = "chr11_TAD390")

