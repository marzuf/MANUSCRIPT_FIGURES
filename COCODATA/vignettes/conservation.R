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

