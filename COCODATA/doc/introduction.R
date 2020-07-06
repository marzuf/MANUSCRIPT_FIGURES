## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(COCODATA)

## ----FCC_example---------------------------------------------------------
foldChanges_values <- c(-1.2, -0.2, 0.1, 0.3, 0.9, 1.33, 0.4)

get_fcc(foldChanges_values)

