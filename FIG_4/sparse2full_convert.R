# Rscript sparse2full_convert.R

sparse_dt <- read.delim("k562_10kb.txt", header=F, stringsAsFactors=F)
head(sparse_dt)

source("sparse2full.R")

full_dt <- sparse2full(sparse.mat=sparse_dt)

dim(full_dt)

outFile <- "k562_10kb_dense.txt"

stopifnot(dim(full_dt)[1] == dim(full_dt)[2])

write.table(full_dt, file=outFile, col.names=F, row.names=F, sep="\t", quote=F)

cat("written: ", outFile, "\n")
