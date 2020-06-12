# Rscript sparse2full_convert.R

options(scipen=199)

sparse_dt <- read.delim("k562_10kb.txt", header=F, stringsAsFactors=F)
head(sparse_dt)

source("sparse2full.R")

full_dt <- sparse2full(sparse.mat=sparse_dt)

dim(full_dt)





outFile <- "k562_10kb_homer2.txt"
chr <- "chr7"
resol <- 10000
full_dt <- data.frame(full_dt)
stopifnot(dim(full_dt)[1]== dim(full_dt)[2])

tmpA <- data.frame(
HiCMatrix = paste0(chr, "-", seq(from=0, length.out=nrow(full_dt), by=resol)),
Regions = paste0(chr, "-", seq(from=0, length.out=nrow(full_dt), by=resol)),
stringsAsFactors=FALSE
)

colnames(full_dt) <- paste0(chr, "-", seq(from=0, length.out=nrow(full_dt), by=resol))

new_dt <- cbind(tmpA, full_dt)

write.table(new_dt, file=outFile, col.names=T, row.names=F, sep="\t", quote=F)

cat("written: ", outFile, "\n")

stopifnot(dim(new_dt)[1]+2 == dim(new_dt)[2])

stop("-ok\n")









outFile <- "k562_10kb_dense.txt"

stopifnot(dim(full_dt)[1] == dim(full_dt)[2])

write.table(full_dt, file=outFile, col.names=F, row.names=F, sep="\t", quote=F)

cat("written: ", outFile, "\n")

outFile <- "k562_10kb_homer.txt"

chr <- "chr7"
resol <- 10000

rownames(full_dt) <- paste0(chr, "-", seq(from=0, length.out=nrow(full_dt), by=resol))
colnames(full_dt) <- rownames(full_dt)
write.table(full_dt, file=outFile, col.names=T, row.names=T, sep="\t", quote=F)
cat("written: ", outFile, "\n")
