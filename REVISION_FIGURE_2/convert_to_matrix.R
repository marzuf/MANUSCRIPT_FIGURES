file="extract_hic/RWPE1_chr12_obs_KR_20kb.txt"
outfile = "extract_hic/RWPE1_chr12_obs_KR_20kb_matrix.txt"
bin.size=20000

# Rscript convert_to_matrix.R

require(Matrix)

options("scipen"=100, "digits"=6)
chr.data = read.delim(file,header=FALSE)

colnames(chr.data) = c("binsA_bp","binsB_bp","counts")
stopifnot(chr.data$binsA_bp <= chr.data$binsB_bp)  # the data give the upper right

chr.data$binsA = chr.data$binsA_bp/bin.size
chr.data$binsB = chr.data$binsB_bp/bin.size

# they are 0-based and should be one based for SparseMatrix
chr.data$binsA = chr.data$binsA+1
chr.data$binsB = chr.data$binsB+1
stopifnot(chr.data$binsA <= chr.data$binsB)  # the data give the upper right

chr.matrix = sparseMatrix(i=chr.data$binsA,j=chr.data$binsB,x=chr.data$counts)
stopifnot(dim(chr.matrix)[1] == dim(chr.matrix)[2])

# create symmetric matrix
# uplo: The default is "U" unless ‘x’ already has a ‘uplo’ slot
sym_chr.matrix = forceSymmetric(chr.matrix, uplo="U")
i=which( chr.data$binsA != chr.data$binsB)[1]
stopifnot(sym_chr.matrix[chr.data$binsA[i],chr.data$binsB[i]] == sym_chr.matrix[chr.data$binsB[i],chr.data$binsA[i]])
stopifnot(sym_chr.matrix[chr.data$binsA[i],chr.data$binsB[i]] == chr.matrix[chr.data$binsA[i],chr.data$binsB[i]])
stopifnot(dim(sym_chr.matrix)[1] == dim(sym_chr.matrix)[2])

sym_chr.matrix2 <- as.data.frame(as.matrix(sym_chr.matrix))
write.table(sym_chr.matrix2, file=outfile, sep="\t", quote=F, col.names=F, row.names=F)
cat(paste0("... written ", outfile, "\n"))