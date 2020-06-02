# !!! COPYRIGHT: from HiCCompare package !!!

# will construct the full matrix from the sparse matrix notation
sparse2full <- function(sparse.mat, hic.table = FALSE, column.name = NA) {
  # get the min and max bins if the matrix is a simple sparse matrix
  if (ncol(sparse.mat) == 3) {
    sparse.mat <- as.matrix(sparse.mat)
    bins <- unique(c(sparse.mat[, 1], sparse.mat[, 2]))
    bins <- as.numeric(bins)
    bins <- bins[order(bins)]
    bin.size <- min(diff(bins))
    min.bin <- min(bins)
    max.bin <- max(bins)
    # set sequence of rows/cols for full matrix
    cols <- seq(min.bin, max.bin, by = bin.size)
    # set up matrix to be correct size and row/col names
    mat <- matrix(nrow = length(cols), ncol = length(cols))
    rownames(mat) <- cols
    colnames(mat) <- cols
  }
  # if the matrix is a hic.table
  if (ncol(sparse.mat) > 3 & !hic.table) {
    stop("Enter a sparse matrix with 3 columns or if you are trying to enter a
         hic.table please set hic.table = TRUE")
  }
  if (ncol(sparse.mat) > 3) {
    sparse.mat <- as.matrix(sparse.mat)
    bins <- unique(c(sparse.mat[, 2], sparse.mat[, 5]))
    bins <- as.numeric(bins)
    bins <- bins[order(bins)]
    bins <- unique(bins)
    bin.size <- min(diff(bins))
    min.bin <- min(bins)
    max.bin <- max(bins)
    # set sequence of rows/cols for full matrix
    cols <- seq(min.bin, max.bin, by = bin.size)
    # set up matrix to be correct size and row/col names
    mat <- matrix(nrow = length(cols), ncol = length(cols))
    rownames(mat) <- cols
    colnames(mat) <- cols
  }
  if (hic.table) {
    if (is.na(column.name)) {
      stop("Enter a value for column.name")
    }
    # reconstruct matrix using by converting bin names to matrix cell
    # locations
    col.val <- which(colnames(sparse.mat) == column.name)
    sparse.mat <- sparse.mat[, c(2, 5, col.val)]
    sparse.mat <- apply(sparse.mat, 2, as.numeric)
    # match bin names to column/row number
    sparse.mat[,1] <- match(sparse.mat[,1], cols)
    sparse.mat[,2] <- match(sparse.mat[,2], cols)
    # populate matrix
    rcix <- sparse.mat[,c(1,2)]
    mat[rcix] <- sparse.mat[, 3]
    rcix <- rcix[, c(2,1)]
    mat[rcix] <- sparse.mat[, 3]

  } else {
    # reconstruct matrix using by converting bin names to matrix cell
    # locations
    # match bin names to column/row number
    sparse.mat[,1] <- match(sparse.mat[,1], cols)
    sparse.mat[,2] <- match(sparse.mat[,2], cols)
    # populate matrix
    rcix <- sparse.mat[, 1:2]
    mat[rcix] <- sparse.mat[,3]
    rcix <- rcix[,c(2,1)]
    mat[rcix] <- sparse.mat[,3]
  }
  # replace NAs in matrix with 0
  mat[is.na(mat)] <- 0
  message(paste0("Matrix dimensions: ", dim(mat)[1], "x", dim(mat)[2]))
  return(mat)
}
