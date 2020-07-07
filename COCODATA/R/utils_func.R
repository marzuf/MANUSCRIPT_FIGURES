#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#' Print text and log
#'
#' Function that "cat" a message and print it to log file.
#'
#' @param mytext The text to be printed.
#' @param mylogf Location where to print the text (e.g. log file or NULL) 
#' @return The ratioFC score
#' @export
printAndLog <- function(mytext, mylogf) {
  cat(mytext)
  cat(mytext, append=T, file = mylogf)
}


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#' Return the first lines/columns of a dataframe
#'
#' Returns the first lines and columns of a dataframe (as many columns as lines)
#'
#' @param dt A two-dimensional dataframe or matrix.
#' @param ntop The maximal number of lines/columns to display
#' @return The first rows/columns of the dataframe
#' @export
head_sq <- function(dt, ntop=5) {
  nRow <- min(c(ntop, nrow(dt)))	
  nCol <- min(c(ntop, ncol(dt)))
  dt[seq_len(nRow),seq_len(nCol)]
}
