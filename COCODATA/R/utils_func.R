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
