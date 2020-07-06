#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
# FLUX: (instead of installing flux package)
#' AUC calculation
#' 
#' Taken from "flux" package (not able to install package after update of R). Compute the Area Under the Curve (AUC) between two vectors of points.
#'
#' @param x Vector of x positions
#' @param y Vector of y positions
#' @param thresh Threshold below which area is not calculated
#' @param By default the data density is artificially increased by adding 100 data points between given adjacent data points.
#' @param By default the vectors in x and y are ordered along increasing x because integration makes no sense with unordered data.
#' @return The AUC
#' @export
auc <- function (x, y, thresh = NULL, dens = 100, sort.x = TRUE) {
    x <- x[!is.na(x)]
    y <- y[!is.na(x)]
    x <- x[!is.na(y)]
    y <- y[!is.na(y)]
    if (sort.x) {
        ord <- order(x)
        x <- x[ord]
        y <- y[ord]
    }
    idx = 2:length(x)
    x <- as.vector(apply(cbind(x[idx - 1], x[idx]), 1, function(x) seq(x[1], 
        x[2], length.out = dens)))
    y <- as.vector(apply(cbind(y[idx - 1], y[idx]), 1, function(x) seq(x[1], 
        x[2], length.out = dens)))
    if (!is.null(thresh)) {
        y.0 <- y <= thresh
        y[y.0] <- thresh
    }
    idx = 2:length(x)
    integral <- as.double((x[idx] - x[idx - 1]) %*% (y[idx] + 
        y[idx - 1]))/2
    integral
}


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

#' Stouffer's method for p-values combination
#'
#' Function that combines p-values using Stouffer's method.
#'
#' @param ps The p-values to be combined
#' @return The normalized dataframe
#' @export

stouffer <- function(ps, two.tails=FALSE) {
	if(two.tails) ps = ps / 2
	# transform p-values into quantiles of standard normal
	qps = qnorm(1-ps, lower.tail = TRUE)
	# take the average of the quantiles
	iqps = sum(qps) / sqrt(length(qps))
	# find p-val of the average
	p = 1 - pnorm(iqps, lower.tail = TRUE)
	if(two.tails) p = 2 * p
	return(p)
}

#the individual p-values are transformed into the quantiles of a standard normal
#the p-val of the average of the quantiles is then found

# Fisher's method treates large and small p-values asymmetrically
# e.g. combining 0.999 and 0.001 gives 0.008
# is asymmetrically sensitive to small p-values

# Z-transform test: one-to-one mapping of the standard normal curve of the p-value
# Z = standard deviate (= a number drawn from from normal distribution with mean 0 and sd 1)
# the test converts p-values into standard normal deviates
# sum of the Z divided by square root of the number of tests has a std normal distribution
# no asymmetry problem

# weighted version of the Z-transform test: a weight can be assigned to each test
# if each test is given equal weight, this reduces to Z-transform test
# ideally, each study is weighted proportional to the inverse of its error variance
# (i.e. by the reciprocal of its squared standard error)
# more generally, the weights should be the inverse of the squared standard error fo the effect
# size estimate for each study


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

#' Fold-change concordance (FCC) score
#'
#' Function that calculates the FCC score for a vector of fold-change (FC) values.
#'
#' @param fc_vect Vector of fold-changes
#' @return The FCC score
#' @export

get_fcc <- function(fc_vect) {
  (2* sum(fc_vect < 0)/length(fc_vect) -1) *  (2* sum(abs(fc_vect[fc_vect<0]))/sum(abs(fc_vect)) -1)
}


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

#' RatioDown score
#'
#' Function that calculates ratio of down-regulated genes for a vector of fold-change (FC) values.
#'
#' @param fc_vect Vector of fold-changes
#' @return The ratioDown score
#' @export
get_ratioDown <- function(fc_vect) {
  sum(fc_vect < 0)/length(fc_vect) 
}


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


#' RatioFC score
#'
#' Function that calculates ratio of negative fold-change (FC) for a vector of FC values.
#'
#' @param fc_vect Vector of fold-changes
#' @return The ratioFC score
#' @export
get_ratioFC <- function(fc_vect) {
  sum(abs(fc_vect[fc_vect<0]))/sum(abs(fc_vect))
}








