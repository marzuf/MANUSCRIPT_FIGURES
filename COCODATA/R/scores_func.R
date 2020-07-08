#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
# FLUX: (instead of installing flux package)
#' AUC calculation
#' 
#' Taken from "flux" package (not able to install package after update of R). Compute the Area Under the Curve (AUC) between two vectors of points.
#'
#' @param x Vector of x positions.
#' @param y Vector of y positions.
#' @param thresh Threshold below which area is not calculated.
#' @param dens By default the data density is artificially increased by adding 100 data points between given adjacent data points.
#' @param By default the vectors in x and y are ordered along increasing x because integration makes no sense with unordered data.
#' @return The AUC.
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
#' @param ps The p-values to be combined.
#' @return The normalized dataframe.
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
#' @param fc_vect Vector of fold-changes.
#' @return The FCC score.
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
#' @param fc_vect Vector of fold-changes.
#' @return The ratioDown score.
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
#' @param fc_vect Vector of fold-changes.
#' @return The ratioFC score.
#' @export
get_ratioFC <- function(fc_vect) {
  sum(abs(fc_vect[fc_vect<0]))/sum(abs(fc_vect))
}



#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


#' AUC ratio 
#'
#' Function that calculates the AUC ratio for FCC cumsum curves (observed/permutation).
#'
#' @param fcc_vect Vector of observed FCCs.
#' @param fcc_permDT Dataframe with permutation FCCs (each column is a permutation).
#' @param permQt Quantile of the permutation (default: 0.95) used for computing the AUC ratio.
#' @param doPlot If true, plot the cumsum wave curves.
#' @param plotCex Cex value for the axis labs and titles (default: 1.2).
#' @param ploygonPermutCol Color for the area of permutation value (default: grey color).
#' @param qt95PermutCol Color for the line for the permQt permutation values (default: darker grey).
#' @param lwdObs Width of the line of the observed data (default: 1.2).
#' @param colObs Color for the observed data (default: "darkred").
#' @param ... Other arguments are passed to plot() function.
#' @return (invisible) A list containing the observed AUC (observed_auc), the AUC for the permutation (permut_auc) and the AUC ratio (auc_ratio).
#' @export

get_auc_ratio <- function(fcc_vect, fcc_permDT, permQt=0.95, doPlot=FALSE,
  							plotCex=1.2, polygonPermutCol=rgb(188/255,188/255,188/255, 0.3), qt95PermutCol=rgb(140/255,140/255,140/255),
								lwdObs=1.2, colObs="darkred", ...) {
  pointObsCol <- colObs
  # sort decreasing the FCC
  obs_fcc <- sort(fcc_vect, decreasing = TRUE)
  permut_FCC_unsort <- fcc_permDT
  permut_FCC <- apply(permut_FCC_unsort, 2, sort, decreasing=TRUE)
  rownames(permut_FCC) <- NULL
  cat(paste0("... found obs. TADs:\t", length(obs_fcc), "\n"))
  cat(paste0("... found permut. TADs:\t", nrow(permut_FCC), "\n"))
  maxTADs <- max(c(length(obs_fcc), nrow(permut_FCC)))
  
  maxRankPlot <- ceiling(maxTADs/1000)*1000
  x_val <- c(1:length(obs_fcc))
  cumsum_permut_dt <- apply(permut_FCC, 2, function(x) cumsum(x))
  qt95Permut_cumsum <- apply(cumsum_permut_dt, 1, quantile, probs=permQt) # for each rank, take the quantile
  
  obs_cumsum <- cumsum(obs_fcc)
  
  auc_obs <- auc(x = x_val, y = obs_cumsum)
  auc_permutQt <- auc(x = x_val, y = qt95Permut_cumsum)
  auc_ratioQt <- auc_obs/auc_permutQt
  
  
  fcc_auc_ratios <- list(
    observed_auc=auc_obs,
    permut_auc=auc_permutQt,
    auc_ratio = auc_ratioQt
  )
  
  if(doPlot){
    pct_inc_qt <- round(auc_obs/auc_permutQt,2)
    my_xlab <- paste0("TADs ranked by FCC")
	my_ylab <- paste0("FCC cumulative sum")
    par(bty="l")
    plot(obs_cumsum ~ x_val,
         type = "l",
         pch = 16, cex = 0.7,
         xlab= my_xlab, 
         ylab= my_ylab,
         cex.main = plotCex,
         cex.lab = plotCex,
         cex.axis = plotCex,
         col = colObs,
         axes=FALSE,
         lwd = lwdObs,
         bty="l", ...)
    box()
    axis(2, cex.axis=plotCex, cex.lab=plotCex)
    axis(1, cex.axis=plotCex, cex.lab=plotCex, at = seq(from=0, to=maxRankPlot, by=maxRankPlot/10)) # we used 2000 and 200 hard-coded
    polygon(x = c(x_val, rev(x_val)), 
            y = c( apply(cumsum_permut_dt, 1, function(x) min(x)), rev(apply(cumsum_permut_dt, 1, function(x) max(x)))),
            border=NA,
            col = polygonPermutCol)
    lines(
      x = x_val,
      y= qt95Permut_cumsum,
      col = qt95PermutCol
    )
    legend("topleft",
           xjust=0.5, yjust=0,
           pch = c(16, 15, -1 ),
           lty=c(-1,-1,1),
           legend = c(paste0("observed (n=", length(obs_fcc), ")"), "min-max permut.", paste0(permQt, "-qt permut.")), 
           pt.cex = c(0.7, 2),
           col = c(pointObsCol, polygonPermutCol, qt95PermutCol),
           bty="n")
    legtxt <- as.expression(bquote(frac(AUC[obs.], AUC[permut.]) == .(pct_inc_qt)))
    legend("right", legend=c(legtxt), bty="n")
  }
  invisible(fcc_auc_ratios)
}








