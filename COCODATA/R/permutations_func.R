

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#' Gene-to-TAD permutations
#' 
#' Performs the permutation (parallelized).
#'
#' @param g2TADdt Gene-to-TAD dataframe
#' @param RNAdt Gene expression dataframe
#' @param geneIDlist List of gene IDs that should be used
#' @param nClass The number of classes of expression in which genes are shuffled
#' @param withExprClass If shuffling should take place within classes
#' @param TADonly If only genes from TADs should be used
#' @param nSimu The number of permutations
#' @param nCpu The number of CPUs to use
#' @param aggregFun With which function gene expression should be aggregated (across samples)
#' @return The permtuation data
#' @export

get_multiShuffledPositions_vFunct <- function(g2TADdt, RNAdt, geneIDlist, nClass, withExprClass, TADonly, nSimu, nCpu, aggregFun) {
  if(withExprClass) {
    stopifnot(!is.null(RNAdt) & !is.null(nClass))
  }
  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
  g2TADdt$region <- as.character(g2TADdt$region)
  doMC::registerDoMC(nCpu)
  # need to ensure that I get the same order for the genes
  # do the first one
  allT <- get_ShuffledPositions_vFunct(g2TADdt = g2TADdt, RNAdt = RNAdt, geneIDlist = geneIDlist, 
                                      nClass = nClass, TADonly = TADonly, withExprClass = withExprClass, aggregFun=aggregFun) 
  colnames(allT) <- c(colnames(allT)[1], paste0(colnames(allT)[2], "1")) # region1
  genes1 <- allT$entrezID
  if(nSimu >1){
    tmpDT <- foreach::foreach(i=2:nSimu, .combine='cbind') %dopar% {
      if(withExprClass) {
        cat(paste0("... WITH CLASS ", aggregFun, " - shuffle: ", i, "/", nSimu, "\n"))
      } else{
        cat(paste0("... NO CLASS - shuffle: ", i, "/", nSimu, "\n"))
      }
      x <- get_ShuffledPositions_vFunct(g2TADdt = g2TADdt, RNAdt = RNAdt, geneIDlist = geneIDlist, 
                                       nClass = nClass, TADonly = TADonly, withExprClass = withExprClass, aggregFun=aggregFun) 
      stopifnot(all(genes1 == x[,1]))
      x[,2]
    }
    colnames(tmpDT) <- paste0("region", 2:nSimu)
    allT <- cbind(allT, tmpDT)
  }
  return(allT)  
}  



#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#' Gene-to-TAD permutations
#' 
#' Performs the permutation (for 1 permutation).
#'
#' @param g2TADdt Gene-to-TAD dataframe
#' @param RNAdt Gene expression dataframe
#' @param geneIDlist List of gene IDs that should be used
#' @param nClass The number of classes of expression in which genes are shuffled
#' @param withExprClass If shuffling should take place within classes
#' @param TADonly If only genes from TADs should be used
#' @param aggregFun With which function gene expression should be aggregated (across samples)
#' @return The permtuation data
#' @export

### same as _vJune but can pass aggregFun for aggregating the expression values

get_ShuffledPositions_vFunct <- function(g2TADdt, RNAdt, geneIDlist, nClass, withExprClass, TADonly, aggregFun) {
  warning("geneIDlist argument should correspond to rownames of RNAdt")
  warning("duplicated - ambiguous - are removed !!! ")
  duplicatedID <- geneIDlist[duplicated(geneIDlist)]
  RNAdt <- RNAdt[which(! geneIDlist %in% duplicatedID),]
  geneIDlist <- geneIDlist[which(! geneIDlist %in% duplicatedID)]
  rownames(RNAdt) <- geneIDlist
  if(withExprClass) {
    stopifnot(!is.null(RNAdt) & !is.null(nClass))
    stopifnot(!is.null(aggregFun))
  }
  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
  g2TADdt$region <- as.character(g2TADdt$region)
  # take only the genes for which we have their positions
  # and subset the rnaseq data for these genes only
  if(TADonly) {
    # take only the genes that are in TAD
    geneListTAD <- geneIDlist[geneIDlist  %in% g2TADdt$entrezID[grep("_TAD", g2TADdt$region)] ]
    RNAdt <- RNAdt[geneListTAD,]
  } else {
    geneListTAD <- geneIDlist[geneIDlist  %in% g2TADdt$entrezID]
    RNAdt <- RNAdt[geneListTAD,]
  }
  ##########
  ##### DO IT BY SHUFFLING THE LABELS BY CLASS OF EXPRESSION
  ##########
  if(withExprClass) {
    # define expression classes based on median expression
    geneAggregExpression <- data.frame(gene = geneListTAD, expValue = apply(RNAdt, 1, aggregFun))
    rownames(geneAggregExpression) <- NULL
    # rank by expression value
    geneAggregExpression <- geneAggregExpression[order(geneAggregExpression$expValue),]
    # split into 'nClass' groups 
    nR <- nrow(geneAggregExpression)
    nGeneByClass <- rep(nR %/% nClass, nClass)
    # add the remainder (1 to each, to avoid having a lot more in one class)
    if(nR%%nClass > 0)
      nGeneByClass[1:(nR%%nClass)] <- nR %/% nClass + 1
    stopifnot(sum(nGeneByClass) == nR)
    geneClass <- rep(c(1:nClass), nGeneByClass)
    stopifnot(length(geneClass) == nR)
    # add a column to DF with their corresponding class
    geneAggregExpression$class <- geneClass
    #save(geneAggregExpression, file="geneAggregExpression.Rdata")
    stopifnot(length(unique(geneAggregExpression$class)) == nClass)
    # add a column with their initial TAD
    geneAggregExpression$initRegion <- sapply(geneAggregExpression$gene, function(x) g2TADdt$region[g2TADdt$entrezID==x])
    # now, for each class, reshuffle the TAD -> new column with the reshuffled positions
    geneAggregExpression$shuffRegion <- foreach(i_cl = 1:nClass, .combine='c') %do% {
      subDT <- geneAggregExpression[geneAggregExpression$class == i_cl,]
      initPos <- subDT$initRegion
      newPos <- sample(initPos, size=length(initPos), replace=F)  
      stopifnot(all(unique(as.character(initPos)) %in% unique(as.character(newPos))))
      newPos
    }
    shuffGenePosDT <- data.frame(entrezID = geneAggregExpression$gene, 
                                 region = as.character(geneAggregExpression$shuffRegion),
                                 stringsAsFactors = F)
  } else {
    ##########
    ##### JUST SHUFFLE THE LABELS, IRRESPECTIVE OF GENE EXPRESSION
    ##########
    newPos <- sample(g2TADdt$region[g2TADdt$entrezID %in% geneListTAD], replace=F)
    shuffGenePosDT <- data.frame(entrezID = geneListTAD, 
                                 region = as.character(newPos),
                                 stringsAsFactors = F)
  }
  # to be compatible with the functions that use gene2tadDT, return a similar DF
  return(shuffGenePosDT)
}



