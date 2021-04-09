#################################################################### 
################################## HiTC
#################################################################### 
createHTC = function(file, bin.size, chr, dim = -1, reindex = FALSE, header=TRUE){
  cat("...... start createHTC \n")
  options("scipen"=100, "digits"=6)
  chr.data = read.delim(file,header=header)
  colnames(chr.data) = c("binsA","binsB","counts")
  # data from Rao et al. indicate bins by genome coordinates, they need to be turned into indeces
  if(reindex){
    chr.data$binsA = chr.data$binsA/bin.size
    chr.data$binsB = chr.data$binsB/bin.size
  }
  
  chr.data$binsA = chr.data$binsA+1
  chr.data$binsB = chr.data$binsB+1
  
  chr.matrix = sparseMatrix(i=chr.data$binsA,j=chr.data$binsB,x=chr.data$counts)
  
  # resize matrix
  if(dim == -1)
    dim = max(ncol(chr.matrix),nrow(chr.matrix))
  
  if(ncol(chr.matrix) < dim){
    nadd = dim-ncol(chr.matrix)
    for(i in 1:nadd)
      chr.matrix = cbind(chr.matrix,rep(0,nrow(chr.matrix)))
  }
  
  if(nrow(chr.matrix) < dim){
    nadd = dim-nrow(chr.matrix)
    for(i in 1:nadd)
      chr.matrix = rbind(chr.matrix,rep(0,ncol(chr.matrix)))
  }
  
  cat(c(dim(chr.matrix),"\t"))
  cat("\n")
  
  stopifnot(ncol(chr.matrix) == dim & nrow(chr.matrix) == dim)
  
  # create symmetric matrix
  # uplo: The default is "U" unless ‘x’ already has a ‘uplo’ slot
  chr.matrix = forceSymmetric(chr.matrix)
  # before taking the upper as reference, just ensure that 3 column format stores the upper triangle matrix
  stopifnot(all(chr.data$binsA <= chr.data$binsB))
  
  #************ TMP
  #tmp_matrix <- as.data.frame(as.matrix(chr.matrix))
  #write.table(tmp_matrix, file = "foo_test_matrix_list.txt", quote=F, sep="\t", row.names=F, col.names=F)
  #*****************
  
  ranges = IRanges(start=seq(1,(nrow(chr.matrix)*bin.size),bin.size),width=rep(bin.size,nrow(chr.matrix)))
  chr.granges = GRanges(seqnames=Rle(c(chr),c(nrow(chr.matrix))), ranges = ranges)
  
  # give a name to each row, here just a consecutive number
  rownames(chr.matrix) = c(1:nrow(chr.matrix))
  colnames(chr.matrix) = c(1:ncol(chr.matrix))
  
  names(chr.granges) = rownames(chr.matrix)
  
  HTC = HTCexp(chr.matrix,chr.granges,chr.granges)
  
  return(HTC)
}
require(Matrix)
require(HiTC)

# htc_dt <- createHTC(file = "extract_hic/RWPE1_chr12_obs_KR_20kb.txt", bin.size=20000, chr="chr12", dim = -1, reindex = TRUE, header=FALSE)
# 
# toptad1  <- extractRegion(htc_dt, chr="chr12", from=54160001, to=54440000)
# 
# png("chr12_RWPE1_hitc_plot.png")
# plot(toptad1, col.pos=c("white", "orange", "red", "black"))
# dev.off()

#################################################################### 
##################################  Sushi
#################################################################### 

htc_mat <- read.delim( "extract_hic/RWPE1_chr12_obs_KR_20kb_matrix.txt", header=F, stringsAsFactors = FALSE)
dim(htc_mat)

require(Sushi)
# plotHic(hicdata, chrom, chromstart, chromend, max_y = 30, zrange = NULL,
#         palette = SushiColors(7), flip = FALSE)
# 
# data(Sushi_HiC.matrix)
# 
# chrom            = "chr11"
# chromstart       = 500000
# chromend         = 5050000

# phic = plotHic(Sushi_HiC.matrix,chrom,chromstart,chromend,max_y = 20,zrange=c(0,28),palette = topo.colors,flip=FALSE)

# labelgenome(chrom,chromstart,chromend,side=1,scipen=20,n=4,scale="Mb",edgeblankfraction=0.20,line=.18,chromline=.5,scaleline=0.5)

# addlegend(phic[[1]],palette=phic[[2]],title="score",side="right",bottominset=0.4,topinset=0,xoffset=-.035,labelside="left",width=0.025,title.offset=0.035)

colnames(htc_mat) <- seq(from=0, by=20000, length.out=ncol(htc_mat))
rownames(htc_mat) <- seq(from=0, by=20000, length.out=nrow(htc_mat))

png("chr12_RWPE1_sushi_plot.png")
phic = plotHic(htc_mat,chrom="chr12",
               chromstart=54160000,chromend=54440000)
dev.off()

#################################################################### 
#################################################################### 
require(dryhic)
png("chr12_RWPE1_dryhic_plot.png")
plot_matrix(
  mat=htc_mat,
  coord=c(54160000,54440000),
  resolution=20000)
dev.off()

#################################################################### 
#################################################################### 
source("my_plot_matrix.R")
my_plot_matrix(mat = htc_mat, 
               tad_coord = c(54160000,54440000),
               resolution = 20000, 
               bins_around=NULL) 
                            
                            
