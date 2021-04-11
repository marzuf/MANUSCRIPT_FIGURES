options(scipen = 999)

### adapted from dryhic package





require(dryhic)
require(magrittr)
my_plot_matrix <- function (mat, tad_coord, 
                            resolution, 
                            bins_around=NULL, 
                            transformation = logfinite,
                            color = colorRampPalette(c("white", "red"))(100), 
                            na.col = "white",
                            sym = FALSE,
                            
                            borderCol = "darkgrey",
                            borderWidth=2,
                            withDiagBorder=TRUE,
                            plotBorder=TRUE,
                            plotOnly=NULL,
                            checkSim=TRUE,
                            plotWithLeg=FALSE,
                            trim = 0.01, rotate = FALSE, unit_x_axis = 1000000, label_x_axis = "Genomic Position / Mbp",
                            
                            # for plotting the legend
                            minColRange=NULL, maxColRange=NULL, 
                            subTit="", legTit ="",
                            matrixMargins=c(5.1,4.1,4.1,2.1),
                            legMargins=c(5.1,0.5,4.1,0.5),
                            legCoords=  c(1,1,1.5,2), #  xl,yb,xr,yt
                            legWidth =0.1,
                            
                            saveMatFile = NULL,
                            saveMatcolFile = NULL,
                            
                             ...)
{
  
  options(scipen = 999)
  
    
  # ADDED >>>
  # take 0-based start
  tad_coord[1] <- tad_coord[1]-1
  stopifnot(tad_coord %% resolution == 0)
  # ## if I give 1-80000 -> I want to plot the 4 following bins: 0,20000,40000,60000
  stopifnot(tad_coord[1] %in% colnames(mat))
  stopifnot(tad_coord[1] %in% rownames(mat))
  
  cat(paste0("tad_coord[2]=",tad_coord[2], "\n"))
cat(paste0("head col=",paste0(head(colnames(mat)), collapse=","), "\n"))
cat(paste0("tail col=",paste0(tail(colnames(mat)), collapse=","), "\n"))
cat(paste0("head col=",paste0(head(rownames(mat)), collapse=","), "\n"))
cat(paste0("tail col=",paste0(tail(rownames(mat)), collapse=","), "\n"))

cat(paste0("resolution=",resolution, "\n"))

  stopifnot( (tad_coord[2]-resolution) %in% colnames(mat))
  stopifnot( (tad_coord[2]-resolution) %in% rownames(mat))
  # <<<

  rownames(mat) <- colnames(mat) <- rownames(mat) %>% gsub("^.*:", 
                                                           "", .)
  
  nBins <- tad_coord[2]/resolution -  1-tad_coord[1]/resolution+1 
  
  
  if(!is.null(bins_around)) {
    coord <- tad_coord+c(-bins_around[1]*resolution, bins_around[2]*resolution)  
    cat(paste0("coord = ",paste0(coord, collapse=","),"\n") )
    # startX_line <- bins_around[1]
    # endX_line <- (tad_coord[2]-tad_coord[1])/resolution+1-bins_around[2]
    # endX_line <- (tad_coord[2]-resolution-tad_coord[1])/resolution+1-bins_around[2]
    
    # endX_line <- (tad_coord[2]-tad_coord[1])/resolution + startX_line
    
    # cat(paste0("startX_line = ",paste0(startX_line, collapse=","),"\n") )
    # cat(paste0("endX_line = ",paste0(endX_line, collapse=","),"\n") )
    
    
    bl_x <- bins_around[1]
    bl_y <- bins_around[2]
    
    tl_x <- bl_x
    tl_y <- bins_around[2] + nBins
    
    tr_x <- bins_around[1] + nBins
    tr_y <- tl_y
    
    br_x <- tr_x
    br_y <- bl_y
    cat(paste0("nBins=", nBins, "\n"))
    
  }else {
    coord <- tad_coord
    # startX_line <- 0
    # endX_line <- (tad_coord[2]-tad_coord[1])/resolution+1
    # endX_line <- (tad_coord[2]-resolution-tad_coord[1])/resolution+1
    
    bl_x <- 0
    bl_y <- 0
    
    tl_x <- bl_x
    tl_y <-  nBins
    
    tr_x <- nBins
    tr_y <- tl_y
    
    br_x <- tr_x
    br_y <- bl_y
    cat(paste0("nBins=", nBins, "\n"))
    
  }
  # ADDED >>>
  # tad input 40001-120000 -> 2-5
  bin_start <- coord[1]/resolution
  bin_end <- coord[2]/resolution -  1
  
  
  # <<<
  # 1->40'0000 => 0/20000->40'000/20000 => 0->20 => to plot: 0,1
  lims <- seq(floor(bin_start), ceiling(bin_end)) * 
    resolution
  
  cat(paste0("# plotted: ", length(lims), "\n"))
  cat(paste0("plotted: ", paste0(lims/resolution, collapse=","), "\n"))
  
  i <- rownames(mat) %in% lims
  mat <- as.matrix(mat[i, i])
  # mat <- mat[match(lims, rownames(mat)), match(lims, rownames(mat))]
  mat <- mat[match(lims, rownames(mat)), match(lims, rownames(mat)), drop=FALSE]
  rownames(mat) <- colnames(mat) <- lims
  
  cat(paste0("colnames(mat): ", paste0(colnames(mat), collapse=","), "\n"))
  
  # ADDED >>>
  if(is.null(bins_around)) {
    stopifnot(min(lims) == (tad_coord[1]))
    stopifnot(max(lims) == tad_coord[2]-resolution)
    stopifnot(length(lims) == (tad_coord[2]-tad_coord[1])/resolution)
    
  } else {
    cat("TODOHERE") # Rscript revision_fig5_hicPlot.R RWPE1 chr12_CD194 chr12 40001 120000 0 20000 pngand test 
    stopifnot(min(lims) == (tad_coord[1]-bins_around[1]*resolution))
    cat(paste0("max(lims) = ", max(lims), "\n"))
    cat(paste0("tad_coord[2] = ", tad_coord[2], "\n"))
    cat(paste0("bins_around[2] = ", bins_around[2], "\n"))
    stopifnot(max(lims) == (tad_coord[2]-resolution+bins_around[2]*resolution))
    cat(paste0("length(lims) = ", length(lims), "\n"))
    cat(paste0("sum(bins_around)  = ", sum(bins_around) , "\n"))
    cat(paste0(" (tad_coord[2]-tad_coord[1])/resolution) = ",  (tad_coord[2]-tad_coord[1])/resolution), "\n")
    stopifnot(length(lims) == sum(bins_around) + (tad_coord[2]-tad_coord[1])/resolution)
  }
  
  # <<<
  if(checkSim) stopifnot(isSymmetric(mat))
  
  if(!is.null(saveMatFile)) {
    dir.create(dirname(saveMatFile), recursive = TRUE)
    save(mat, file=saveMatFile, version=2)
    cat(paste0("... written: ", saveMatFile, "\n"))
  }
  
  if (rotate) 
    mat <- mat[nrow(mat):1, ]
  guides <- pretty(x = rownames(mat) %>% as.numeric)
  guides_pos <- data.frame(y = 1:nrow(mat), x = rownames(mat) %>% 
                             as.numeric) %>% lm(y ~ x, .) %>% predict(newdata = data.frame(x = guides))
  # par(mar = c(4, 0, 0, 0), pty = "s")
  par(mar = c(4, 0, 3, 0), pty = "s")
  if (trim > 0) {
    trim <- as.matrix(mat) %>% c %>% quantile(c(trim/2, 1 - 
                                                  trim/2), na.rm = T)
    mat[mat < trim[1]] <- trim[1]
    mat[mat > trim[2]] <- trim[2]
  }
  x <- as.matrix(mat) %>% transformation %>% as.matrix
  if (sym) {
    upper <- max(abs(x), na.rm = T)
    lower <- -upper
  }  else {
    lower <- min(x, na.rm = T)
    upper <- max(x, na.rm = T)
  }
  if (max(x, na.rm = T) == min(x, na.rm = T)) {
    x[] <- color[round(length(color)/2)]
  }  else {
    if(is.null(minColRange)) 
      minColRange <- lower
    if(is.null(maxColRange)) 
      maxColRange <- upper
    
    
    x[] <- color[cut(c(x), seq(minColRange, maxColRange, len = length(color) + 
                                 1), include = T)]
  }
  
  if(!is.null(saveMatcolFile)) {
    dir.create(dirname(saveMatcolFile), recursive = TRUE)
    save(x, file=saveMatcolFile, version=2)
    cat(paste0("... written: ", saveMatcolFile, "\n"))
  }
  
  
  if(!is.null(plotOnly)) {
    stopifnot(plotOnly %in% c("upper_tri", "upper_tri_noDiag", "lower_tri", "lower_tri_noDiag"))  
    
    if(grepl( "upper_tri", plotOnly )) {
      x[lower.tri(x, diag = plotOnly=="upper_tri_noDiag")] <- NA
    }
    if(grepl( "lower_tri", plotOnly )) {
      x[upper.tri(x, diag = plotOnly=="lower_tri_noDiag")] <- NA
    }
  }
  
  
  x[is.na(x)] <- na.col
  range_pos <- as.numeric(rownames(x)) %>% range
  nr <- nrow(x)
  nc <- ncol(x)
  d <- sqrt(nr^2 + nc^2)
  d2 <- 0.5 * d
  
  if(plotWithLeg) {
    layout(matrix(1:2,nrow=1),widths=c(1-legWidth, legWidth))
    colfunc <- colorRampPalette(c("blue", "white", "red"))
    par(mar=matrixMargins,xpd=T)
    
  }
  
  plot(NA, type = "n", xlim = c(0, nr), ylim = c(0, nc), xlab = label_x_axis, 
       ylab = "", asp = 1, axes = F, cex.lab = 1.5, ...)
  rasterImage(as.raster(unclass(x)), xleft = 0, xright = nc, 
              ybottom = 0, ytop = nr, interpolate = FALSE)
  axis(1, at = guides_pos, labels = guides/unit_x_axis, cex.axis = 1.5)
  
  if(plotBorder) {
    if(is.null(plotOnly) || grepl("lower_tri", plotOnly))
      # segments(x0=startX_line, y0=startX_line, x1=endX_line, y1=startX_line, col=borderCol,lwd=borderWidth) # bottom
    # segments(x0=startX_line, y0=startX_line, x1=endX_line, y1=startX_line, col=borderCol,lwd=borderWidth)
    segments(x0=bl_x, y0=bl_y, x1=br_x, y1=br_y, col=borderCol,lwd=borderWidth)
    if(is.null(plotOnly) || grepl("lower_tri", plotOnly))
      # segments(x0=startX_line, y0=startX_line, x1=startX_line, y1=endX_line, col=borderCol,lwd=borderWidth) # left
    segments(x0=bl_x, y0=bl_y, x1=tl_x, y1=tl_y, col=borderCol,lwd=borderWidth)
    if(is.null(plotOnly) || grepl("upper_tri", plotOnly))
      # segments(x0=startX_line, y0=endX_line, x1=endX_line, y1=endX_line, col=borderCol, lwd=borderWidth) # top
    segments(x0=tl_x, y0=tl_y, x1=tr_x, y1=tr_y, col=borderCol, lwd=borderWidth) 
    if(is.null(plotOnly) || grepl("upper_tri", plotOnly))
      # segments(x0=endX_line, y0=endX_line, x1=endX_line, y1=startX_line, col=borderCol,lwd=borderWidth) 
    segments(x0=tr_x, y0=tr_y, x1=br_x, y1=br_y, col=borderCol,lwd=borderWidth) # right
    if(withDiagBorder){
      # segments(x0=startX_line, y0=endX_line, x1=endX_line, y1=startX_line, col=borderCol,lwd=borderWidth)
      if(rotate) {
        segments(x0=bl_x, y0=bl_y, x1=tr_x, y1=tr_y, col=borderCol,lwd=borderWidth)
      }  else {
        segments(x0=tl_x, y0=tl_y, x1=br_x, y1=br_y, col=borderCol,lwd=borderWidth)
      }
    }
  }
  
  
  if(plotWithLeg) {
    mtext(side=3, text=subTit)
    
    xl <- legCoords[1]
    yb <- legCoords[2]
    xr <- legCoords[3]
    yt <- legCoords[4]
    
    par(mar=legMargins)
    plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(1,2),xaxt="n",yaxt="n",bty="n")
    rect(
      xl,
      head(seq(yb,yt,(yt-yb)/10),-1),
      xr,
      tail(seq(yb,yt,(yt-yb)/10),-1),
      col=colfunc(10)
    )
    
    # abline(h=1, v=1)
    
    leg_txt <- round(seq(from=minColRange, to=maxColRange, length.out=10),2)
    
    # leg_labelpos <- tail(seq(yb,yt,(yt-yb)/10),-1)-0.05
    leg_labelpos <- tail(seq(yb,yt,(yt-yb)/10),-1)-0.15
    
    mtext(leg_txt,side=2,at=leg_labelpos,las=2,cex=0.7)
    mtext(bquote(paste("", bold(.(legTit)), "")), side=2,at=max(leg_labelpos)+(1*leg_labelpos[2]-leg_labelpos[1]),las=2,cex=0.9, line=-2.5)
    # mtext(bquote(paste("", bold(.(legTit)), "")), side=2,at=max(leg_labelpos)+(2*leg_labelpos[2]-leg_labelpos[1]),las=2,cex=0.9)
    
  }
  
  invisible()
}
