
### adapted from dryhic package

require(dryhic)
require(magrittr)
my_plot_matrix <- function (mat, tad_coord, resolution, bins_around=NULL, transformation = logfinite,
                            color = colorRampPalette(c("white", "red"))(100), sym = FALSE,
                            borderCol = "darkgrey",
                            borderWidth=2,
                            withDiagBorder=TRUE,
                            plotBorder=TRUE,
                            plotOnly=NULL,
                            trim = 0.01, rotate = FALSE, unit_x_axis = 1000000, label_x_axis = "Genomic Position / Mbp",
                            na.col = "white", ...)
{
  stopifnot(tad_coord %in% colnames(mat))
  stopifnot(tad_coord %in% rownames(mat))
  
  options(scipen = 999)
  rownames(mat) <- colnames(mat) <- rownames(mat) %>% gsub("^.*:", 
                                                           "", .)
  
  
  
  if(!is.null(bins_around)) {
    coord <- tad_coord+c(-bins_around*resolution, bins_around*resolution)  
    startX_line <- bins_around[1]
    endX_line <- (tad_coord[2]-tad_coord[1])/resolution+1-bins_around[2]
  }else {
    coord <- tad_coord
    startX_line <- 0
    endX_line <- (tad_coord[2]-tad_coord[1])/resolution+1
  }
  
  
  lims <- seq(floor(coord[1]/resolution), ceiling(coord[2]/resolution)) * 
    resolution
  i <- rownames(mat) %in% lims
  mat <- as.matrix(mat[i, i])
  mat <- mat[match(lims, rownames(mat)), match(lims, rownames(mat))]
  rownames(mat) <- colnames(mat) <- lims
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
    x[] <- color[cut(c(x), seq(lower, upper, len = length(color) + 
                                 1), include = T)]
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
  plot(NA, type = "n", xlim = c(0, nr), ylim = c(0, nc), xlab = label_x_axis, 
       ylab = "", asp = 1, axes = F, cex.lab = 1.5, ...)
  rasterImage(as.raster(unclass(x)), xleft = 0, xright = nc, 
              ybottom = 0, ytop = nr, interpolate = FALSE)
  axis(1, at = guides_pos, labels = guides/unit_x_axis, cex.axis = 1.5)
  
  if(plotBorder) {
    if(is.null(plotOnly) || grepl("lower_tri", plotOnly))
      segments(x0=startX_line, y0=startX_line, x1=endX_line, y1=startX_line, col=borderCol,lwd=borderWidth)
    if(is.null(plotOnly) || grepl("lower_tri", plotOnly))
      segments(x0=startX_line, y0=startX_line, x1=startX_line, y1=endX_line, col=borderCol,lwd=borderWidth)
    if(is.null(plotOnly) || grepl("upper_tri", plotOnly))
      segments(x0=startX_line, y0=endX_line, x1=endX_line, y1=endX_line, col=borderCol,lwd=borderWidth)
    if(is.null(plotOnly) || grepl("upper_tri", plotOnly))
      segments(x0=endX_line, y0=endX_line, x1=endX_line, y1=startX_line, col=borderCol,lwd=borderWidth)
    if(withDiagBorder)
      segments(x0=startX_line, y0=endX_line, x1=endX_line, y1=startX_line, col=borderCol,lwd=borderWidth)
  }
  
  
  invisible()
}
