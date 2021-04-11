# tmp_matcol1 <- matrix(1:25,5)
# tmp_matcol2 <- matrix(51:75,5)
# tmp_matcol1[upper.tri(tmp_matcol1, diag=F)] <- tmp_matcol2[upper.tri(tmp_matcol2)]

# Rscript revision_fig5_combinedPlot.R

outFolder <-file.path("REVISION_FIG5_COMBINEDPLOT")
dir.create(outFolder,recursive = TRUE)

plotType <- "svg"
if(plotType=="png"){
  myWidth <- 400
  myHeight <- 400
}else {
  myWidth <- 6
  myHeight <- 6
}

hicds1 <- "22Rv1"
hicds2 <- "RWPE1"
binsizekb <- "20"
# tad_to_plot <- "chr17_CD174_46720001_46880000"
# tad_to_plot <- "chr12_CD194_54160001_54440000"
tad_to_plot <- "chr7_CD424_116080001_116320000"
binaround <- 0

outFile <- file.path(outFolder, paste0(hicds1, "_", hicds2, "_",tad_to_plot, "_", binsizekb, "kb_combinedPlot.", plotType))
outFile2 <- file.path(outFolder, paste0(hicds1, "_", hicds2, "_",tad_to_plot, "_", binsizekb, "kb_combinedPlot_nodiag.", plotType))
outFile3 <- file.path(outFolder, paste0(hicds1, "_", hicds2, "_",tad_to_plot, "_", binsizekb, "kb_combinedPlot_NAdiag.", plotType))

plotBorder <- TRUE
plotDiag <- TRUE
borderCol <- "darkgrey" 
borderWidth <- 2

diagCol <- rgb(matrix(col2rgb(borderCol), ncol=3), maxColorValue = 255)

plotTit <- gsub("(.+)_(.+)_(.+)_(.+)", paste0("\\1_\\2:\\3-\\4"), tad_to_plot )

subTit <- paste0("lowerTri=", hicds1, "; upperTri=", hicds2)

inFolder <- "REVISION_FIG5_HICPLOT"

matcol1 <- get(load(file.path(inFolder, paste0(hicds1, "_", tad_to_plot, "_", binsizekb, "kb_", binaround, "around_hicplot_matcol.Rdata"))))
matcol2 <- get(load(file.path(inFolder, paste0(hicds2, "_", tad_to_plot, "_", binsizekb, "kb_", binaround, "around_hicplot_matcol.Rdata"))))

matcol <- matcol1
matcol[upper.tri(matcol, diag=F)] <- matcol2[upper.tri(matcol2)]
diag(matcol) <- NA
stopifnot(matcol[lower.tri(matcol, diag=F)] == matcol1[lower.tri(matcol1, diag=F)])
stopifnot(matcol[upper.tri(matcol, diag=F)] == matcol2[upper.tri(matcol2, diag=F)])


x <- matcol

nr = nrow(matcol)
nc = ncol(matcol)
label_x_axis <- ""

do.call(plotType, list(outFile, height=myHeight, width=myWidth))

plot(NA, type = "n", xlim = c(0, nr), ylim = c(0, nc), xlab = label_x_axis, 
     ylab = "", asp = 1, axes = F, cex.lab = 1.5,
     main=plotTit)
rasterImage(as.raster(unclass(x)), xleft = 0, xright = nc, 
                    ybottom = 0, ytop = nr, interpolate = FALSE)
mtext(side=3, text=subTit)

i=1
for(i in 1:nr) {
  tl_x <- i-1
  tl_y <- nr-(i-1)
  br_x <- i
  br_y <- nr-i
  cl_x <- i-1 # corner lower
  cl_y <- nr-i # corner lower
  cu_x <- cl_x+1
  cu_y <- cl_y+1
  polygon(x=c(tl_x, cl_x, br_x),  border=matcol1[i,i],
          y = c(tl_y, cl_y, br_y), col = matcol1[i,i])
  polygon(x=c(tl_x, br_x, cu_x), border=matcol2[i,i],
          y = c(tl_y, br_y, cu_y), col = matcol2[i,i])
}


if(plotBorder) {
  start_x <- binaround
  start_y <- nr-binaround
  end_x <- nc - binaround
  end_y <- binaround
  
  # top segment 
  segments(x0=start_x, y0=start_y, 
           x1=end_x, y1=start_y, 
           col=borderCol,lwd=borderWidth)
  # right
  segments(x0=end_x, y0=start_y, 
           x1=end_x, y1=end_y, 
           col=borderCol,lwd=borderWidth)
  # bottom
  segments(x0=end_x, y0=end_y,
           x1=start_x, y1=end_y, 
           col=borderCol,lwd=borderWidth)
  # bottom
  segments(x0=start_x, y0=end_y, 
           x1=start_x, y1=start_y,
           col=borderCol,lwd=borderWidth)
  
  if(plotDiag)
    segments(x0=start_x, y0=start_y, 
             x1=end_x, y1=end_y,
             col=borderCol,lwd=borderWidth)
}

foo <- dev.off()

cat(paste0("... written: ", outFile, "\n"))

diag(matcol) <- rgb(matrix(col2rgb(borderCol), ncol=3), maxColorValue = 255)
x <- matcol

do.call(plotType, list(outFile2, height=myHeight, width=myWidth))

plot(NA, type = "n", xlim = c(0, nr), ylim = c(0, nc), xlab = label_x_axis, 
     ylab = "", asp = 1, axes = F, cex.lab = 1.5,
     main=plotTit)
rasterImage(as.raster(unclass(x)), xleft = 0, xright = nc, 
            ybottom = 0, ytop = nr, interpolate = FALSE)
mtext(side=3, text=subTit)

if(plotBorder) {
  start_x <- binaround
  start_y <- nr-binaround
  end_x <- nc - binaround
  end_y <- binaround
  
  # top segment 
  segments(x0=start_x, y0=start_y, 
           x1=end_x, y1=start_y, 
           col=borderCol,lwd=borderWidth)
  # right
  segments(x0=end_x, y0=start_y, 
           x1=end_x, y1=end_y, 
           col=borderCol,lwd=borderWidth)
  # bottom
  segments(x0=end_x, y0=end_y,
           x1=start_x, y1=end_y, 
           col=borderCol,lwd=borderWidth)
  # bottom
  segments(x0=start_x, y0=end_y, 
           x1=start_x, y1=start_y,
           col=borderCol,lwd=borderWidth)
  
  if(plotDiag)
    segments(x0=start_x, y0=start_y, 
             x1=end_x, y1=end_y,
             col=borderCol,lwd=borderWidth)
}

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


diag(matcol) <- NA
x <- matcol

do.call(plotType, list(outFile3, height=myHeight, width=myWidth))

plot(NA, type = "n", xlim = c(0, nr), ylim = c(0, nc), xlab = label_x_axis, 
     ylab = "", asp = 1, axes = F, cex.lab = 1.5,
     main=plotTit)
rasterImage(as.raster(unclass(x)), xleft = 0, xright = nc, 
            ybottom = 0, ytop = nr, interpolate = FALSE)
mtext(side=3, text=subTit)

if(plotBorder) {
  start_x <- binaround
  start_y <- nr-binaround
  end_x <- nc - binaround
  end_y <- binaround
  
  # top segment 
  segments(x0=start_x, y0=start_y, 
           x1=end_x, y1=start_y, 
           col=borderCol,lwd=borderWidth)
  # right
  segments(x0=end_x, y0=start_y, 
           x1=end_x, y1=end_y, 
           col=borderCol,lwd=borderWidth)
  # bottom
  segments(x0=end_x, y0=end_y,
           x1=start_x, y1=end_y, 
           col=borderCol,lwd=borderWidth)
  # bottom
  segments(x0=start_x, y0=end_y, 
           x1=start_x, y1=start_y,
           col=borderCol,lwd=borderWidth)
  
  if(plotDiag)
    segments(x0=start_x, y0=start_y, 
             x1=end_x, y1=end_y,
             col=borderCol,lwd=borderWidth)
}

foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

# 1er ligne i=1
# - lower: t_l=0,nr / 0,nr-1 / b_r=1,nr-1
# - upper: t_l / b_r / 1, nr
# 
# 2ème ligne i=2
# - lower:  t_l 1, nr-1 / 1, nr-2 / b_r = 2, nr-2
# - upper: t_l/b_r/2, nr-1
# 
# ligne i
# - lower: t_l = (i-1), (nr-(i-1)) / (i-1), nr-i / b_r=i, nr-i
# - upper; t_l/b_r/i,(nr-(i-1))


