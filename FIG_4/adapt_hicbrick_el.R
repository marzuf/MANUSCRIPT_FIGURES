prep_dt <- function(init_dt, binResol, plotDist, cap_val = NULL) {
  
  stopifnot(ncol(init_dt) == 3)
  colnames(init_dt) <- c("binA", "binB", "count")
  
  out_df <- init_dt
  out_df$row <- out_df$binA/binResol
  out_df$col <- out_df$binB/binResol
  out_df$dist <- out_df$col-out_df$row
  out_df$keep <- out_df$dist >= 0 
  out_df$val <- out_df$count
  out_df <- out_df[abs(out_df$dist) <= plotDist,]
  out_df <- na.omit(out_df)
  
  stopifnot(out_df$row %% 1 == 0)
  stopifnot(out_df$col %% 1 == 0)
  
  Matrix.df <- get(load("Matrix.df.Rdata"))
  
  capped.val <- quantile(out_df$val,cap_val)
  
  out_df$val[out_df$val > capped.val] <- capped.val

  stopifnot(!is.na(out_df))
  
  stopifnot(nrow(out_df) > 0)
  
  return(out_df[,c("row", "col", "val", "dist", "keep")])
  
  
}



# points_score = data.frame(
#   start1=c(0,0,0,0,0,0,
#            10000, 10000,10000,10000,10000,
#            15000, 15000, 15000, 15000,
#            20000, 20000, 20000,
#            25000, 25000, 
#            30000),
#   start2=c(0,10000,15000,20000,25000,30000,
#            10000,15000,20000,25000,30000,
#            15000,20000,25000,30000,
#            20000,25000,30000,
#            25000,30000,
#            30000),
#   score=c(16,20,65,2,2,1,
#           35,4,2,1,4,
#           66,12,3,4,
#           17,1,2,
#           10,1,
#           15),
#   stringsAsFactors = FALSE
# )
# 
# distance <- 60
# binResolution <- 5000
# my_Matrix.df <- sparse_dt
# my_Matrix.df <- points_score
# my_Matrix.df$row <- my_Matrix.df$start1/binResolution
# my_Matrix.df$col <- my_Matrix.df$start2/binResolution
# my_Matrix.df$dist <- my_Matrix.df$col-my_Matrix.df$row
# my_Matrix.df$keep <- my_Matrix.df$dist >= 0 
# my_Matrix.df$val <- my_Matrix.df$score
# my_Matrix.df <- my_Matrix.df[abs(my_Matrix.df$dist) <= distance,]
# my_Matrix.df <- na.omit(my_Matrix.df)
# head(my_Matrix.df)
# my_Matrix.df$val[1] <- 100
# 
# my_Matrix.df <- get(load("Matrix.df.Rdata"))
# Matrix.df_s <- my_Matrix.df
# 
# my_tads <- data.frame(chr="chr1", start=c(1,2900001,6200001,7800001), end=c(800000,4200000,6500000,8300000), stringsAsFactors = FALSE)
# my_tads <- data.frame(chr="chr1", start=c(1,800001), end=c(800000,2900000), stringsAsFactors = FALSE)
# my_tads <- data.frame(chr="chr1", start=c(800001, 2900001), end=c(2900000, 4200000), stringsAsFactors = FALSE)
# 
# my_resol <- 100000

convert_tad_coord <- function(tad_dt, binResol, plotStart){
  stopifnot(c("start", "end", "chr") %in% colnames(tad_dt))
  do.call(rbind, lapply(1:nrow(tad_dt), function(x){
  
  chr <- tad_dt$chr[x]
  newStart <- (tad_dt$start[x]-1)/binResol
  newEnd <- (tad_dt$end[x])/binResol - 1
  
  xshift <- (plotStart-1)/binResol
  # xshift=0
  
  data.frame(
    x=c((newStart+newEnd)/2-xshift, newEnd-xshift, newStart-xshift, (newStart+newEnd)/2-xshift),
    y=c((newEnd-newStart)/2, 0, 0, (newEnd-newStart)/2),
    line.group = c(paste0(chr, "_tad_", x, "_end"), paste0(chr, "_tad_", x, "_end"), 
                   paste0(chr, "_tad_", x, "_start"), paste0(chr, "_tad_", x, "_start")),
    colours="My_Group",
    group="Group.1",
  stringsAsFactors = F
  )
  
}))
}

# mz_Brick_vizart_plot_heatmap(File = file.path(tempdir(),
#                                                 "chr3R-1-10MB-normal-colours-log10-rotate-3-tads.pdf"),
#                                Matrix.df = my_Matrix.df,
#                                tad_dt = my_tads,
#                                # x_coords = "chr3R:1:2900000",
#                                # y_coords = "chr3R:1:2900000",
#                              x_coords = "chr3R:800001:4200000",
#                              y_coords = "chr3R:800001:4200000",
#                                resolution = 100000,
#                                colours = "#230C0F",
#                                FUN = Failsafe_log10,
#                                value_cap = 0.99,
#                                distance = 60,
#                                legend_title = "Log10 Hi-C signal",
#                                palette = "Reds",
#                                width = 15,
#                                height = 5,
#                                line_width = 0.8,
#                                cut_corners = TRUE,
#                                rotate = TRUE,
#                                return_object=TRUE)

  
  

mz_Brick_vizart_plot_heatmap<-function (File, Matrix.df, resolution, x_coords, y_coords,
                                        tad_dt=NULL,
                                        
                                        FUN = NULL, 
                                     value_cap = NULL, distance = NULL, rotate = FALSE, x_axis = TRUE, 
                                     x_axis_title = NULL, y_axis = TRUE, y_axis_title = NULL, 
                                     title = NULL, legend_title = NULL, return_object = FALSE, 
                                     x_axis_num_breaks = 5, y_axis_num_breaks = 5, palette, col_direction = 1, 
#                                     extrapolate_on = NULL, 
x_axis_text_size = 10, y_axis_text_size = 10, 
                                     text_size = 10, legend_title_text_size = 8, legend_text_size = 8, 
                                     title_size = 10, 
                                     # tad_ranges = NULL, 
                                     
                                     group_col = NULL, tad_colour_col = NULL, 
                                     colours = NULL, colours_names = NULL, cut_corners = FALSE, 
                                     highlight_points = NULL, width = 10, height = 6, line_width = 0.5, 
                                     units = "cm", legend_key_width = unit(3, "cm"), legend_key_height = unit(0.5, 
                                                                                                              "cm")) 
{

  
  list_of_coords <- list(x_coords = x_coords, y_coords = y_coords)
  Parsed_string <- HiCBricks:::._Parse_genomic_coordinates(list_of_coords)
  x.coord.parsed <- Parsed_string[["x_coords"]]
  y.coord.parsed <- Parsed_string[["y_coords"]]
  
  
  plotStart <-  as.numeric(Parsed_string[["x_coords"]][["start"]])
  stopifnot(!is.na(plotStart))
  xshift <- (plotStart-1)/resolution
  
  
  
  x.coord.breaks <- HiCBricks:::make_axis_coord_breaks(from = min(Matrix.df$row),
                                                       to = max(Matrix.df$row), how.many = x_axis_num_breaks, two.sample = FALSE)
  
  save(x.coord.breaks, file="toy_x.coord.breaks.Rdata", version=2)
  
  x_axis.coord.labs <- paste0( (x.coord.breaks+xshift)*resolution/10^6, "mb")
  
  # x_axis.coord.labs <- HiCBricks:::Make_axis_labels(Brick = Bricks[[1]], 
  #                                                   resolution = resolution, chr = x.coord.parsed["chr"], 
  #                                                   positions = x.coord.breaks)
  # two.sample <- (rotate & length(Bricks) == 2)
  two.sample <- FALSE
  y.coord.breaks <- HiCBricks:::make_axis_coord_breaks(from = min(Matrix.df$row),
                                                       to = max(Matrix.df$row), how.many = x_axis_num_breaks,
                                                       two.sample = two.sample)
  # y_axis.coord.labs <- HiCBricks:::Make_axis_labels(Brick = Bricks[[1]], 
  #                                                   resolution = resolution, chr = y.coord.parsed["chr"], 
  #                                                   positions = abs(y.coord.breaks))
  Colours <- HiCBricks:::Make_colours(palette = palette, 
#extrapolate_on = extrapolate_on, 
                                      direction = col_direction)
  # two.sample <- (length(Bricks) == 2)
  two.sample <- FALSE
  Matrix.df$rescale <- HiCBricks:::rescale_values_for_colours(Object = Matrix.df, 
                                                              two.sample = two.sample)
  Value.dist <- HiCBricks:::make_colour_breaks(Object = Matrix.df, how.many = length(Colours), 
                                               two.sample = two.sample)
  Legend.breaks.list <- HiCBricks:::get_legend_breaks(Object = Matrix.df, 
                                                      how.many = 5, value_cap = value_cap, colours = Colours, 
                                                      two.sample = two.sample)
  Colour.breaks <- Legend.breaks.list[["col.breaks"]]
  Colour.labs <- Legend.breaks.list[["col.labs"]]
  Colours <- Legend.breaks.list[["cols"]]
  if (rotate) {
    # y.coord.breaks <- y.coord.breaks - min(y.coord.breaks)
    # x.coord.breaks <- x.coord.breaks - min(x.coord.breaks)
    # if (length(Bricks) == 2) {
    #   Upper.tri.map <- Matrix.df[Matrix.df$dist >= 0, ]
    #   Lower.tri.map <- Matrix.df[Matrix.df$dist <= 0, ]
    #   Lower.tri.map$dist <- abs(Lower.tri.map$dist)
    #   Upper.rotated.map <- HiCBricks:::RotateHeatmap(Matrix = Upper.tri.map, 
    #                                                  value.var = "rescale", upper = TRUE)
    #   Lower.rotated.map <- HiCBricks:::RotateHeatmap(Matrix = Lower.tri.map, 
    #                                                  value.var = "rescale", upper = FALSE)
    #   Entire.rotated.map <- rbind(Upper.rotated.map, Lower.rotated.map)
    #   y.coord.breaks <- y.coord.breaks/2
    #   y.coord.breaks <- c(rev(y.coord.breaks) * -1, y.coord.breaks)
    #   y_axis.coord.labs <- c(rev(y_axis.coord.labs), y_axis.coord.labs)
    # }
    # else {
      Upper.tri.map <- Matrix.df[Matrix.df$dist >= 0, ]
      Entire.rotated.map <- HiCBricks:::RotateHeatmap(Matrix = Upper.tri.map, 
                                                      value.var = "rescale", upper = TRUE)
      # y.coord.breaks <- y.coord.breaks/2
    # }
  }
  Brick_theme <- HiCBricks:::Get_heatmap_theme(x_axis = x_axis, y_axis = y_axis, 
                                               text_size = text_size, x_axis_text_size = x_axis_text_size, 
                                               y_axis_text_size = y_axis_text_size, legend_title_text_size = legend_title_text_size, 
                                               legend_text_size = legend_text_size, title_size = title_size, 
                                               legend_key_width = legend_key_width, legend_key_height = legend_key_height)
  Labels <- HiCBricks:::Get_heatmap_titles(title = title, x_axis_title = x_axis_title, 
                                           y_axis_title = y_axis_title, legend_title = legend_title, 
                                           x_coords = x_coords, y_coords = y_coords, rotate = rotate)
  Boundaries.obj <- NULL
  if (!is.null(tad_dt)) {
    # Boundaries.obj <- HiCBricks:::Format_boundaries_normal_heatmap(Bricks = Bricks, 
    #                                                                resolution = resolution, Ranges = tad_ranges, group_col = group_col, 
    #                                                                cut_corners = cut_corners, colour.col = tad_colour_col, 
    #                                                                colours = colours, colours_names = colours_names, 
    #                                                                region.chr = x.coord.parsed["chr"], region.start = as.numeric(x.coord.parsed["start"]), 
    #                                                                region.end = as.numeric(x.coord.parsed["end"]), distance = distance, 
    #                                                                rotate = rotate)
    Boundaries.obj <- convert_tad_coord(tad_dt, resolution, plotStart)
  }
  if (rotate) {
    ids <- xcoords <- ycoords <- NULL
    
    save(Entire.rotated.map, file="toy_Entire.rotated.map.Rdata",version=2)
    
    ThePlot <- ggplot(Entire.rotated.map, aes(x = xcoords, 
                                              y = ycoords))
    ThePlot <- ThePlot + geom_polygon(aes(fill = values, 
                                          group = ids))
    xlims <- c(0, max(Entire.rotated.map[, "xcoords"]))
    ylims <- c(min(Entire.rotated.map[, "ycoords"]), max(Entire.rotated.map[, 
                                                                            "ycoords"]))
    y.coord.breaks <- seq(ceiling(min(Entire.rotated.map[,
                                                         "ycoords"])), ceiling(max(Entire.rotated.map[, "ycoords"])),
                          length.out = y_axis_num_breaks)
    y_axis.coord.labs <- y.coord.breaks * 2
  }
  else {
    Matrix.df$row <- Matrix.df$row - 0.5
    Matrix.df$col <- Matrix.df$col - 0.5
    ThePlot <- ggplot(Matrix.df, aes(x = row, y = col))
    ThePlot <- ThePlot + geom_tile(aes(fill = rescale))
    xlims <- c(min(Matrix.df$row) - 0.5, max(Matrix.df$row) + 
                 0.5)
    ylims <- c(min(Matrix.df$col) - 0.5, max(Matrix.df$col) + 
                 0.5)
  }
  if (!is.null(tad_dt)) {
    line.group <- x <- y <- NULL
    ThePlot <- ThePlot + geom_line(data = Boundaries.obj, 
                                   aes(x = x, y = y, group = line.group, colour = colours), 
                                   size = line_width)
    ThePlot <- ThePlot + scale_colour_manual(values = colours)
  }
  save(xlims, file="toy_xlims.Rdata", version=2)
  save(x_axis.coord.labs, file="toy_x_axis.coord.labs.Rdata", version=2)
  save(x.coord.breaks, file="toy_x.coord.breaks", version=2)
  ThePlot <- ThePlot + scale_x_continuous(limits = xlims, expand = c(0,
                                                                     0), breaks = x.coord.breaks, labels = x_axis.coord.labs)
  ThePlot <- ThePlot + scale_y_continuous(limits = ylims, expand = c(0,
                                                                     0), breaks = y.coord.breaks) #, labels = y_axis.coord.labs)
  ThePlot <- ThePlot + scale_fill_gradientn(legend_title, values = Value.dist, 
                                            breaks = Colour.breaks, labels = Colour.labs, colors = Colours)
  ThePlot <- ThePlot + Brick_theme
  ThePlot <- ThePlot + labs(title = Labels["title"], x = Labels["x_axis"], 
                            y = Labels["y_axis"])
  ggsave(filename = File, plot = ThePlot, width = width, height = height, 
         units = units)
  if (return_object) {
    return(ThePlot)
  }
  else {
    return(TRUE)
  }
}
