Matrix.df <- get(load("Matrix.df.Rdata"))

capped.val <- quantile(Matrix.df$val,value_cap)
Matrix.df$val[Matrix.df$val > capped.val] <- capped.val
}
Matrix.df$dist <- Matrix.df$col - Matrix.df$row
Matrix.df$keep <- FALSE
if(i == 1){
  Matrix.df$keep[Matrix.df$dist >= 0] <- TRUE
}else{
  Matrix.df$keep[Matrix.df$dist <= 0] <- TRUE
}
Matrix.df.list[[i]] <- Matrix.df
}

init_Brick_vizart_plot_heatmap<-function (File, Bricks, resolution, x_coords, y_coords, FUN = NULL, 
                                        value_cap = NULL, distance = NULL, rotate = FALSE, x_axis = TRUE, 
                                        x_axis_title = NULL, y_axis = TRUE, y_axis_title = NULL, 
                                        title = NULL, legend_title = NULL, return_object = FALSE, 
                                        x_axis_num_breaks = 5, y_axis_num_breaks = 5, palette, col_direction = 1, 
                                        extrapolate_on = NULL, x_axis_text_size = 10, y_axis_text_size = 10, 
                                        text_size = 10, legend_title_text_size = 8, legend_text_size = 8, 
                                        title_size = 10, tad_ranges = NULL, group_col = NULL, tad_colour_col = NULL, 
                                        colours = NULL, colours_names = NULL, cut_corners = FALSE, 
                                        highlight_points = NULL, width = 10, height = 6, line_width = 0.5, 
                                        units = "cm", legend_key_width = unit(3, "cm"), legend_key_height = unit(0.5, 
                                                                                                                 "cm")) 
{
  if (!is.list(Bricks)) {
    stop("Bricks expects an argument of type list.", " Please refer to the vignette to understand the parameter.")
  }
  Matrix.df <- HiCBricks:::Get_one_or_two_brick_regions(Bricks = Bricks, 
                                                        resolution = resolution, x_coords = x_coords, y_coords = y_coords, 
                                                        distance = distance, value_cap = value_cap, FUN = FUN)
  save(Matrix.df, file="Matrix.df.Rdata", version=2)
  
  if (nrow(Matrix.df) == 0) {
    stop("The matrix was empty!")
  }
  list_of_coords <- list(x_coords = x_coords, y_coords = y_coords)
  Parsed_string <- HiCBricks:::._Parse_genomic_coordinates(list_of_coords)
  
  save(Parsed_string, file="Parsed_string.Rdata",version=2)
  x.coord.parsed <- Parsed_string[["x_coords"]]
  y.coord.parsed <- Parsed_string[["y_coords"]]
  x.coord.breaks <- HiCBricks:::make_axis_coord_breaks(from = min(Matrix.df$row), 
                                                       to = max(Matrix.df$row), how.many = x_axis_num_breaks, 
                                                       two.sample = FALSE)
  
  save(x.coord.breaks, file="x.coord.breaks.Rdata", version=2)
  
  x_axis.coord.labs <- HiCBricks:::Make_axis_labels(Brick = Bricks[[1]], 
                                                    resolution = resolution, chr = x.coord.parsed["chr"], 
                                                    positions = x.coord.breaks)
  two.sample <- (rotate & length(Bricks) == 2)
  y.coord.breaks <- HiCBricks:::make_axis_coord_breaks(from = min(Matrix.df$row), 
                                                       to = max(Matrix.df$row), how.many = x_axis_num_breaks, 
                                                       two.sample = two.sample)
  y_axis.coord.labs <- HiCBricks:::Make_axis_labels(Brick = Bricks[[1]], 
                                                    resolution = resolution, chr = y.coord.parsed["chr"], 
                                                    positions = abs(y.coord.breaks))
  Colours <- HiCBricks:::Make_colours(palette = palette, extrapolate_on = extrapolate_on, 
                                      direction = col_direction)
  two.sample <- (length(Bricks) == 2)
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
    y.coord.breaks <- y.coord.breaks - min(y.coord.breaks)
    x.coord.breaks <- x.coord.breaks - min(x.coord.breaks)
    if (length(Bricks) == 2) {
      Upper.tri.map <- Matrix.df[Matrix.df$dist >= 0, ]
      Lower.tri.map <- Matrix.df[Matrix.df$dist <= 0, ]
      Lower.tri.map$dist <- abs(Lower.tri.map$dist)
      Upper.rotated.map <- HiCBricks:::RotateHeatmap(Matrix = Upper.tri.map, 
                                                     value.var = "rescale", upper = TRUE)
      Lower.rotated.map <- HiCBricks:::RotateHeatmap(Matrix = Lower.tri.map, 
                                                     value.var = "rescale", upper = FALSE)
      Entire.rotated.map <- rbind(Upper.rotated.map, Lower.rotated.map)
      y.coord.breaks <- y.coord.breaks/2
      y.coord.breaks <- c(rev(y.coord.breaks) * -1, y.coord.breaks)
      y_axis.coord.labs <- c(rev(y_axis.coord.labs), y_axis.coord.labs)
    }
    else {
      Upper.tri.map <- Matrix.df[Matrix.df$dist >= 0, ]
      Entire.rotated.map <- HiCBricks:::RotateHeatmap(Matrix = Upper.tri.map, 
                                                      value.var = "rescale", upper = TRUE)
      y.coord.breaks <- y.coord.breaks/2
    }
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
  if (!is.null(tad_ranges)) {
    Boundaries.obj <- HiCBricks:::Format_boundaries_normal_heatmap(Bricks = Bricks, 
                                                                   resolution = resolution, Ranges = tad_ranges, group_col = group_col, 
                                                                   cut_corners = cut_corners, colour.col = tad_colour_col, 
                                                                   colours = colours, colours_names = colours_names, 
                                                                   region.chr = x.coord.parsed["chr"], region.start = as.numeric(x.coord.parsed["start"]), 
                                                                   region.end = as.numeric(x.coord.parsed["end"]), distance = distance, 
                                                                   rotate = rotate)
  }
  if (rotate) {
    ids <- xcoords <- ycoords <- NULL
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
  if (!is.null(tad_ranges)) {
    line.group <- x <- y <- NULL
    ThePlot <- ThePlot + geom_line(data = Boundaries.obj, 
                                   aes(x = x, y = y, group = line.group, colour = colours), 
                                   size = line_width)
    ThePlot <- ThePlot + scale_colour_manual(values = colours)
  }
  
  save(xlims, file="xlims.Rdata", version=2)
  save(x_axis.coord.labs, file="x_axis.coord.labs.Rdata", version=2)
  save(x.coord.breaks, file="x.coord.breaks", version=2)
  
  ThePlot <- ThePlot + scale_x_continuous(limits = xlims, expand = c(0, 
                                                                     0), breaks = x.coord.breaks, labels = x_axis.coord.labs)
  ThePlot <- ThePlot + scale_y_continuous(limits = ylims, expand = c(0, 
                                                                     0), breaks = y.coord.breaks, labels = y_axis.coord.labs)
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