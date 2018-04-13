#' getOverlap
#' Calculate the number of cells and the area in km2 of overlap between map1
#' and map2.
#' @param map1 a presence/absence (1/0) raster.
#' @param map2 a presence/absence (1/0) raster.

getOverlap <- function(map1, map2) {
  # Find the cells of presence in map2 that are also present in
  # map1.

  # Add the rasters - this will produce a map of minimum common extent.
  # It will also give three values - 0 where no species are present, 1
  # where only one species is present, and 2 where both species are present.
  combo <- suppressWarnings(map1 + map2)
  # crop the map1 and map2 maps to smallest common extent.
  map1 <- crop(map1, combo)
  map2 <- crop(map2, combo)


  # get cells that equal 2 - these are where overlap occurs.
  overlap_cells <- raster::xyFromCell(combo, which(raster::getValues(combo) == 2), spatial = FALSE)
  map1_cells <- raster::xyFromCell(map1, which(raster::getValues(map1) == 1), spatial = FALSE)
  map2_cells <- raster::xyFromCell(map2, which(raster::getValues(map2) == 1), spatial = FALSE)

  # calculate the km2 area of each cell in the combined map.
  combo_area <- area(combo)

  overlap_area <- sum(raster::extract(combo_area, overlap_cells))
  map1_area <- sum(raster::extract(combo_area, map1_cells))
  map2_area <- sum(raster::extract(combo_area, map2_cells))

  # create results object.
  res <- list(map1_ncells = nrow(map1_cells),
              map2_ncells = nrow(map2_cells),
              overlap_ncells = nrow(overlap_cells),
              map1_area = map1_area,
              map2_area = map2_area,
              overlap_area = overlap_area)
  return(res)
}

plotOverlap <- function(host, cuckoo) {
  host[host == 1] <- 2
  combo <- suppressWarnings(host + cuckoo)
  breakpoints <- c(-1, 0, 1, 2, 3)
  colours <- c("lightgrey", "#FC8D59", "#FFFFBF", "#91BFDB")
  plot(combo, breaks = breakpoints, col = colours)
}
