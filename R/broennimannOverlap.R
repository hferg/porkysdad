
#' broennimannEcospat
#' 
#' Takes raster layers and occurence data and fits the broennimann
#' niche similarity, equivalency and overlap models.
#' This is simply a wrapper to functions in the R package ecospat (https://cran.r-project.org/web/packages/ecospat/index.html)
#' @param native_stack A raster stack of the environmental data for the 
#' range of the native species.
#' @param invasive_stack A raster stack of the environmental data for the
#' range being invaded
#' @param natice_occ Occurence points for the native species (i.e. where the 
#' species that will be invading into the invasive range is in it's native
#' range). A two-column matrix with column headings x and y (corresponding to
#' longitude and lattitude)
#' @param invasice_occ Occurence points for the species in the invasive range
#' (either the invading species in it's invasive range, OR the species that 
#' might compete with the invader).A two-column matrix with column headings 
#' x and y (corresponding to longitude and lattitude)
#' @param kept_axes The nf parameter in the function dudi.pca - is the number
#' of kept axes - defaults to 2.
#' @name broennimannEcospat
#' @export

broennimannEcospat <- function(native_stack, invasive_stack, 
  native_occ, invasive_occ, kept_axes = 2) {
  # first make the tables to go into the vignette stuff - 
  # it will be x, y, all the environmental stuff, and then 
  # an occurence yes/no column. Everywhere it ISN'T is a no.

  # make sure projections match.
  if (sp::proj4string(native_stack) != sp::proj4string(invasive_stack)) {
    stop("Projection of native and invasive stacks do not match.")
  }

  if (ncol(native_occ) != 2) {
    stop("Native occurences must be just two columns, x (long) and y (lat).")
    if (paste(colnames(native_occ), collapse = "") != "xy") {
      stop("Column names of native occurences must be x (long) and y (lat)")
    }
  }

  if (ncol(invasive_occ) != 2) {
    stop("Invasive occurences must be just two columns, x (long) and y (lat).")
    if (paste(colnames(invasive_occ), collapse = "") != "xy") {
      stop("Column names of invasive occurences must be x (long) and y (lat)")
    }
  }

  sp::coordinates(native_occ) <- ~x+y
  sp::coordinates(invasive_occ) <- ~x+y

  # turn points into spatial object and match projection to stacks.
  latlonproj <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
  sp::proj4string(native_occ) <- latlonproj
  sp::proj4string(invasive_occ) <- latlonproj

  native_occ <- sp::spTransform(native_occ, CRS = sp::proj4string(native_stack))
  invasive_occ <- sp::spTransform(invasive_occ, CRS = sp::proj4string(invasive_stack))

  # Add in the occurence points to these tables.
  native_tmp <- native_stack[[1]]
  native_tmp[!is.na(native_tmp)] <- 0
  native_tmp[raster::extract(native_tmp, native_occ, cellnumbers = TRUE)[ , "cells"]] <- 1
  names(native_tmp) <- "occurence"
  invasive_tmp <- invasive_stack[[1]]
  invasive_tmp[!is.na(invasive_tmp)] <- 0
  invasive_tmp[raster::extract(invasive_tmp, invasive_occ, cellnumbers = TRUE)[ , "cells"]] <- 1
  names(invasive_tmp) <- "occurence"

  native_stack <- raster::stack(native_stack, native_tmp)  
  invasive_stack <- raster::stack(invasive_stack, invasive_tmp)  

  # get all values from the native and raster stacks.
  print("Extracting environmental data...")
  nat <- raster::rasterToPoints(native_stack)
  inv <- raster::rasterToPoints(invasive_stack)

  # remove NA
  print("Removing NAs...")
  nat <- nat[!is.na(nat[ , 3]), ]
  inv <- inv[!is.na(inv[ , 3]), ]


  # Compute the environmental PCA before adding species occurence.
  print("Calculating PCA...")
  pca.env <- ade4::dudi.pca(
    rbind(nat, inv)[ , 3:(ncol(nat) - 1)], 
    scannf = FALSE, 
    nf = kept_axes
    )

  # Now add to this species occurence.
  # predict scores on axes - this is where the occurence comes in.
  scores.globclim <- pca.env$li

  # scores for the species in native and invaded ranges.
  scores.sp.nat <- ade4::suprow(pca.env, nat[which(nat[ , ncol(nat)] == 1), 3:(ncol(nat) - 1)])$li
  scores.sp.inv <- ade4::suprow(pca.env, inv[which(inv[ , ncol(inv)] == 1), 3:(ncol(nat) - 1)])$li

  # and now scores for the entire native and study area.
  scores.clim.nat <- ade4::suprow(pca.env, nat[ , 3:(ncol(nat) - 1)])$li
  scores.clim.inv <- ade4::suprow(pca.env, inv[ , 3:(ncol(inv) - 1)])$li

  # from this, calculate the occurence density grid (an ecospat function)
  print("Calculating density grids...")
  grid.clim.nat <- ecospat::ecospat.grid.clim.dyn(
    glob = scores.globclim,
    glob1 = scores.clim.nat,
    sp = scores.sp.nat,
    R = 100, 
    th.sp = 0
    )

  grid.clim.inv <- ecospat::ecospat.grid.clim.dyn(
    glob = scores.globclim,
    glob1 = scores.clim.inv,
    sp = scores.sp.inv,
    R = 100, 
    th.sp = 0
    )

  # calculate niche overlap
  print("Calculating niche overlap...")
  D.overlap <- ecospat::ecospat.niche.overlap(grid.clim.nat, grid.clim.inv, cor = T)$D

  # nice equivalency
  print("Calculating niche equivalency...")
  eq.test <- ecospat::ecospat.niche.equivalency.test(
    grid.clim.nat,
    grid.clim.inv,
    rep = 1000,
    alternative = "greater"
    )

  # niche similarit test
  print("Calculating niche similarity...")
  sim.test <- ecospat::ecospat.niche.similarity.test(
    grid.clim.nat,
    grid.clim.inv,
    rep = 1000,
    alternative = "greater"
    )

  scores <- list(
    global = scores.globclim,
    native = scores.clim.nat,
    invasive = scores.clim.inv
    )

  grids <- list(
    native = grid.clim.nat,
    invasive = grid.clim.inv
    )

  res <- list(
    overlap = D.overlap,
    eq.test = eq.test,
    sim.test = sim.test,
    pca = pca.env,
    scores = scores,
    grids = grids
    )

  return(res)
}

# # I need to find the part in the broenniman functions where the two
# # species, and locations, are seperated from each other - what if one
# # species is in cliamte 1, and the other climate 2?

# # Start with climate data for each othe study sites. What are these data?
# # One trait, or many?
# # climate data is x, y (coords), then the values are columns. Data for all
# # sites of the study area (i.e. all pixels).
# # first read in some climate data, and also just get the bee data (while I am
# # # at it...)
climate_root <- "/home/hfg/Documents/projects/advent/data/bioclim"
bioclim_now <- raster::stack(paste0(climate_root, "/bioclim_now_cropped.grd"))
layers <- c("bio_1", "bio_4", "bio_5", "bio_13", "bio_15")
bioc <- bioclim_now[[which(names(bioclim_now) %in% layers)]]
native_stack <- invasive_stack <- bioc

# # since these bees are at the same place clim a and clim b are the same.

bees <- read.csv("/home/hfg/Documents/projects/advent/data/d1c/step_bumblebees/CANPOLIN_2014_05_13_ungrided.csv")
native_occ <- bees[bees$taxon == "Bombus lucorum", 4:3]
invasive_occ <- bees[bees$taxon == "Bombus monticola", 4:3]
colnames(native_occ) <- colnames(invasive_occ) <- c("x", "y")

# x <- broennimannEcospat(native_stack, invasive_stack, native_occ, invasive_occ)

# # # Start with the functions that are provided in the ESM
