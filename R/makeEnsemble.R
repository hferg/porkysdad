#### makeEnsemble
# A series of functions to fit models, build an ensemble and return
# the ensemble model and it's evaulation.
# Arguments: Presence points.
#     Predictor layers stack.
#     Extent.
#     Projection.
#     Vector of models to fit
#     Buffer radius.

#' functions included
#' makeEnsemble
#' futurePresence
#' evaluateEnsemble

#' futurePresence
#' This function takes a suitability surface for the future of some species
#' (either the output of an ensemble model, or a single model prediction) and
#' then uses that surface to weight the probability of draws of presences. In
#' this way, presence points for the future can be obtained, and then models
#' that rely on presence points can be run (e.g. niche overlap models).
#' @param suitability The suitability surface that weights the probability of
#' presence draws.
#' @param npres The number of presence points to sample
#' @param nreps The number of replicates of future presence sampling to draw.

futurePresence <- function(suitability, npres, nreps) {
  # first get the values out of the suitability - that's the weights.
  # then make a vector of cell numbers?
  # then sample npres from that vector weighted by the values.
  probs <- x <- raster::getValues(suitability)
  # Add in a 0 when NA.
  probs[is.na(probs)] <- 0
  probs[probs < 0.5] <- 0
  futures <- vector(mode = "list", length = nreps)
  for (i in seq_len(nreps)) {
    pres <- sample(1:length(probs), npres, prob = (probs / max(probs)), replace = TRUE)
    # Now pres is the number of cells sampled as present...
    x[!is.na(x)] <- 0
    x[pres] <- 1
    futures[[i]] <- setValues(suitability, x)
  }
}


####################################################################################
# TESTING

# Set working directory for Laura's future projections.

setwd("/home/hfg/Documents/students/kuurne_laura/sa_bees_future/sa_species_futures")
root <- getwd()

# Then the different year and RCP scenarios are in seperate folders.
# e.g.

setwd(file.path(root, "RCP262050/Distribution"))
suitability <- x <- readRDS(list.files()[1])
npres <- 500
nreps <- 12