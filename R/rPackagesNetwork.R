# Scraping CRAN for package info to make a 
# social network of R packages.

# First get the names of all the packages.
library(magrittr)

getPackageInfo <- function(pack_name) {
  pack <- xml2::read_html(paste0("https://cran.ma.imperial.ac.uk/web/packages/", pack_name, "/index.html")) %>%
    rvest::html_nodes("table:nth-child(3) td") %>%
    rvest::html_text()
  return(pack)
}

getMassPackage <- function(packages, wait_secs) {
  res <- vector(mode = "list", length = length(packages))
  for (i in seq_along(packages)) {
    print(i)
    res[[i]] <- getPackageInfo(packages[i])
    Sys.sleep(wait_secs)
  }
  
  return(res)
}

parsePackageInfo <- function(pack) {
  res <- vector(mode = "list", length = 5)
  names(res) <- c("depends", "imports", "suggests", "version", "published")
  if ("Depends:" %in% pack) {
    deps <- pack[which(pack == "Depends:") + 1] %>%
      strsplit(", ")
    if (length(deps) > 1) {
      res$depends <- deps[[1]][2:length(deps[[1]])]
    }
  }

  if ("Suggests:" %in% pack) {
    sugs <- pack[which(pack == "Suggests:") + 1] %>%
      strsplit(", ")
    res$suggests <- sugs[[1]]
  }

  if ("Imports:" %in% pack) {
    imps <- pack[which(pack == "Imports:") + 1] %>%
      strsplit(", ")
    res$imports <- imps[[1]]
  }

  res$version <- pack[which(pack == "Version:") + 1]
  res$published <- pack[which(pack == "Published:") + 1]
  return(res)
}

massParsePackage <- function(pack_info) {
  res <- lapply(pack_info, parsePackageInfo)
  return(res)
}

makeNetworkTable <- function(pack_info) {

}

packs <- xml2::read_html("https://cran.ma.imperial.ac.uk/web/packages/available_packages_by_name.html") %>% 
  rvest::html_nodes("td a") %>%
  rvest::html_text() %>%
  getMassPackage(wait_secs = 0)
pack_parsed <- massParsePackage(packs)
saveRDS(pack_parsed, file = "/home/hfg/pCloudDrive/hfg_projects/rPackages/packageInfo.RDS")





# How do I put the resulting information into something useful for
# a network plot?
# I think I want a matrix with a row and column for each package,
# and then I want to put something in the cell if the row imports
# or depends on the column...
# This will be a bit of a pain with DiagrammeR - making the edge_df is
# going to be awkward...

# this times out after 100 queries, so build in a wait-time...


x <- lapply(packs, getPackageInfo)


xx <- lapply(x, parsePackageInfo)

t <- getPackageInfo(packs[[1]]) %>%
  parsePackageInfo()



