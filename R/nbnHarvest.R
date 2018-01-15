################################################################################
#
# A series of functions that are designed to collect occurence records from the
# NBN Atlas.

################################################################################
#' isDataset
#'
#' Takes a number, and probes the URL with that number as the extension to see
#' if that number leads to a dataset, or not. This basically just works off the
#' fact that NBN records are all indexed at URLs with a varying number at the
#' end.
#' @name isDataset
#' @param num The number to test.
#' @export

isDataset <- function(num) {

  dset_root <- "https://registry.nbnatlas.org/public/showDataResource/"
  main <- xml2::read_html(paste0(dset_root, "dr", num))

  hd <- rvest::html_nodes(main, "h1")
  typ <- rvest::html_nodes(main, ".section h3")
  hd_tx <- rvest::html_text(hd)
  typ_tx <- rvest::html_text(typ)

  if (hd_tx == "United Kingdom's collections & data partners") {
    res <- "Not a dataset"
  } else {
    # This tests if the dataset is a proper dataset, or a species list.
    if ("Species lists" %in% typ_tx) {
      res <- "Species list"
    } else {
      res <- hd_tx
    }

  }
  return(res)
}

################################################################################
#' findDatasets
#'
#' A shortcut function to apply isDataset over a large range of potential
#' dataset numbers - returns a 2D array containing the numbers of actual
#' datasets (i.e. excluding species lists and things that aren't anything).
#' @name findDatasets
#' @param nums A vector of numbers to test.
#' @export

findDatasets <- function(nums) {
  xx <- data.frame(num = nums, name = NA)
  xx$name <- sapply(xx$num, function(x) isDataset(x))
  xx <- xx[xx$name != "Not a dataset", ]
  xx <- xx[xx$name != "Species list", ]
  return(xx)
}

################################################################################
#' getLinks
#'
#' Collects links to individual records from valid dataset IDs. It does this
#' by collecting the number of records, breaking that up into chunks of 1000,
#' and then going to the URL that displays links to individual records for the
#' first 1-1000 records, 1001-2000 records etc. and then saves those links,
#' writing them to the working directory. It writes them in chunks of 1000. I
#' write these links to the disk since a) it can be a super slow process and I
#' don't want to repeat it if something goes wrong, and b) I don't want to
#' repeat it over and over and over-query the servers at NBN. If a reference
#' list of species names is provided then only species records that match that
#' list will be saved.
#' @name getLinks
#' @param num A valid dataset number.
#' @param reference_list A list of species names that links should be collected
#' for.
#' @param print If TRUE then the dataset number if printed to screen - this is
#' mostly just for testing...
#' @param filename This is the prefix of the filename that the saved links are
#' written to.
#' @export

getLinks <- function(num, reference_list, print = FALSE, filename) {
  if (print) {
    print(num)
  }
  recs_hub_url <- paste0("https://records.nbnatlas.org/occurrence/search?q=data_resource_uid:dr", num)

  num_recs <- xml2::read_html(recs_hub_url) %>%
    rvest::html_nodes("#returnedText strong") %>%
    rvest::html_text()

  num_recs <- as.numeric(gsub(",", "", num_recs))

  breaks <- data.frame(offset = c(0:floor(num_recs / 1000)), max = 1000)
  for (i in seq_len(nrow(breaks))) {
    offset <- breaks$offset[i]
    max <- breaks$max[i]
    recs_url <- paste0("https://records.nbnatlas.org/occurrences/search?taxa=&q=data_resource_uid%3Adr",
                       num, "&fq=&wkt=&lat=&lon=&radius=&offset=",
                       offset, "&max=",
                       max)

    main <- xml2::read_html(recs_url)

    links <- main %>%
      rvest::html_nodes(".occurrenceLink") %>%
      rvest::html_attr("href")

    species <- main %>%
      rvest::html_nodes("#results i") %>%
      rvest::html_text()

    if (!is.null(reference_list)) {
      links <- links[tolower(species) %in% tolower(reference_list)]
    }

    ifelse (i == 1, h <- TRUE, h <- FALSE)
    write.table(links, file = filename, sep = ",", append = !h, col.names = h, row.names = FALSE)
  }
}

################################################################################
#' getLinksAllRecords
#' I just noticed that the link https://records.nbnatlas.org/occurrence/search
#' seems to display all the records in the NBN - it seems like this might just
#' be a better way to go about this - move through all the records collecting
#' links for the species on the C4 list, and then download them. Since this
#' misses out the step where I find out the dataset information I think that
#' I need to make sure that I include this when downloading the actual data,
#' but that's another function to go hand-in hand with this one.

getLinksAllRecords <- function(range, filename, max) {
  recs_hub_url <- paste0("https://records.nbnatlas.org/occurrence/search")

  num_recs <- xml2::read_html(recs_hub_url) %>%
    rvest::html_nodes("#returnedText strong") %>%
    rvest::html_text()

  num_recs <- as.numeric(gsub(",", "", num_recs))
  breaks <- data.frame(offset = c(0:floor(num_recs / max)), max = max)
  breaks$offset <- breaks$offset * breaks$max

  for (i in min(range):max(range)) {
    offset <- breaks$offset[i]
    max <- breaks$max[i]
    recs_url <- paste0("https://records.nbnatlas.org/occurrences/search?taxa=&q=&fq=&wkt=&lat=&lon=&radius=&dir=&sort=&offset=",
                       offset, "&max=",
                       max)

    main <- xml2::read_html(recs_url)

    links <- main %>%
      rvest::html_nodes(".occurrenceLink") %>%
      rvest::html_attr("href")

    # Have to take the whole row, since some of these don't have a species - 
    # just a genus or class or whatever...
    info <- main %>%
      rvest::html_nodes(".rowA") %>%
      rvest::html_text()

    # Extract species from this, if present.
    getSpecies <- function(info, reference_list) {
      if (grepl("species:", info)) {
        x <- strsplit(info, "species:")[[1]][2]
        x <- strsplit(gsub("\n", "", x), "\\|")[[1]][1]
        res <- stringr::str_trim(x)
      } else {
        res <- FALSE
      }
      return(res)
    }

    species <- sapply(info, getSpecies)
    # grepl doesn't work here.
    if (!is.null(reference_list)) {
      links <- links[tolower(species) %in% tolower(reference_list)]
    }

    ifelse (i == 1, h <- TRUE, h <- FALSE)
    write.table(links, file = filename, sep = ",", append = !h, col.names = h, row.names = FALSE)
  }
}


################################################################################
#' getRecordData
#'
#' Takes a set of links, and downloads the data for the individual records at
#' each of those links, saving that data to disk link by link. This is a little
#' slower than just collecting it all then saving it, but it means that a) we
#' don't have to store millions of records worth of data in memory, and b) that
#' if for whatever reason the connection is lost, we have all the data up to
#' that point in time. It also has a built-in sleep every 1000 links (approx.
#' half an hour) so that it stops querying the server for around ten minutes.
#' @name getRecordData
#' @param links Either a vector of links, or the filename of some links file
#' that was written out by \link[porkysdad]{getLinks}.
#' @param filename The name of the file that the data is going to be written
#' into.
#' @export

getRecordData <- function(links, filename) {
  root_url <- "https://records.nbnatlas.org"
  for (i in seq_along(links)) {
    progress(i / length(links) * 100)
    main <- xml2::read_html(paste0(root_url, links[i]))
    ds <-  main %>%
      rvest::html_nodes("#datasetTable") %>%
      rvest::html_table()
    ev <- main %>%
      rvest::html_nodes("#eventTable") %>%
      rvest::html_table()
    tx <- main %>%
      rvest::html_nodes("#taxonomyTable") %>%
      rvest::html_table()
    gs <- main %>%
      rvest::html_nodes("#geospatialTable") %>%
      rvest::html_table()
    ds <- ds[[1]]
    ev <- ev[[1]]
    tx <- tx[[1]]
    gs <- gs[[1]]

    # now compile the stuff I want to keep and return a transposed table
    # of just that stuff...
    d <- data.frame(
      binomial = ifelse(any(tx$X1 == "Scientific name"),
                        tx[tx$X1 == "Scientific name", 2], NA),
      kingdom = ifelse(any(tx$X1 == "Kingdom"),
                       tx[tx$X1 == "Kingdom", 2], NA),
      phylum = ifelse(any(tx$X1 == "Phylum"),
                      tx[tx$X1 == "Phylum", 2], NA),
      class = ifelse(any(tx$X1 == "Class"),
                     tx[tx$X1 == "Class", 2], NA),
      order = ifelse(any(tx$X1 == "Order"),
                     tx[tx$X1 == "Order", 2], NA),
      family = ifelse(any(tx$X1 == "Family"),
                      tx[tx$X1 == "Family", 2], NA),
      genus = ifelse(any(tx$X1 == "Genus"),
                     tx[tx$X1 == "Genus", 2], NA),
      species = ifelse(any(tx$X1 == "Species"),
                       tx[tx$X1 == "Species", 2], NA),
      authority = ifelse(any(tx$X1 == "Name according to"),
                         tx[tx$X1 == "Name according to", 2], NA),
      lat = ifelse(any(gs$X1 == "Latitude"),
                   gs[gs$X1 == "Latitude", 2], NA),
      lon = ifelse(any(gs$X1 == "Longitude"),
                   gs[gs$X1 == "Longitude", 2], NA),
      proj = ifelse(any(gs$X1 == "Geodetic datum"),
                    gs[gs$X1 == "Geodetic datum", 2], NA),
      error = ifelse(any(gs$X1 == "Coordinate uncertainty in metres"),
                     gs[gs$X1 == "Coordinate uncertainty in metres", 2], NA),
      date = ifelse(any(ev$X1 == "Record date"),
                    strsplit(ev[ev$X1 == "Record date", 2], "\\n")[[1]][1], NA),
      provider = ifelse(any(ds$X1 == "Data provider"),
                        ds[ds$X1 == "Data provider", 2], NA),
      nbn_id = ifelse(any(ds$X1 == "Occurrence ID"),
                      ds[ds$X1 == "Occurrence ID", 2], NA),
      license = ifelse(any(ds$X1 == "License"),
                       ds[ds$X1 == "License", 2], NA),
      dataset_id = ifelse(any(ds$X1 == "Dataset id"),
                          ds[ds$X1 == "Dataset id", 2], NA),
      collection_code = ifelse(any(ds$X1 == "Collection Code"),
                               ds[ds$X1 == "Collection Code", 2], NA)
    )
    ifelse (i == 1, h <- TRUE, h <- FALSE)
    write.table(d, file = filename, sep = ",", append = !h, col.names = h, row.names = FALSE)
    # then if i is a factor of something or other take a break for five minutes.
    # This should take a ten minute break every half hour.

    if (i %% 3000 == 0) {
      Sys.sleep(rnorm(1, mean = 600, sd = 15))
    }

  }
}

################################################################################
#' getAllData
#'
#' Wraps all the above functions together - give it a range of numbers and it
#' will check if they are datasets, get the links then download the data.
#' Can set to "links only" mode, which will download the links, but not get the
#' data. To go hand in hand with this, it can also be supplied with some
#' links, in which case it won't bother testing to see if a dataset is a dataset
#' and so on...
#' @name getAllData
#' @param numrange A range of numbers to test to see if they are datasets, and
#' get data from.
#' @param filename The name of the file that the data is going to be written
#' into - prefixes _links.csv for the links file.
#' @param reference_list A list of species to get data for.
#' @param links If specified then data is downloaded from these links.
#' @param links_only If TRUE then no data is downloaded, but the links are
#' gathered.
#' @export

getAllData <- function(numrange, filename, reference_list = NULL, links = NULL,
                       links_only = FALSE, ...) {
  if (is.null(links)) {
    print("Checking datasets...")
    x <- findDatasets(numrange)
    print("Finding links...")
    l_file <- paste0(filename, "_links.csv")
    lapply(x$num, function(x) getLinks(x, reference_list = reference_list,
                                            filename = l_file, ...))
  } else {
    links <- links
  }

  if (!links_only) {
    if (is.null(links)) {
      links <- read.csv(paste0(filename, "_links.csv"))
    } else {
      links <- links
    }
    print("Gathering data...")
    getRecordData(links, filename)
    print("Done.")
  }
}
