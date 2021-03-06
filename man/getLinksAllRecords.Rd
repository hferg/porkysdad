% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nbnHarvest.R
\name{getLinksAllRecords}
\alias{getLinksAllRecords}
\title{getLinksAllRecords
I just noticed that the link https://records.nbnatlas.org/occurrence/search
seems to display all the records in the NBN - it seems like this might just
be a better way to go about this - move through all the records collecting
links for the species on the C4 list, and then download them. Since this
misses out the step where I find out the dataset information I think that
I need to make sure that I include this when downloading the actual data,
but that's another function to go hand-in hand with this one.}
\usage{
getLinksAllRecords(range, batch_size, filename, reference_list)
}
\arguments{
\item{range}{The number of batches to get links for (i.e. range = 2:7 will 
get the second to seventh batches).}

\item{batch_size}{The size of each batch of records to look at.}

\item{filename}{This is the prefix of the filename that the saved links are
written to.}

\item{reference_list}{A list of species names that links should be collected
for. Function compares the species in a record to the reference list, and 
retains it if there is a match. Case doesn't matter... Should be a single 
vector of binomials.}

\item{print}{If TRUE then the dataset number if printed to screen - this is
mostly just for testing...}
}
\description{
getLinksAllRecords
I just noticed that the link https://records.nbnatlas.org/occurrence/search
seems to display all the records in the NBN - it seems like this might just
be a better way to go about this - move through all the records collecting
links for the species on the C4 list, and then download them. Since this
misses out the step where I find out the dataset information I think that
I need to make sure that I include this when downloading the actual data,
but that's another function to go hand-in hand with this one.
}
