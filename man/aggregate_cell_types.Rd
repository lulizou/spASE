% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CSIDE_utils.R
\name{aggregate_cell_types}
\alias{aggregate_cell_types}
\title{Aggregates the pixel occurrences for each cell type in the \code{\linkS4class{RCTD}} object}
\usage{
aggregate_cell_types(myRCTD, barcodes, doublet_mode = T)
}
\arguments{
\item{barcodes}{the barcodes, or pixel names, of the \code{\linkS4class{SpatialRNA}} object to be used when counting cell types.}

\item{doublet_mode}{(default TRUE) if TRUE, uses RCTD doublet mode weights. Otherwise, uses RCTD full mode weights}

\item{RCTD}{an \code{\linkS4class{RCTD}} object with annotated cell types e.g. from the \code{\link{run.RCTD}} function.}
}
\value{
a named vector of number of pixel occurrences for each cell type
}
\description{
The difference with \code{\link{count_cell_types}} is that this function does not filter out pixels
based on total cell type weight, as occurs in the CSIDE algorithm.
}
