% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getGenePOS.R
\name{getGeneRanges}
\alias{getGeneRanges}
\title{Get Gene Ranges}
\usage{
getGeneRanges(gff= gff,
             geneID = geneID,
             type = type,
             gffTermContaingeneID = "Parent")
}
\arguments{
\item{gff}{imported gff}

\item{geneID}{target geneID}

\item{type}{vector consist with one or more types in gff}

\item{gffTermContaingeneID}{which term contains the geneID in your gff,
defalt is Parent}
}
\value{
GRanges
}
\description{
Get Gene Ranges
}
\examples{
data("geneHapR_test")
geneRanges <- getGeneRanges(gff = gff,
                            geneID = "test1G0387",
                            type = "CDS",
                            gffTermContaingeneID = "Parent")
}
