% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis_geography.R
\name{hapDistribution}
\alias{hapDistribution}
\title{Display of Geography Distribution}
\usage{
hapDistribution(hap, AccINFO, LON.col, LAT.col, hapNames,
                database = "world", regions = ".",
                zColours = zColours,
                legend = TRUE, symbolSize = 1,
                ratio = 1, cex.legend = 0.8,
                lwd.pie = 1,
                ...)
}
\arguments{
\item{hap}{an object of hapResult class}

\item{AccINFO}{a data.frame contains accession information}

\item{LON.col, LAT.col}{column names of
longitude(\code{LON.col}) and latitude(\code{LAT.col})}

\item{hapNames}{haplotype names used for display}

\item{database}{
    character string naming a geographical database, a list of
    \code{x}, \code{y}, and \code{names} obtained from a previous call to
    \code{map} or a spatial object of class \code{SpatialPolygons} or \code{SpatialLines}.
    The string choices include a \code{\link[maps]{world}} map,
    three USA databases (\code{\link[maps]{usa}}, \code{\link[maps]{state}},
    \code{\link[maps]{county}}), and more (type \code{help(package='maps')} to see the package index). If the requied database is in a different package that has not been attached, the string may be started with "packagename::". The
    location of the map databases may be overridden by setting the
    \code{R_MAP_DATA_DIR} environment variable.
  }

\item{regions}{
    character vector that names the polygons to draw.
    Each database is composed of a collection of polygons, and each polygon has
    a unique name.
    When a region is composed of more than one polygon, the individual polygons
    have the name of the region, followed by a colon and a qualifier,
    as in \code{michigan:north} and \code{michigan:south}.
    Each element of \code{regions} is matched
    against the polygon names in the database and, according to
    \code{exact}, a subset is selected for drawing. The regions may also be 
    defined using (perl) regular expressions. This makes it possible to use 
    'negative' expressions like \code{"Norway(?!:Svalbard)"}, which means Norway
    and all islands except Svalbard. All entries are case insensitive.
    The default selects all polygons in the database.
  }

\item{zColours}{colours to apply to the pie section for each attribute column}

\item{legend}{a keyword specified the position of legend, one of
"bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center";
or a numeric vector of length two contains x,y coordinate of the legend}

\item{symbolSize}{a numeric specified the symbol size}

\item{ratio}{the ratio of Y to N in the output map, set to 1 as default}

\item{cex.legend}{character expansion factor for legend relative to current \code{par("cex")}}

\item{lwd.pie}{line width of the pies}

\item{...}{Extra arguments passed to \code{polygon} or \code{lines}. Of particular interest may be the options \code{border} and\code{lty} that control the color and line type of the polygon borders when \code{fill = TRUE}.}
}
\value{
No return value
}
\description{
show distribution of intereted haplotypes on maps
}
\examples{
\donttest{
data("geneHapR_test")
hapDistribution(hapResult,
                AccINFO = AccINFO,
                LON.col = "longitude",
                LAT.col = "latitude",
                hapNames = c("H001", "H002", "H003"))
}
}
