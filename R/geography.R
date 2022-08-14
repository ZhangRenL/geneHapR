#' @name hapDistribution
#' @title Display of Geography Distribution
#' @description
#' show distribution of intereted haplotypes on maps
#' @importFrom maps map
#' @importFrom graphics polygon
#' @usage
#' hapDistribution(hap, AccINFO, LON.col, LAT.col, hapNames,
#'                 database = "world", regions = ".",
#'                 zColours = zColours,
#'                 legend = TRUE, symbolSize = 1,
#'                 ratio = 1, cex.legend = 0.8,
#'                 ...)
#' @examples
#' \donttest{
#' data("geneHapR_test")
#' hapDistribution(hapResult,
#'                 AccINFO = AccINFO,
#'                 LON.col = "longitude",
#'                 LAT.col = "latitude",
#'                 hapNames = c("H001", "H002", "H003"))
#' }
#' @param hap an object of hapResult class
#' @param AccINFO a data.frame contains accession information
#' @param LON.col,LAT.col column names of
#' longitude(`LON.col`) and latitude(`LAT.col`)
#' @param hapNames haplotype names used for display
#' @param zColours colours to apply to the pie section for each attribute column
#' @param symbolSize a numeric specified the symbol size
#' @param legend a keyword specified the position of legend, one of
#' "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center";
#' or a numeric vector of length two contains x,y coordinate of the legend
#' @param cex.legend character expansion factor for legend relative to current `par("cex")`
#' @param ratio the ratio of Y to N in the output map, set to 1 as default
#' @inheritParams maps::map
#' @return No return value
#' @export
hapDistribution <-
    function(hap,
             AccINFO,
             LON.col,
             LAT.col,
             hapNames,
             database = "world",
             regions = ".",
             zColours = zColours,
             legend = TRUE,
             symbolSize = 1,
             ratio = 1,
             cex.legend = 0.8,
             ...) {
        hap <- na.omit(hap)
        if (missing(zColours))
            zColours <- rainbow(length(hapNames))

        if (missing(hapNames))
            hapNames <- hap$Hap

        if (inherits(hap, 'hapResult')) {
            acc2hap <- hap[, "Hap"]
            names(acc2hap) = hap[, "Accession"]
        } else if (inherits(hap, "hapSummary")) {
            acc2hap <- rep(hap[, "Hap"], hap[, "freq"])
            names(acc2hap) <- unlist(strsplit(hap[, "Accession"], ";"))
        }

        # extract geo data
        AccINFO$Hap <- acc2hap[row.names(AccINFO)]
        geoData <- AccINFO[, c("Hap", LON.col, LAT.col)]
        geoData$value <- 1

        # reshape the data
        formu <- paste0(LON.col, "+", LAT.col, "~Hap")
        dF <-
            reshape2::dcast(
                geoData,
                formula = formu,
                value.var = "value",
                fun.aggregate = length
            )


        symbolScale <- 1
        maps::map(database = database, regions = regions, ...)
        for (locationNum in 1:length(dF[, hapNames[1]])) {
            sliceValues <- as.numeric(dF[locationNum, hapNames])
            if (sum(sliceValues, na.rm = TRUE) == 0)
                next
            cumulatProps <-
                c(0,
                  cumsum(sliceValues) / sum(sliceValues,
                                            na.rm = TRUE))
            pointsInCircle = 360
            radius <-
                sqrt(sum(sliceValues, na.rm = TRUE)) * symbolScale
            radius <- radius * symbolSize
            for (sliceNum in 1:length(sliceValues)) {
                n <- max(2, floor((
                    pointsInCircle *
                        (cumulatProps[sliceNum + 1] -
                             cumulatProps[sliceNum])
                )))
                P <-
                    list(
                        x = ratio * radius * cos(
                            2 * pi * seq(cumulatProps[sliceNum],
                                         cumulatProps[sliceNum + 1], length = n)
                        ) + dF[locationNum,
                               LON.col],
                        y = radius * sin(
                            2 * pi * seq(cumulatProps[sliceNum],
                                         cumulatProps[sliceNum + 1], length = n)
                        ) + dF[locationNum,
                               LAT.col]
                    )
                graphics::polygon(c(P$x, dF[locationNum, LON.col]),
                                  c(P$y, dF[locationNum, LAT.col]),
                                  col = zColours[sliceNum])
            }
        }


        # add legend
        if (legend[1] != FALSE) {
            if (legend[1] == TRUE)
                legend = "bottomleft"

            if(length(legend) == 1) {
                legend(legend,
                       legend = hapNames,
                       fill = zColours,
                       cex = cex.legend)
            } else if(length(legend) == 2) {
                if(!is.numeric(legend))
                    stop("'legend' should a TRUE/FALSE indicate wherther plot the legend;
or keword indicate where to plot the legend;
or a numeric vector of length two contains x,y coordinate of the legend.")
                legend(legend[1],
                       legend[2],
                       legend = hapNames,
                       fill = zColours,
                       cex = cex.legend)
            }

        }
    }
