#' @name hapDistribution
#' @title display geography distribution
#' @importFrom maps map
#' @importFrom graphics polygon
#' @import mapdata
#' @usage
#' hapDistribution(hap, AccINFO, Acc.col, LON.col, LAT.col, hapNames,
#'                 mapRegion = "world", map.package = "maps",
#'                 zColours = zColours,
#'                 legend = "leftbottom", symbolSize = 1,
#'                 maxZVal = 1, oceanCol="white",
#'                 landCol="grey90", borderCol = "black", cex.legend = 0.8,
#'                 ...)
#' @param hap an object of `hapResult` class
#' @param AccINFO a data.frame contains accession information
#' @param Acc.col,LON.col,LAT.col column names of accession(`Acc.col`),
#' longitude(`LON.col`) and latitude(`LAT.col`)
#' @param hapNames haplotype names used for display
#' @param legend a keyword specified the position of legend, one of
#' "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center"
#' @param map.package which map packages should be used, one of `maps` and `rworldmap`
#' @param cex.legend character expansion factor for legend relative to current par("cex")
#' @param mapRegion if `map.package` set as "rworldmap", a country name from
#' `rworldmap::getMap()[['NAME']]` or 'world','africa','oceania','eurasia','uk'
#'  sets map extents, overrides `xlim`,`ylim`; if `map.package` set as 'maps',
#'
#' @inheritParams rworldmap::mapPies
#' @export
hapDistribution <-
    function(hap,
             AccINFO,
             Acc.col,
             LON.col,
             LAT.col,
             hapNames,
             mapRegion = "world",
             map.package = "maps",
             zColours = zColours,
             legend = "leftbottom",
             symbolSize = 1,
             maxZVal = 1,
             oceanCol = "white",
             landCol = "grey90",
             borderCol = "black",
             cex.legend = 0.8,
             ...) {
        hap <- na.omit(hap)
        if (missing(zColours))
            zColours <- rainbow(length(hapNames))

        if (missing(hapNames))
            hapNames <- hap$Hap

        if (inherits(hap, "hapResult")) {
            acc2hap <- hap[, "Hap"]
            names(acc2hap) = hap[, "Accession"]
        } else if (inherits(hap, "hapSummary")) {
            acc2hap <- rep(hap[, "Hap"], hap[, "freq"])
            names(acc2hap) <- unlist(strsplit(hap[, "Accession"], ";"))
        }

        # extract geo data
        AccINFO$Hap <- acc2hap[AccINFO[, Acc.col]]
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
            ) #对数据进行预处理
        nameX = LON.col
        nameY = LAT.col
        ratio <- 1
        symbolScale <- 1
        maps::map(mapRegion, ...)
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
                    pointsInCircle * (cumulatProps[sliceNum +
                                                       1] - cumulatProps[sliceNum])
                )))
                P <-
                    list(
                        x = ratio * radius * cos(
                            2 * pi * seq(cumulatProps[sliceNum],
                                         cumulatProps[sliceNum + 1], length = n)
                        ) + dF[locationNum,
                               nameX],
                        y = radius * sin(
                            2 * pi * seq(cumulatProps[sliceNum],
                                         cumulatProps[sliceNum + 1], length = n)
                        ) + dF[locationNum,
                               nameY]
                    )
                graphics::polygon(c(P$x, dF[locationNum, nameX]), c(P$y, dF[locationNum,
                                                                            nameY]), col = zColours[sliceNum])
            }
        }


        # add legend
        if (legend[1] != FALSE) {
            if (legend[1] == TRUE)
                legend = "leftbottom"
            if(length(legend) == 2)

            legend(legend,
                   legend = hapNames,
                   fill = zColours,
                   cex = cex.legend)
        }
    }
