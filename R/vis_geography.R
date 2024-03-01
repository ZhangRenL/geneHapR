#' @name hapDistribution
#' @title Display of Geography Distribution
#' @description
#' show distribution of intereted haplotypes on maps
#' @importFrom maps map
#' @importFrom graphics polygon
#' @examples
#' \donttest{
#' data("geneHapR_test")
#' data(china)
#' hapDistribution(hapResult,
#'                 AccINFO = AccINFO,
#'                 LON.col = "longitude",
#'                 LAT.col = "latitude",
#'                 hapNames = c("H001", "H002", "H003"))
#' }
#' @param hap an object of hapResult class
#' @param lty.pie the line type of pie border
#' @param AccINFO a data.frame contains accession information
#' @param LON.col,LAT.col column names of
#' longitude(`LON.col`) and latitude(`LAT.col`)
#' @param hapNames haplotype names used for display
#' @param hap.color,zColours colors to apply to the pie section for each attribute column,
#' "zColours" will be detached in future.
#' @param symbolSize a numeric specified the symbol size.
#' It will be detached in future. Please use "symbol.lim" instead.
#' @param symbol.lim a numeric vector give the maximum and minimum size  of each symbol
#' @param legend a keyword specified the position of legend, one of
#' "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center";
#' or a numeric vector of length two contains x,y coordinate of the legend
#' @param cex.legend character expansion factor for legend relative to current `par("cex")`
#' @param ratio the ratio of Y to N in the output map, set to 1 as default
#' @param lwd.pie line width of the pies
#' @param showlabel a bool vector indicates whether show the labels
#'  which represens number of individuals. Default as TRUE.
#' @param borderCol.pie The color of pie's border, default is NA, which means no border will be plotted
#' @inheritParams maps::map
#' @return No return value
#' @export
#' @param label.cex a number indicates the text size in label, default as 0.8
#' @param label.col color of the labels, default as "black"
#' @param label.font Font of label, 1 for normal, 2 for bold, 3 for italica, 4 for bold-italica
#' @param label.adj the position of label, default as c(0.5, 0.5)
#' @param map.fill.color vector of colors. If fill is FALSE, the first color is
#'  used for plotting all lines, and any other colors are ignored. Otherwise,
#'  the colors are matched one-one with the polygons that get selected by the
#'  region argument (and are reused cyclically, if necessary). If fill = TRUE,
#'  the default boundary line colour is given by par("fg"). To change this,
#'  you can use the border argument (see '...'). A color of NA causes the
#'  corresponding region to be deleted from the list of polygons to be drawn.
#'  Polygon colors are assigned after polygons are deleted due to values of the
#'   xlim and ylim arguments
hapDistribution <- function(hap,
                            AccINFO,
                            LON.col,
                            LAT.col,
                            hapNames,
                            database = "world",
                            regions = ".",
                            hap.color = hap.color,
                            zColours = zColours,
                            legend = TRUE,
                            symbolSize = 1,
                            symbol.lim = c(1,10),
                            ratio = 1,
                            cex.legend = 0.8,
                            lwd.pie = 1,
                            borderCol.pie = NA,
                            lty.pie = 1,
                            showlabel = TRUE,
                            label.col = "black",
                            label.cex = 0.8,
                            label.font = 1,
                            label.adj = c(0.5,0.5),
                            map.fill.color = 1,
                            ...) {
    hap <- na.omit(hap)
    if(! missing(hap.color)){
        zColours <- hap.color
    } else {
        if (missing(zColours))
            zColours <- rainbow(length(hapNames))
    }

    if(is.function(zColours))
        zColours <- zColours(length(hapNames))

    if (missing(hapNames))
        hapNames <- hap$Hap

    if (inherits(hap, 'hapResult')) {
        acc2hap <- hap[, "Hap"]
        names(acc2hap) <- hap[, "Accession"]
    } else if (inherits(hap, "hapSummary")) {
        acc2hap <- rep(hap[, "Hap"], hap[, "freq"])
        names(acc2hap) <- unlist(strsplit(hap[, "Accession"], ";"))
    }

    # extract geo data
    AccINFO$Hap <- acc2hap[row.names(AccINFO)]
    geoData <- AccINFO[, c("Hap", LON.col, LAT.col)]
    geoData$value <- 1
    # Check the column data format
    if(! inherits(geoData[,LAT.col], "numeric")){
        warning("The '", LAT.col, "' column should be numeric")
        geoData[,LAT.col] <- suppressWarnings(as.numeric(geoData[,LAT.col]))
        if(! inherits(geoData[,LAT.col], "numeric"))
            stop("We couldn't conver the '", LAT.col,
                 "' column as 'numeric'")
    }
    if(! inherits(geoData[,LON.col], "numeric")){
        warning("The '", LON.col, "' column should be numeric")
        geoData[,LON.col] <- suppressWarnings(as.numeric(geoData[,LON.col]))
        if(! inherits(geoData[,LON.col], "numeric"))
            stop("We couldn't conver the '", LON.col,
                 "' column as 'numeric'")
    }


    # reshape the data
    formu <- paste0(LON.col, "+", LAT.col, "~Hap")
    dF <-
        reshape2::dcast(
            geoData,
            formula = formu,
            value.var = "value",
            fun.aggregate = length
        )
    # Check wether all hapnames in Accinfo
    if( ! all(hapNames %in% names(dF))) {
        n_names <- hapNames[! hapNames %in% names(dF)]
        n_names <- paste0(n_names, collapse = ",", sep = "")
        m <- paste0("We could not find any individual belong to '",
                  n_names, "' in your AccINFO, Procceed?")
        proceed <- askYesNo(m, default = FALSE)
        if(proceed) {
            hapNames <- hapNames[hapNames %in% names(dF)]
        } else return(NULL)
    }
    if (any(tolower(database) == "china", tolower(regions) == "china" )){
        # usethis::use_data(china, internal = F)
        requireNamespace("sf")
        data("china")
        plot(china, border = map.fill.color,...)
    } else maps::map(database = database, regions = regions, col = map.fill.color, ...)
    # maps::map(database = database, regions = regions)
    # scale the circle size using symbol.lim instead of symbolSize
    symbolSize <- 1
    freqs <- apply(dF[,hapNames], 1, sum)
    mx <- sqrt(max(freqs))
    mn <- sqrt(min(freqs[freqs > 0]))
    dif <- mx - mn

    # plotting
    for (locationNum in 1:nrow(dF)) {
        sliceValues <- as.numeric(dF[locationNum, hapNames])
        if (sum(sliceValues, na.rm = TRUE) == 0)
            next

        cumulatProps <-
            c(0,
              cumsum(sliceValues) / sum(sliceValues,
                                        na.rm = TRUE))
        pointsInCircle <- 360
        radius <-
            sqrt(sum(sliceValues, na.rm = TRUE)) * symbolSize
        # scale the circle size
        radius <- (radius - mn) / dif * abs(diff(symbol.lim)) + min(symbol.lim)
        for (sliceNum in seq_len(length(sliceValues))) {
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
                              col = zColours[sliceNum], lty = lty.pie,
                              lwd = lwd.pie, border = borderCol.pie)
        }
        if(isTRUE(showlabel))
            text(dF[locationNum, LON.col],
                 dF[locationNum, LAT.col],
                 sum(dF[locationNum, hapNames]),
                 col = label.col,
                 cex = label.cex,
                 font = label.font,
                 adj = label.adj)

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

