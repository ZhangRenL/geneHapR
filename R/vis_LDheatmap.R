# ldheatmap - Plots measures of pairwise linkage disequilibria for SNPs
# Copyright (C) 2004 J.Shin, S. Blay, N. Lewin-Koh, J.Graham, B.McNeney

# To cite LDheatmap in publications use:
# Shin J-H, Blay S, McNeney B and Graham J (2006). LDheatmap: An R
# Function for Graphical Display of Pairwise Linkage Disequilibria
# Between Single Nucleotide Polymorphisms. J Stat Soft, 16 Code Snippet 3

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA


###########################################################################
#' @name LDheatmap
#' @aliases LDheatmap
#' @title This function produces a pairwise LD plot.
#' @description \code{LDheatmap()} is used to produce a graphical display, as a heat map,
#' of pairwise linkage disequilibrium (LD) measurements for SNPs.
#' The heat map is a false color image in the upper-left diagonal of a square plot.
#' Optionally, a line parallel to the diagonal of the image indicating
#' the physical or genetic map positions of the SNPs may be added, along
#' with text reporting the total length of the genomic region considered.
#'
#' @param hap R object of hapSummary or hapResult class.
#' @param gff annotations
#' @param Chr,start,end,geneID chromosome, start and end position, gene ID for
#' extract annotation in target range.
#' @param distances A character string to specify whether the provided map locations
#' are in physical or genetic distances.
#' If \code{distances = "physical"} (default), the text describing the total
#' length of the region will be \dQuote{Physical Length:XXkb} where XX is the
#' length of the region in kilobases. If \code{distances = "genetic"}, the
#' text will be \dQuote{Genetic Map Length:YYcM} where YY is
#' the length of the region in centiMorgans.
#' If \code{gdat} is an object of class LDheatmap,
#' \code{distances} is taken from \code{gdat}.
#'
#' @param LDmeasure A character string specifying the measure of LD
#' - either allelic correlation \eqn{r^2} or Lewontin's
#' |D\eqn{'}|; default = \code{"r"} for \eqn{r^2};
#' type \code{"D'"} for |D\eqn{'}|. This argument is ignored when the user has already
#' supplied calculated LD measurements through \code{gdat} (i.e., when \code{gdat}
#' is a matrix of pairwise LD measurements or an object of class \code{"LDheatmap"}).
#'
#' @param title A character string for the main title of the plot.
#' Default is \dQuote{Pairwise LD}.
#' @param add.map If \code{TRUE} (default) a diagonal line indicating
#' the physical or genetic map positions of the SNPs will be added to
#' the plot, along with text indicating the total length of the
#' genetic region.
#' @param map.height the height of gene map, default is 0.02
#' @param colorLegend If \code{TRUE} (default) the color legend is drawn.
#' @param geneMapLocation A numeric value specifying the position of the line
#' parallel to the diagonal of the matrix; the larger the value, the
#' farther it lies from the matrix diagonal. Ignored when \code{add.map = FALSE}.
#' @param geneMapLabelX A numeric value specifying the x-coordinate
#' of the text indicating the total length of the genomic region
#' being considered. Ignored when \code{add.map = FALSE}.
#' @param geneMapLabelY A numeric value specifying the y-coordinate
#' of the text indicating the total length of the genomic region
#' being considered. Ignored when \code{add.map = FALSE}.
#' @param color A range of colors to be used for drawing the heat map. Default
#' is \code{grDevices::colorRampPalette(c("red", "grey"))(30)}.
#' @param color_gmodel,color_snp color of genemodel and snp, default as grey80.
#' @param newpage If \code{TRUE} (default), the heat map will be drawn on a new page.
#' @param name A character string specifying the name of the LDheatmap
#' graphical object (\code{grob}) to be produced.
#' @param vp.name A character string specifying the name of the viewport
#' where the heat map is going to be drawn.
#' @param pop If \code{TRUE}, the viewport where the heat map is drawn is
#' \code{pop}ped (i.e. removed) from the viewport tree after drawing.
#' Default = \code{FALSE}.
#' @param text If \code{TRUE}, the LD measurements are printed on each cell.
#' @importFrom genetics is.genotype LD nallele as.genotype
#' @details
#' The input object \code{gdat} can be a data frame of \code{genotype} objects
#' (a data structure from the \pkg{genetics} package), a \code{SnpMatrix} object
#' (a data structure from the \pkg{snpStats} package), or
#' any square matrix with values between 0 and 1 inclusive.
#' LD computation is much faster for \code{SnpMatrix} objects than for
#' \code{genotype} objects.
#' In the case of a matrix of LD values between 0 and 1,
#' the values above the diagonal will be plotted.
#' In the display of LD, SNPs appear in the order supplied by the
#' user as the horizontal and vertical coordinates are increased
#' and one moves along the off-diagonal line, from the bottom-left
#' to the top-right corner. To achieve this, the conventions of
#' the \code{image()} function have been adopted, in which horizontal
#' coordinates correspond to the rows of the matrix and vertical coordinates
#' correspond to columns, and vertical coordinates are indexed in increasing
#' order from bottom to top.
#' For the argument \code{color}, an appropriate
#' color palette for quantitative data is recommended,
#' as outlined in the help page of
#' the \code{\link[RColorBrewer:ColorBrewer]{brewer.pal}()} function of
#' the
#' \pkg{RColorBrewer} package.
#' See the package vignette \code{LDheatmap} for more examples and details
#' of the implementation. Examples of adding ``tracks'' of genomic
#' annotation above a flipped heatmap are in the package vignette
#' \code{addTracks}.
#'
#'
#' @return An object of class \code{"LDheatmap"} which contains the following components:
#' \item{LDmatrix}{ The matrix of pairwise LD measurements plotted in the heat map. }
#' \item{LDheatmapGrob}{ A grid graphical object (grob) representing the produced heat map. }
#' \item{heatmapVP}{ The viewport in which the heat map is drawn. See \link[grid:viewport]{viewport}.}
#' \item{genetic.distances}{The vector of the supplied physical or
#' genetic map locations, or the vector of equispaced marker distances
#' when no distance vector is supplied.}
#' \item{distances}{ A character string specifying whether the provided map
#' distances are physical or genetic. }
#' \item{color}{ The range of colors used for drawing the heat map. }
#' The \code{grob} \code{LDheatmapGrob} has three \code{grob}s as its children (components).
#' They are listed below along with their own children and respectively represent
#' the color image with main title, genetic map and color key of the heat map:
#' \code{"heatMap"} - \code{"heatmap"}, \code{"title"};
#' \code{"geneMap"} - \code{"diagonal"}, \code{"segments"},
#' \code{"title"}, \code{"symbols"}, \code{"SNPnames"}; and
#' \code{"Key"} - \code{"colorKey"}, \code{"title"}, \code{"labels"},
#' \code{"ticks"}, \code{"box"}.
#'
#' @references Shin J-H, Blay S, McNeney B and Graham J (2006). LDheatmap:
#' An R Function for Graphical Display of Pairwise Linkage
#' Disequilibria Between Single Nucleotide Polymorphisms.
#' Journal of Statistical Software, \bold{16} Code Snippet 3
#'
#' @note The produced heat map can be modified in two ways.
#' First, it is possible to edit \emph{interactively} the grob components of the heat map,
#' by using the function \code{\link[grid:grid.edit]{grid.edit}};
#' the function will not work if there is no
#' open graphical device showing the heat map.
#' Alternatively, the user can use the function
#' \code{\link[grid:grid.edit]{editGrob}} and work with
#' the grob \code{LDheatmapGrob} returned by \code{LDheatmap}.
#' See Examples for usage.

#' \code{LDheatmap()} uses \code{\link[grid:Grid]{Grid}}, which
#' does not respond to \code{par()} settings.
#' Hence modifying \code{par()} settings of \code{mfrow} and \code{mfcol}
#' will not work with \code{LDheatmap()}. The Examples section shows how to
#' display multiple heat maps on one plot without the use
#' of \code{par()}.
#'
#'
#' @examples # Pass LDheatmap a SnpMatrix object
#' data(geneHapR_test)
#' plot_LDheatmap(hap = hapResult,
#'                gff = gff,
#'                Chr = hapResult[1,2],
#'                start = 4000, end = 8200)
#' @export
plot_LDheatmap <- function (hap,gff,Chr,start,end,geneID,
                            distances = "physical",
                            LDmeasure = "r",
                            title = "Pairwise LD",
                            add.map = TRUE,
                            map.height = 1,
                            colorLegend = TRUE,
                            geneMapLocation = 0.15,
                            geneMapLabelX = NULL,
                            geneMapLabelY = NULL,
                            color = NULL,
                            color_gmodel="grey",
                            color_snp = "grey",
                            newpage = TRUE,
                            name = "ldheatmap",
                            vp.name = NULL,
                            pop = FALSE,
                            text = FALSE){
    # TO remove flip option
    flip = TRUE
    if(inherits(hap, "hapResult")){
        gdat <- data.frame(hap)
        SNP.name <- t(hap[4,-c(1, ncol(hap))])[,1]
        genetic.distances <- suppressWarnings(as.numeric(gdat[2,])) %>%
            na.omit()
        gdat <- gdat[-c(1:4), -c(1,ncol(gdat))]
        for(i in seq_len(ncol(gdat))){
            gdat[,i] <- sapply(gdat[,i],
                           function(x) {
                               if(stringr::str_detect(x,"[|]")) {
                                   stringr::str_replace(x,"[|]","/")
                               } else {
                                   paste0(x,"/",x)
                               }
                           }) %>% unlist()
            gdat[,i] <- genetics::as.genotype(gdat[,i])
        }


    } else return("Use hapResult instead of hapSummary")

    LDheatmap(gdat,gff=gff,Chr=Chr,start=start,end=end,geneID=geneID,
              genetic.distances = genetic.distances,
              distances = distances,
              LDmeasure = LDmeasure,
              title = title,
              add.map = add.map,map.height = map.height,
              colorLegend = colorLegend,
              geneMapLocation = geneMapLocation,
              geneMapLabelX = geneMapLabelX,
              geneMapLabelY = geneMapLabelY,
              SNP.name = SNP.name,
              color = color,
              color_gmodel = color_gmodel,
              color_snp = color_snp,
              newpage = newpage,
              name = "ldheatmap",
              vp.name = vp.name,
              pop = pop,
              flip = flip,
              text = text)
}


#' @import grid
#' @importFrom graphics plot.new
LDheatmap <- function (gdat,gff,Chr,start,end,geneID,
                       genetic.distances = NULL,
                       distances = "physical",
                       LDmeasure = LDmeasure,
                       title = "Pairwise LD",
                       add.map = TRUE, map.height = 0.02,
                       colorLegend = TRUE,
                       geneMapLocation = 0.15,
                       geneMapLabelX = NULL,
                       geneMapLabelY = NULL,
                       SNP.name = NULL,
                       color = NULL,
                       color_gmodel = color_gmodel,
                       color_snp = color_snp,
                       newpage = TRUE,
                       name = "ldheatmap",
                       vp.name = NULL,
                       pop = FALSE,
                       flip = TRUE,
                       text = FALSE)
{
    # requireNamespace("grid")

    # Draw the Color Key
    if (is.null(color)) {
        color <- grDevices::colorRampPalette(c("red", "grey"))(30)
    } else if(length(color) < 10) {
        color <- grDevices::colorRampPalette(c(color))(30)
    }


    # adds a genetic map to the heatmap plot along the diagonal
    # This part only identifies if the genemap will need to be flipped, does not do anything else yet
    if (is.null(flip)) {
        flip <- FALSE
    }

    ## Calculate or extract LDmatrix, then stored in LDMatrix as an upper triangular matrix
    # using a data.frame
    if (inherits(gdat, "data.frame")) {
        for (i in 1:ncol(gdat)) {
            if (!genetics::is.genotype(gdat[, i]))
                stop("column ", i, " is not a genotype object\n")
        }

        ## Exclude SNPs with more or less than 2 alleles:
        gvars <-
            unlist(sapply(gdat, function(x)
                genetics::nallele(x) == 2))
        if(any(! gvars))
            warning("Only bi-alleles supported,",
                    "Variables with less or more than 2 allels will be omitted.")

        if(sum(gvars) < 2)
            stop("Variants number is less than two after removed non-bialleles")
        genetic.distances <- genetic.distances[gvars]
        gdat <- gdat[gvars]

        ## Sort data in ascending order of SNPs map position:
        if (!is.vector(genetic.distances))
        {
            stop("Distance should be in the form of a vector")
        }
        o <- order(genetic.distances)
        genetic.distances <- genetic.distances[o]
        gdat <- gdat[, o]
        myLD <- genetics::LD(gdat)
        if (LDmeasure == "r")
            LDmatrix <- myLD[[LDmeasure]] ^ 2
        else if (LDmeasure == "D'")
            LDmatrix <- abs(myLD[[LDmeasure]])
        else
            stop("Invalid LD measurement, choose r or D'.")
    }


    ## Draw the heat map
    heatmapVP <-
        grid::viewport(
            width = unit(.8, "snpc"),
            height = unit(.8, "snpc"),
            name = vp.name
        )
    flipVP <-
        grid::viewport(
            width = unit(.8, "snpc"),
            height = unit(.8, "snpc"),
            y = 0.6,
            angle = -45,
            name = "flipVP"
        )

    # Colour selection scaling
    if (color[1] == "blueToRed")
        color = rainbow(20,
                        start = 4 / 6,
                        end = 0,
                        s = .7)[20:1]
    if (newpage)
        grid::grid.newpage()
    mybreak <- 0:length(color) / length(color)

    imgLDmatrix <- LDmatrix

    # Flip or not, determines way data is read into the display
    byrow <- ifelse(flip, FALSE, TRUE) #FALSE if flip = TRUE

    colcut <-
        as.character(cut(
            1 - imgLDmatrix,
            mybreak,
            labels = as.character(color),
            include.lowest = TRUE
        ))


    # Determines if colour is done as an integer or as a colour code, updates accordingly
    if (is.numeric(color))
        colcut <- as.integer(colcut)
    ImageRect <-
        makeImageRect(dim(LDmatrix)[1], dim(LDmatrix)[2], colcut, name = "heatmap", byrow)


    # Controls text placement
    ImageText <- NULL
    if (text)
        ImageText <-
        makeImageText(dim(LDmatrix)[1],
                      dim(LDmatrix)[2],
                      round(imgLDmatrix, digits = 2),
                      name = "heatmaptext")
    title <-
        grid::textGrob(title, 0.5, 1.05, gp = grid::gpar(cex = 1.0), name = "title")

    if (flip) {
        ImageRect <- grid::editGrob(ImageRect, vp = flipVP)
        if (text) {
            # Added flip = TRUE parameter to better utilize makeImageText() in the flipped case
            ImageText <-
                makeImageText(
                    dim(LDmatrix)[1],
                    dim(LDmatrix)[2],
                    round(imgLDmatrix, digits = 2),
                    name = "heatmaptext",
                    flip = TRUE
                )
            textVal <- ImageText
            ImageText <-
                grid::editGrob(
                    ImageText,
                    vp = flipVP,
                    rot = 0,
                    just = c("right", "top")
                )
        }
    }

    # Updates heatmap in the gTree
    heatMap <-
        grid::gTree(children = grid::gList(ImageRect, ImageText, title),
              name = "heatMap")

    # Draw a diagonal line indicating the physical or genetic map positions of the SNPs
    nsnps <- ncol(LDmatrix)
    step <- 1 / (nsnps - 1)
    ind <- match(SNP.name, row.names(LDmatrix), nomatch = 0)
    geneMapVP <- NULL
    if (flip)
        geneMapVP <- flipVP
    geneMap <-
        LDheatmapMapNew.add (
            nsnps,
            genetic.distances = genetic.distances,
            geneMapLocation = geneMapLocation,
            add.map=add.map,
            geneMapLabelX = geneMapLabelX,
            geneMapLabelY = geneMapLabelY,
            distances = distances,
            vp = geneMapVP,
            SNP.name = SNP.name,
            ind = ind,
            color_gmodel = color_gmodel,color_snp = color_snp,
            flip = flip,
            gff=gff,Chr=Chr,
            start=start,end=end,geneID=geneID,map.height = map.height
        )

    # Draw the Color Key
    if (colorLegend)
        Key <- LDheatmapLegend.add(color, LDmeasure, heatmapVP)
    else
        Key <- NULL

    # Assemble the heatmap, genetic map and color key into a grob and draw it
    LDheatmapGrob <- grid::gTree(
        children = grid::gList(heatMap, geneMap, Key),
        vp = heatmapVP,
        name = name,
        cl = "ldheatmap"
    )
    grid::grid.draw(LDheatmapGrob)
    if (pop) {
        grid::downViewport(heatmapVP$name)
        grid::popViewport()
    } #pop the heat map viewport

    ldheatmap <-
        list(
            LDmatrix = LDmatrix,
            LDheatmapGrob = LDheatmapGrob,
            heatmapVP = heatmapVP,
            flipVP = geneMapVP,
            genetic.distances = genetic.distances,
            distances = distances,
            color = color
        )
    class(ldheatmap) <- "LDheatmap"
    invisible(ldheatmap)
} # function LDheatmap ends



preDrawDetails.ldheatmap <- function(x) {
    fontsize <- grid::convertX(unit(1 / 20, "grobwidth", grid::rectGrob()), "points")
    grid::pushViewport(grid::viewport(gp = grid::gpar(fontsize = fontsize)))
}


postDrawDetails.ldheatmap <- function(x) {
    grid::popViewport()
}


preDrawDetails.symbols <- function(x) {
    fontsize <- grid::convertX(unit(1 / 20, "grobwidth", grid::rectGrob()), "points")
    grid::pushViewport(grid::viewport(gp = grid::gpar(fontsize = fontsize)))
}

postDrawDetails.symbols <- function(x) {
    grid::popViewport()
}


# Functions below are used in conjunction with LDHeatmap.R
# Provide supporting functionality for major function

makeImageRect <- function(nrow, ncol, cols, name, byrow = TRUE) {
    xx <- (1:ncol) / ncol
    yy <- (1:nrow) / nrow
    # Creates coordinate pairs, repeating either column numbers (if byrow = TRUE) or row numbers (if byrow = FALSE) to force that type of fill
    if (byrow) {
        right <- rep(xx, nrow)
        top <- rep(yy, each = ncol)
    } else {
        right <- rep(xx, each = nrow)
        top <- rep(yy, ncol)
    }
    grid::rectGrob(
        x = right,
        y = top,
        width = 1 / ncol,
        height = 1 / nrow,
        just = c("right", "top"),
        gp = grid::gpar(col = NA, fill = cols),
        name = name
    )
}

makeImageText <- function(nrow, ncol, cols, name, flip = FALSE) {
    cols <- as.character(cols)
    cols[is.na(cols)] <- ""
    cols <- paste(" ", cols)
    xx <- (1:ncol) / ncol
    yy <- (1:nrow) / nrow

    # Need to fill cells in different order, as was done to generate image
    if (flip) {
        right <- rep(xx, each = nrow)
        top <- rep(yy, ncol)
    }
    else{
        right <- rep(xx, nrow)
        top <- rep(yy, each = ncol)
    }
    grid::textGrob(
        cols,
        x = right,
        y = top,
        gp = grid::gpar(cex = 0.3),
        just = c("right", "top"),
        name = name
    )
}

LDheatmapLegend.add <- function(color, LDmeasure, vp) {
    ImageRect <-
        makeImageRect(2, length(color), cols = c(rep(NA, length(color)), color[length(color):1]),
                      "colorKey")
    keyVP <-
        grid::viewport(
            x = 1.1,
            y = -.10,
            height = .10,
            width = .5,
            just = c("right", "bottom"),
            name = "keyVP"
        )
    #Adding the label 'Color key'
    if (LDmeasure == "r") {
        ttt <- expression(paste(R ^ 2, " Color Key"))
    } else {
        ttt <- "D' Color Key"
    }
    title <-
        grid::textGrob(
            ttt,
            x = 0.5,
            y = 1.25,
            name = "title",
            gp = grid::gpar(cex = 0.7)
        )

    #Adding labels to the color key
    labels <-
        grid::textGrob(
            paste(0.2 * 0:5),
            x = 0.2 * 0:5,
            y = 0.25,
            gp = grid::gpar(cex = 0.6),
            name = "labels"
        )

    #Drawing ticks at the bottom axis of the color key
    ticks <-
        grid::segmentsGrob(
            x0 = c(0:5) * 0.2 ,
            y0 = rep(0.4, 6),
            x1 = c(0:5) * 0.2 ,
            y1 = rep(0.5, 6),
            name = "ticks"
        )

    #Drawing a box around the color key
    box <-
        grid::linesGrob(
            x = c(0, 0, 1, 1, 0),
            y = c(0.5, 1, 1, 0.5, 0.5),
            name = "box"
        )

    key <-
        grid::gTree(
            children = grid::gList(ImageRect, title, labels, ticks, box),
            name = "Key",
            vp = keyVP
        )
    key
}

LDheatmapMapNew.add <- function(nsnps,
                                add.map,map.height = 0.02,
                                genetic.distances,
                                geneMapLocation = 0.15,
                                geneMapLabelX = NULL,
                                geneMapLabelY = NULL,
                                distances = "physical",
                                vp = NULL,
                                SNP.name = NULL,
                                ind = 0,
                                color_gmodel = color_gmodel,
                                color_snp = color_snp,
                                flip = FALSE,
                                gff,Chr,start,end,geneID) {
    snp <- ((1:nsnps - 1) + 0.5) / nsnps
    #####################
    if (add.map) {
        if(missing(start) | missing(end)){
            min.dist <- min(genetic.distances)
            max.dist <- max(genetic.distances)
        }else{
            min.dist <- min(genetic.distances, start, end)
            max.dist <- max(genetic.distances, start, end)
        }
        total.dist <- max.dist - min.dist

        if (flip)
            geneMapLocation <-
            (-geneMapLocation) # geneMapLocation gets flipped, reflects the gene bar

        # Drawing the diagonal line
        i <- 0.01
        seq.x <- c(0.5 * geneMapLocation + 1 / (nsnps * 2),
                   1 + 0.5 * geneMapLocation - 1 / (nsnps * 2))
        seq.y <- c(-0.5 * geneMapLocation + 1 / (nsnps * 2),
                   1 - 0.5 * geneMapLocation - 1 / (nsnps * 2))
        diagonal <-
            grid::linesGrob(
                seq.x,
                seq.y,
                gp = grid::gpar(lty = 1,col=color_snp),
                name = "diagonal",
                vp = vp
            ) # we draw the line with linesGrob, based on geneMapLocation seq

        xs <- ys <- c()
        w <- map.height
        if(!missing(Chr) | !missing(gff)){
            gene <- GenomicRanges::GRanges(Chr,
                                           IRanges::IRanges(
                                               start = min.dist,
                                               end = max.dist))
            gff <- gff[gff %over% gene]
            type <- tolower(gff$type)
            probe1 <- stringr::str_detect(type, "cds") | stringr::str_detect(type, "utr")
            gff <- gff[probe1]
            if(!missing(geneID)){
                ids <- tolower(gff$ID)
                nms <- tolower(gff$Name)
                geneID <- tolower(geneID)
                probe2 <- stringr::str_detect(nms, geneID) | stringr::str_detect(ids, geneID)
                gff <- gff[probe2]
            }
            gff$Parent <- unlist(gff$Parent)
            strand <- unique(gff@strand)[1]
            ps <- unique(gff$Parent)
            for(p in seq_len(length(ps))){
                gffp <- gff[gff$Parent == ps[p]]
                for(i in seq_len(length(gffp@ranges))){
                    start <- gffp@ranges[i]@start
                    wid <- gffp@ranges[i]@width
                    end <- start + wid
                    r1 <- (start - min.dist) / total.dist
                    r2 <- (end - min.dist) / total.dist
                    x1 <- seq.x[1] + (seq.x[2] - seq.x[1]) * r1
                    x2 <- seq.x[1] + (seq.x[2] - seq.x[1]) * r2
                    if(stringr::str_detect(tolower(gff$type[i]),"cds"))
                        v <- 1 else v <- 1/2
                    xsp <- c(x1-w*v,x1+w*v,x2+w*v,x2-w*v) - w * (p * 3 - 1.5)
                    ysp <- c(x1+w*v,x1-w*v,x2-w*v,x2+w*v) + w * (p * 3 - 1.5)
                    xs <- c(xs, xsp, NA)
                    ys <- c(ys, ysp, NA)
                }
                xsp <- c(seq.x[1]-0.002,seq.x[1]+0.002,seq.x[2]+0.002,seq.x[2]-0.002)
                ysp <- c(seq.x[1]+0.002,seq.x[1]-0.002,seq.x[2]-0.002,seq.x[2]+0.002)
                xs <- c(xs,NA,xsp - w * (p * 3 - 1.5),NA)
                ys <- c(ys,NA,ysp + w * (p * 3 - 1.5),NA)
            }
            gmodel <- grid::polygonGrob(x = xs,
                                        y = ys + 0.15,
                                        gp = gpar(fill = color_gmodel, lty = 0), vp = vp)
        } else  gmodel <- NULL



        ## Adding line segments to the plot: (point1 <- > point2)
        ## point1: relative position of a SNP on the scaled line
        ## point2: position of that SNP on the LD image
        regionx <- seq.x[1] +
            ((genetic.distances - min.dist) / total.dist) * (seq.x[2] -
                                                                 seq.x[1])
        regiony <- seq.y[1] +
            ((genetic.distances - min.dist) / total.dist) * (seq.y[2] -
                                                                 seq.y[1])
        segments <-
            grid::segmentsGrob(snp,
                         snp,
                         regionx,
                         regiony,
                         name = "segments",
                         vp = vp,
                         gp = gpar(col = color_snp))

        ## Adding the text indicating Physical length of the region under study
        if (distances == "physical")
            mapLabel <-
            paste("Physical Length:", round((total.dist / 1000), 1),
                  "kb", sep = "")
        else
            mapLabel <-
            paste("Genetic Map Length:", round(total.dist, 1), "cM", sep = "")


        if(strand == "+") mapLabel <- paste0(mapLabel, "\n5' -> 3'") else
            mapLabel <- paste0(mapLabel, "\n3' <- 5'")
        if (!flip) {
            if (is.null(geneMapLabelY))
                geneMapLabelY <- 0.3
            if (is.null(geneMapLabelX))
                geneMapLabelX <- 0.5
        }
        else {
            if (is.null(geneMapLabelY))
                geneMapLabelY <- 0.8
            if (is.null(geneMapLabelX))
                geneMapLabelX <- 0.4
        }
        title <- grid::textGrob(
            mapLabel,
            geneMapLabelX - w * length(ps) * 4,
            geneMapLabelY + w * length(ps) * 4,
            gp = grid::gpar(cex = 0.9),
            just = "left",
            name = "title"
        )

        geneMap <-
            grid::gTree(children = grid::gList(diagonal,gmodel,segments, title),
                  name = "geneMap")

        ## Labelling some SNPs
        if (!is.null(SNP.name) && (any(ind!= 0))) {
            if (flip) {
                length_SNP_name <- max(nchar(SNP.name))
                long_SNP_name <-
                    paste(rep(8, length_SNP_name), collapse = "")
                name_gap <-
                    grid::convertWidth(grid::grobWidth(grid::textGrob(long_SNP_name)), "npc", valueOnly = TRUE) /
                    sqrt(2)
                diagonal <-
                    grid::linesGrob(
                        seq.x,
                        seq.y,
                        gp = grid::gpar(lty = 1),
                        name = "diagonal",
                        vp = vp
                    )
                #diagonal <- grid::linesGrob(seq.x+name_gap, seq.y-name_gap, gp = grid::gpar(lty = 1), name = "diagonal", vp = vp)
                segments <-
                    grid::segmentsGrob(snp,
                                 snp,
                                 regionx,
                                 regiony,
                                 name = "segments",
                                 vp = vp)
                #segments <- segmentsGrob(snp+name_gap, snp-name_gap, regionx+name_gap, regiony-name_gap, name = "segments", vp = vp)

                ############################################
                # Bug: symbols was set to NULL here for some reason
                symbols <- grid::pointsGrob(
                    snp[ind],
                    snp[ind],
                    pch = "*",
                    gp = grid::gpar(
                        cex = 1.25,
                        bg = "blue",
                        col = "blue"
                    ),
                    name = "symbols",
                    vp = vp
                )
                ############################################
                # Figure out exact necessary coefficient for regionx and regiony with name_gap
                SNPnames <-
                    grid::textGrob(
                        SNP.name,
                        just = "left",
                        rot = -45,
                        regionx[ind] - sqrt(2 + 0.5) * name_gap,
                        regiony[ind] + sqrt(2 + 0.5) * name_gap,
                        gp = grid::gpar(cex = 0.6, col = "blue"),
                        name = "SNPnames",
                        vp = vp
                    )
                # Think of better reason to use the +0.5
                # snp[ind], snp[ind], gp = grid::gpar(cex = 0.6, col = "blue"), name = "SNPnames", vp = vp)
                title <-
                    grid::editGrob(title, y = unit(geneMapLabelY + name_gap, "npc"))
            }
            else{
                symbols <- grid::pointsGrob(
                    snp[ind],
                    snp[ind],
                    pch = "*",
                    gp = grid::gpar(
                        cex = 1.25,
                        bg = "blue",
                        col = "blue"
                    ),
                    name = "symbols",
                    vp = vp
                )
                SNPnames <-
                    grid::textGrob(
                        paste(" ", SNP.name),
                        just = "left",
                        rot = -45,
                        regionx[ind],
                        regiony[ind],
                        gp = grid::gpar(cex = 0.6, col = "blue"),
                        name = "SNPnames",
                        vp = vp
                    )
            }
            geneMap <-
                grid::gTree(
                    children = grid::gList(diagonal, segments, title, symbols, SNPnames),
                    name = "geneMap"
                )
        }
    } # if(add.map) end

    else if (!add.map && !is.null(SNP.name) && (any(ind!= 0))) {
        geneMap <- grid::textGrob(
            paste(" ", SNP.name),
            just = "left",
            rot = -45,
            snp[ind],
            snp[ind],
            gp = grid::gpar(cex = 0.6, col = "blue"),
            name = "SNPnames"
        )
        if (flip)
            geneMap <- grid::editGrob(geneMap, vp = vp)
    }
    else
        geneMap <- NULL

    geneMap
}
