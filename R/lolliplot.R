#' #' Lolliplots
#' #' @description Plot variants and somatic mutations
#' #' @param SNP.gr A object of \link[GenomicRanges:GRanges-class]{GRanges},
#' #' \link[GenomicRanges:GRangesList-class]{GRangesList}
#' #' or a list of \link[GenomicRanges:GRanges-class]{GRanges}.
#' #' All the width of GRanges must be 1.
#' #' @param features A object of \link[GenomicRanges:GRanges-class]{GRanges},
#' #' \link[GenomicRanges:GRangesList-class]{GRangesList}
#' #' or a list of \link[GenomicRanges:GRanges-class]{GRanges}.
#' #' The metadata 'featureLayerID' are used for drawing features in different layers.
#' #' See details in vignette.
#' #' @param ranges A object of \link[GenomicRanges:GRanges-class]{GRanges} or
#' #' \link[GenomicRanges:GRangesList-class]{GRangesList}.
#' #' @param type character. Could be circle, pie, pin, pie.stack or flag.
#' #' @param newpage Plot in the new page or not.
#' #' @param ylab Plot ylab or not. If it is a character vector,
#' #' the vector will be used as ylab.
#' #' @param yaxis Plot yaxis or not.
#' #' @param xaxis Plot xaxis or not. If it is a numeric vector with length greater than 1,
#' #' the vector will be used as the points at which tick-marks are to be drawn.
#' #' And the names of the vector will be used to as labels to be placed at the tick
#' #' points if it has names.
#' #' @param ylab.gp,xaxis.gp,yaxis.gp An object of class gpar for ylab, xaxis or yaxis.
#' #' @param legend If it is a list with named color vectors, a legend will be added.
#' #' @param cex cex will control the size of circle.
#' #' @param dashline.col color for the dashed line.
#' #' @param jitter jitter the position of nodes or labels.
#' #' @param rescale logical(1), character(1), numeric vector,
#' #' or a dataframe with rescale from and to.
#' #' Rescalse the x-axis or not.
#' #' if dataframe is used, colnames must be from.start, from.end,
#' #' to.start, to.end. And the from scale must cover the whole plot region.
#' #' The rescale parameter can be set as "exon" or "intron" to
#' #' emphasize "exon" or "intron" region. The "exon" or "intron"
#' #' can be followed with an integer e.g. "exon_80", or "intron_99".
#' #' The integer indicates the total percentage of "exon" or "intron" region.
#' #' Here "exon" indicates all regions in features.
#' #' And "intron" indicates all flank regions of the features.
#' #' @param label_on_feature Labels of the feature directly on them.
#' #' Default FALSE.
#' #' @param ... not used.
#' #' @return NULL
#' #' @details
#' #' In SNP.gr and features, metadata of the GRanges object will be used to control the
#' #' color, fill, border, alpha, shape, height, cex, dashline.col, data source of pie if the type is pie.
#' #' And also the controls for labels by name the metadata start as
#' #' label.parameter.<properties>
#' #' such as label.parameter.rot, label.parameter.gp. The parameter is used for
#' #' \link[grid]{grid.text}. The metadata 'featureLayerID' for features are used
#' #' for drawing features in different layers. The metadata 'SNPsideID' for SNP.gr
#' #' are used for determining the side of lollipops. And the 'SNPsideID' could only
#' #' be 'top' or 'bottom'.
#' #' @return NULL
#' #' @importFrom IRanges IRanges start end NumericList subsetByOverlaps
#' #' elementNROWS countOverlaps disjoin end<- findOverlaps gaps
#' #' reduce start<- tile width
#' #' @importFrom GenomicRanges GRangesList mcols<- GRanges strand<- mcols seqnames
#' #' @importFrom methods is
#' #' @importFrom grid convertHeight convertX convertY gpar grid.legend grid.lines
#' #' grid.newpage grid.polygon grid.rect grid.text grid.xaxis
#' #' grid.yaxis popViewport pushViewport stringHeight stringWidth unit viewport
#' #' @importFrom scales rescale
#' #' @importClassesFrom grImport Picture
#' #' @importFrom grImport readPicture grid.picture
#' #' @importFrom grDevices rgb col2rgb
#' #' @export
#' #' @examples
#' #' SNP <- c(10, 100, 105, 108, 400, 410, 420, 600, 700, 805, 840, 1400, 1402)
#' #' x <- sample.int(100, length(SNP))
#' #' SNP.gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(SNP, width=1, names=paste0("snp", SNP)),
#' #'                                  value1=x, value2=100-x)
#' #' SNP.gr$color <- rep(list(c("red", 'blue')), length(SNP))
#' #' SNP.gr$border <- sample.int(7, length(SNP), replace=TRUE)
#' #' features <- GenomicRanges::GRanges(
#' #'   "chr1", IRanges::IRanges(
#' #'     c(1, 501, 1001),
#' #'     width=c(120, 500, 405),
#' #'     names=paste0("block", 1:3)),
#' #'   color="black",
#' #'   fill=c("#FF8833", "#51C6E6", "#DFA32D"),
#' #'   height=c(0.1, 0.05, 0.08),
#' #'   label.parameter.rot=45)
#' #' lolliplot(SNP.gr, features, type="pie")
#'
#'
#'
#' lolliplot <- function(SNP.gr, features=NULL, ranges=NULL,
#'                       type="circle",
#'                       newpage=TRUE, ylab=TRUE, ylab.gp=gpar(col="black"),
#'                       yaxis=TRUE, yaxis.gp=gpar(col="black"),
#'                       xaxis=TRUE, xaxis.gp=gpar(col="black"),
#'                       legend=NULL, cex=1,
#'                       dashline.col="gray80",
#'                       jitter=c("node", "label"),
#'                       rescale=FALSE,
#'                       label_on_feature=FALSE,
#'                       ...){
#'     stopifnot(inherits(SNP.gr, c("GRanges", "GRangesList", "list")))
#'     stopifnot(inherits(features, c("GRanges", "GRangesList", "list")))
#'     jitter <- match.arg(jitter)
#'     rescale.old <- rescale
#'     xaxis.old <- xaxis
#'     if(any(type!="circle"&jitter=="label")){
#'         jitter[which(type!="circle"&jitter=="label")] <- "node"
#'         warning("if jitter set to label, type must be cirle.")
#'         message("jitter is set to node.")
#'     }
#'     SNP.gr.name <- deparse(substitute(SNP.gr))
#'     if(is(SNP.gr, "GRanges")){
#'         SNP.gr <- list(SNP.gr)
#'         if(length(SNP.gr.name)==length(SNP.gr)) {
#'             names(SNP.gr) <- SNP.gr.name
#'         }
#'     }
#'     len <- length(SNP.gr)
#'     for(i in seq.int(len)){
#'         stopifnot(is(SNP.gr[[i]], "GRanges"))
#'     }
#'
#'     ## prepare the feature
#'     if(inherits(features, c("GRangesList", "list"))){
#'         for(i in seq_along(features)){
#'             stopifnot("features must be a GRanges or GRangesList object"=
#'                           is(features[[i]], "GRanges"))
#'         }
#'         features <- features[seq.int(len)]
#'     }else{
#'         stopifnot("features must be a GRanges or GRangesList object"=
#'                       is(features, "GRanges"))
#'         #features.name <- deparse(substitute(features))
#'         features <- list(features)[seq.int(len)]
#'     }
#'
#'
#'     TYPES <- c("circle", "pie", "pin", "pie.stack", "flag")
#'     if(any(!type %in% TYPES)){
#'         stop("Error in match argument: ",
#'              paste0("'type' should be one of '",
#'                     paste(TYPES, collapse="', '"), "'."))
#'     }
#'     types <- rep(type, length=len)[seq.int(len)]
#'     rm(type)
#'     ############### handle legend ####################
#'     ## set the legend as a list,
#'     ## if all the legend for different tracks is same
#'     ## set draw legend for last track later
#'     legend <- handleLegend(legend, len, SNP.gr)
#'
#'     ################ handle ranges #####################
#'     ## if !missing(ranges) set ranges as feature ranges
#'     ranges <- handleRanges(ranges, SNP.gr, features, len)
#'
#'     ##cut all SNP.gr by the range
#'     SNP.gr <- cutSNP(SNP.gr, ranges, len)
#'
#'
#'     ################## plot ############################
#'     ## total height == 1
#'     height <- 1/sum(lengths(ranges))
#'     args <- as.list(match.call())
#'     if(length(args$height0)==0){
#'         height0 <- 0
#'     }else{
#'         height0 <- args$height0
#'     }
#'     if(newpage) grid.newpage()
#'     for(i in seq.int(len)){
#'         if(length(ranges[[i]])>1){## split into multilayers
#'             args$newpage <- FALSE
#'             for(j in rev(seq_along(ranges[[i]]))){
#'                 args$ranges <- ranges[[i]][j]
#'                 args$SNP.gr <- SNP.gr[i]
#'                 args$features <- features[[i]]
#'                 args$type <- types[i]
#'                 args$legend <- legend[[i]]
#'                 args$height0 <- height0
#'                 height0 <- do.call(what = lolliplot, args = args)
#'             }
#'         }else{## is GRanges with length==1
#'             type <- match.arg(types[i], TYPES)
#'             if(type=="pin"){ ## read the pin shape file
#'                 pinpath <- system.file("extdata", "map-pin-red.xml", package="trackViewer")
#'                 pin <- readPicture(pinpath)
#'             }else{
#'                 pin <- NULL
#'             }
#'             ## Here we don't know the real height of each tracks
#'             vp <- viewport(x=.5, y=height0 + height*0.5, width=1, height=height)
#'             pushViewport(vp)
#'             LINEW <- as.numeric(convertX(unit(1, "line"), "npc"))
#'             LINEH <- as.numeric(convertY(unit(1, "line"), "npc"))
#'             totalH <- as.numeric(unit(1, "npc"))
#'             if(LINEH > totalH/20){
#'                 LINEH <- totalH/20
#'             }
#'
#'             ## GAP the gaps between any elements
#'             GAP <- .2 * LINEH
#'             ratio.yx <- 1/as.numeric(convertX(unit(1, "snpc"), "npc"))
#'
#'
#'             SNPs <- SNP.gr[[i]]
#'             strand(SNPs) <- "*"
#'             SNPs <- sort(SNPs)
#'
#'             feature <- features[[i]]
#'
#'             ## rescale
#'             rescale <- rescale.old
#'             xaxis <- xaxis.old
#'             if(is.logical(rescale)[1]){
#'                 if(rescale[1]){
#'                     range.tile <- tile(ranges[[i]], n = 5)[[1]]
#'                     if(all(width(range.tile)>2)){
#'                         range.tile.cnt <- countOverlaps(range.tile, SNPs)
#'                         feature.start <- feature.end <- feature
#'                         end(feature.start) <- start(feature.start)
#'                         start(feature.end) <- end(feature.end)
#'                         range.tile.cnt2 <- countOverlaps(range.tile, unique(c(feature.start, feature.end)))
#'                         range.tile.cnt <- range.tile.cnt + range.tile.cnt2
#'                         range.width <- width(ranges[[i]])
#'                         range.tile.width <- log2(range.tile.cnt + 1)
#'                         range.tile.width <- range.tile.width/sum(range.tile.width)
#'                         range.tile.width <- range.width * range.tile.width
#'                         range.tile.width <- cumsum(range.tile.width)
#'                         range.tile.width <- start(ranges[[i]]) + c(0, round(range.tile.width)-1)
#'                         rescale <- data.frame(from.start=start(range.tile), from.end=end(range.tile),
#'                                               to.start=range.tile.width[-length(range.tile.width)],
#'                                               to.end=range.tile.width[-1])
#'                         rescale$to.start[-1] <- rescale$to.start[-1] + 1
#'                     }
#'                 }
#'             }else{
#'                 if(is.numeric(rescale)){## rescale is a percentage, convert it into from.start, from.end, to.start, to.end
#'                     ## reset the scale for exons and introns
#'                     feature.rd <- disjoin(c(feature, ranges[[i]]))
#'                     feature.segment.points <- sort(unique(c(start(feature.rd), end(feature.rd))))
#'                     ## filter the points by viewer port.
#'                     feature.segment.points <-
#'                         feature.segment.points[feature.segment.points>=start(ranges[[i]]) &
#'                                                    feature.segment.points<=end(ranges[[i]])]
#'                     rescale <- rescale/sum(rescale, na.rm = TRUE)
#'                     rescale[is.na(rescale)] <- 0.01
#'                     if(length(rescale)==length(feature.segment.points)-1){
#'                         rescale.ir <- IRanges(feature.segment.points[-length(feature.segment.points)]+1,
#'                                               feature.segment.points[-1])
#'                         start(rescale.ir)[1] <- start(rescale.ir)[1]-1
#'                         rescale.ir.width <- sum(width(rescale.ir))
#'                         rescale.ir.new.width <- cumsum(round(rescale.ir.width*rescale, digits = 0))
#'                         rescale <- data.frame(from.start=start(rescale.ir),
#'                                               from.end=end(rescale.ir),
#'                                               to.start=feature.segment.points[1] +
#'                                                   c(0, rescale.ir.new.width[-length(rescale.ir.new.width)]),
#'                                               to.end=feature.segment.points[1] + rescale.ir.new.width)
#'                     }else{
#'                         stop("The length of rescale is not as same as the number of segments (including features and non-features).")
#'                     }
#'                 }else{
#'                     if(is.character(rescale)){
#'                         if(grepl("exon|intron", rescale[1], ignore.case = TRUE)){
#'                             ## reset the scale for exons and introns
#'                             feature.rd <- disjoin(c(reduce(feature), ranges[[i]]), ignore.strand=TRUE)
#'                             feature.rd <- subsetByOverlaps(feature.rd, ranges[[1]], ignore.strand=TRUE)
#'                             feature.segment.exon <- subsetByOverlaps(feature.rd, feature, ignore.strand=TRUE)
#'                             feature.segment.intron <- subsetByOverlaps(feature.rd, feature, invert=TRUE, ignore.strand=TRUE)
#'                             feature.segment.exon$type <- rep("exon", length(feature.segment.exon))
#'                             feature.segment.intron$type <- rep("intron", length(feature.segment.intron))
#'                             feature.rd <- sort(c(feature.segment.exon, feature.segment.intron))
#'                             feature.segment.points <- sort(unique(c(start(feature.segment.exon), end(feature.segment.exon))))
#'                             ## filter the points by viewer port.
#'                             feature.segment.points <-
#'                                 feature.segment.points[feature.segment.points>=start(ranges[[i]]) &
#'                                                            feature.segment.points<=end(ranges[[i]])]
#'                             ratio.to.full.range <- 9
#'                             if(grepl("exon", rescale[1], ignore.case = TRUE)){#set intron width to 1/10
#'                                 if(grepl("^exon_\\d+$", rescale[1], ignore.case = TRUE)){
#'                                     ratio.to.full.range <-
#'                                         as.numeric(sub("^exon_(\\d+)$", "\\1", rescale[1], ignore.case = TRUE))
#'                                     ratio.to.full.range <- ratio.to.full.range/(100-ratio.to.full.range)
#'                                 }
#'                                 width.full.range <- sum(width(feature.rd)[feature.rd$type=="exon"])
#'                                 rescale <- width(feature.rd)
#'                                 rescale[feature.rd$type=="intron"] <-
#'                                     rep(ceiling(width.full.range/ratio.to.full.range/sum(feature.rd$type=="intron")),
#'                                         length(rescale[feature.rd$type=="intron"]))
#'                             }else{#set exon width to 1/10
#'                                 if(grepl("^intron_\\d+$", rescale[1], ignore.case = TRUE)){
#'                                     ratio.to.full.range <-
#'                                         as.numeric(sub("^intron_(\\d+)$", "\\1", rescale[1], ignore.case = TRUE))
#'                                     ratio.to.full.range <- ratio.to.full.range/(100-ratio.to.full.range)
#'                                 }
#'                                 width.full.range <- sum(width(feature.rd)[feature.rd$type=="intron"])
#'                                 rescale <- width(feature.rd)
#'                                 rescale[feature.rd$type=="exon"] <-
#'                                     rep(ceiling(width.full.range/ratio.to.full.range/sum(feature.rd$type=="exon")),
#'                                         length(rescale[feature.rd$type=="exon"]))
#'                             }
#'                             rescale <- rescale/sum(rescale)
#'                             if(length(rescale)==length(feature.segment.points)-1){
#'                                 rescale.ir <- IRanges(feature.segment.points[-length(feature.segment.points)]+1,
#'                                                       feature.segment.points[-1])
#'                                 start(rescale.ir)[1] <- start(rescale.ir)[1]-1
#'                                 rescale.ir.width <- sum(width(rescale.ir))
#'                                 rescale.ir.new.width <- cumsum(round(rescale.ir.width*rescale, digits = 0))
#'                                 rescale <- data.frame(from.start=start(rescale.ir),
#'                                                       from.end=end(rescale.ir),
#'                                                       to.start=feature.segment.points[1] +
#'                                                           c(0, rescale.ir.new.width[-length(rescale.ir.new.width)]),
#'                                                       to.end=feature.segment.points[1] + rescale.ir.new.width)
#'                             }else{
#'                                 stop("Something wrong with the auto-scale setting. Please report the bug to https://github.com/jianhong/trackViewer/issues")
#'                             }
#'                         }
#'                     }
#'                 }
#'             }
#'             if(is.data.frame(rescale)){
#'                 if(all(c("from.start", "from.end", "to.start", "to.end") %in% colnames(rescale))){
#'                     ## check the from coverage the whole region.
#'                     checkflank <- function(x){
#'                         to <- IRanges::IRanges(x$to.start, x$to.end)
#'                         xol <- findOverlaps(to, drop.self=TRUE,
#'                                             drop.redundant=TRUE,minoverlap=2L)
#'                         if(length(xol)>1){
#'                             stop("There is overlaps of the rescale region for 'to' columns.")
#'                         }
#'                         x <- IRanges::IRanges(x$from.start, x$from.end)
#'                         xol <- findOverlaps(x, drop.self=TRUE,
#'                                             drop.redundant=TRUE,minoverlap=2L)
#'                         if(length(xol)>1){
#'                             stop("There is overlaps of the rescale region for 'from' columns.")
#'                         }
#'                         xgap <- gaps(x, start=min(start(x)), end=max(end(x)))
#'                         if(length(xgap)>0){
#'                             stop("There is gaps of the rescale region for 'from' columns.")
#'                         }
#'                     }
#'                     checkflank(rescale)
#'                     rescale.gr <- function(x){
#'                         if(is(x, "GRanges")){
#'                             x.start <- start(x)
#'                             x.end <- end(x)
#'                             y <- c(x.start, x.end)
#'                             x.cut <- cut(y, breaks=c(rescale$from.start[1], rescale$from.end+1),
#'                                          labels=seq.int(nrow(rescale)), right=FALSE)
#'                             y <- mapply(function(a, b){
#'                                 if(!is.na(b)) {
#'                                     rescale(a, to=c(rescale$to.start[b], rescale$to.end[b]),
#'                                             from=c(rescale$from.start[b], rescale$from.end[b]))
#'                                 }else{
#'                                     a
#'                                 }
#'                             }, y, as.numeric(as.character(x.cut)))
#'                             y <- round(y)
#'                             start(x) <- 1
#'                             end(x) <- y[seq_along(x)+length(x)]
#'                             start(x) <- y[seq_along(x)]
#'                             x
#'                         }else{
#'                             x.cut <- cut(x, breaks=c(rescale$from.start[1], rescale$from.end+1),
#'                                          labels=seq.int(nrow(rescale)), right=FALSE)
#'                             y <- mapply(function(a, b){
#'                                 if(!is.na(b)) {
#'                                     rescale(a, to=c(rescale$to.start[b], rescale$to.end[b]),
#'                                             from=c(rescale$from.start[b], rescale$from.end[b]))
#'                                 }else{
#'                                     a
#'                                 }
#'                             }, x, as.numeric(as.character(x.cut)))
#'                             y <- round(y)
#'                             y
#'                         }
#'                     }
#'                     feature <- rescale.gr(feature)
#'                     SNPs <- rescale.gr(SNPs)
#'                     if(is.logical(xaxis)[1]){
#'                         if(xaxis[1]){
#'                             xaxis <- c(rescale$to.start[1], rescale$to.end)
#'                             names(xaxis) <- c(rescale$from.start[1], rescale$from.end)
#'                         }
#'                     }else{
#'                         xaxis.names <- names(xaxis)
#'                         if(length(xaxis.names)!=length(xaxis)){
#'                             xaxis.names <- as.character(xaxis)
#'                         }
#'                         xaxis <- rescale.gr(xaxis)
#'                         names(xaxis) <- xaxis.names
#'                     }
#'                 }
#'             }
#'
#'             ## convert height to npc number
#'             feature$height <- convertHeight2NPCnum(feature$height)
#'             ## multiple transcripts in one gene could be separated by featureLayerID
#'             feature <- setFeatureLayerID(feature, ranges[[i]])
#'             feature.splited <- split(feature, feature$featureLayerID)
#'
#'             ## bottomblank, the transcripts legend height
#'             bottomblank <- plotFeatureLegend(feature, as.numeric(convertY(unit(1, "line"), "npc")),
#'                                              ranges[[i]], xaxis, xaxis.gp, label_on_feature)
#'
#'             ## get the max score and scoreType
#'             if(length(SNPs$score)>0){
#'                 SNPs$score <- sapply(SNPs$score, mean) ## fix the bug of score is NumericList
#'             }
#'             scoreMax0 <- scoreMax <-
#'                 if(length(SNPs$score)>0) ceiling(max(c(SNPs$score, 1), na.rm=TRUE)) else 1
#'             if(type=="pie.stack") scoreMax <- length(unique(SNPs$stack.factor))
#'             if(!type %in% c("pie", "pie.stack")){
#'                 scoreType <-
#'                     if(length(SNPs$score)>0) all(floor(SNPs$score)==SNPs$score) else FALSE
#'                 if(length(yaxis)>1 && is.numeric(yaxis)){
#'                     if(length(names(yaxis))!=length(yaxis)){
#'                         names(yaxis) <- yaxis
#'                     }
#'                     scoreMax0 <- max(yaxis, scoreMax0)
#'                     scoreMax <- max(yaxis, scoreMax)
#'                 }
#'                 if(scoreMax>10) {
#'                     SNPs$score <- 10*SNPs$score/scoreMax
#'                     scoreMax <- 10*scoreMax0/scoreMax
#'                     scoreType <- FALSE
#'                 }else{
#'                     scoreMax <- scoreMax0
#'                 }
#'             }else{
#'                 scoreType <- FALSE
#'             }
#'
#'             ## if the type is caterpillar, there are lollipop in both sides
#'             ## plot the bottom lollipops first. And push a new viewport
#'
#'             IsCaterpillar <- length(SNPs$SNPsideID) > 0
#'             if(IsCaterpillar){
#'                 if(any(is.na(SNPs$SNPsideID)) ||
#'                    !all(SNPs$SNPsideID %in% c('top', 'bottom'))){
#'                     warning("Not all SNPsideID is top or bottom")
#'                     IsCaterpillar <- FALSE
#'                 }
#'             }
#'
#'             if(IsCaterpillar){
#'                 SNPs.top <- SNPs[SNPs$SNPsideID=='top']
#'                 SNPs.bottom <- SNPs[SNPs$SNPsideID=='bottom']
#'             }else{
#'                 SNPs.top <- SNPs
#'                 SNPs.bottom <- GRanges()
#'             }
#'             if(length(SNPs.bottom)<1) IsCaterpillar <- FALSE
#'             ## viewport of plot region
#'             if(!IsCaterpillar){
#'                 bottomblank <- bottomblank
#'             }
#'             pushViewport(viewport(x=LINEW + .5, y=bottomblank/2 + .5,
#'                                   width= 1 - 7*LINEW,
#'                                   height= 1 - bottomblank,
#'                                   xscale=c(start(ranges[[i]]), end(ranges[[i]])),
#'                                   clip="off"))
#'
#'             ## plot xaxis
#'             bottomHeight <- 0
#'             if(IsCaterpillar){
#'                 ## total height == maxscore + extension + gap + labels
#'                 bottomHeight <- getHeight(SNPs=SNPs.bottom,
#'                                           ratio.yx=ratio.yx,
#'                                           LINEW=LINEW,
#'                                           GAP=GAP,
#'                                           cex=cex,
#'                                           type=type,
#'                                           scoreMax=scoreMax,
#'                                           level="data&labels")
#'                 #bottomHeight <- bottomHeight/len
#'                 vp <- viewport(y=bottomHeight, just="bottom",
#'                                xscale=c(start(ranges[[i]]), end(ranges[[i]])))
#'                 pushViewport(vp)
#'                 xaxis.gp$col <- "gray"
#'                 plot.grid.xaxis(xaxis, gp=xaxis.gp)
#'                 popViewport()
#'             }else{
#'                 plot.grid.xaxis(xaxis, gp=xaxis.gp)
#'             }
#'
#'             ## the baseline, the center of the first transcript
#'             baseline <-
#'                 max(c(feature.splited[[1]]$height/2,
#'                       .0001)) + 0.2 * LINEH
#'             baselineN <-
#'                 max(c(feature.splited[[length(feature.splited)]]$height/2,
#'                       .0001)) + 0.2 * LINEH
#'
#'             ##plot features
#'             feature.height <- plotFeatures(feature.splited, LINEH, bottomHeight, label_on_feature)
#'
#'             if(length(SNPs.bottom)>0){
#'                 plotLollipops(SNPs.bottom, feature.height, bottomHeight, baselineN,
#'                               type, ranges[[i]], yaxis, yaxis.gp, scoreMax, scoreMax0, scoreType,
#'                               LINEW, cex, ratio.yx, GAP, pin, dashline.col,
#'                               side="bottom", jitter=jitter)
#'             }
#'             feature.height <- feature.height + 2*GAP
#'             if(length(SNPs.top)>0){
#'                 plotLollipops(SNPs.top, feature.height, bottomHeight, baseline,
#'                               type, ranges[[i]], yaxis, yaxis.gp, scoreMax, scoreMax0, scoreType,
#'                               LINEW, cex, ratio.yx, GAP, pin, dashline.col,
#'                               side="top", jitter=jitter)
#'             }
#'
#'             ## legend
#'             this.height <- getHeight(SNPs.top,
#'                                      ratio.yx, LINEW, GAP, cex, type,
#'                                      scoreMax=scoreMax,
#'                                      level="data&labels")
#'             this.height <- this.height + bottomHeight + feature.height
#'             this.height <- plotLegend(legend[[i]], this.height, LINEH)
#'
#'             popViewport()
#'
#'             this.height <- bottomblank +
#'                 this.height * (1 - bottomblank)
#'
#'             ## ylab
#'             if(length(yaxis)>1 && is.numeric(yaxis)){
#'                 x <- LINEW
#'             }else{
#'                 x <- unit(3, "lines")
#'                 if(yaxis){
#'                     x <- LINEW
#'                 }
#'             }
#'             vp <- viewport(x=.5, y=this.height*0.5,
#'                            width=1, height=this.height)
#'             pushViewport(vp)
#'             if(is.logical(ylab)){
#'                 if(ylab && length(names(SNP.gr))>0){
#'                     grid.text(names(SNP.gr)[i], x = x,
#'                               y = .5, rot = 90, gp=ylab.gp)
#'                 }
#'             }
#'             if(is.character(ylab)){
#'                 if(length(ylab)==1) ylab <- rep(ylab, len)
#'                 grid.text(ylab[i], x = x,
#'                           y = .5, rot = 90, gp=ylab.gp)
#'             }
#'             popViewport()
#'
#'             popViewport()
#'             height0 <- height0 + this.height*height
#'         }
#'     }
#'     return(invisible(height0))
#' }
#'
#' # SNP <- c(10, 100, 105, 108, 400, 410, 420, 600, 700, 805, 840, 1400, 1402)
#' # x <- sample.int(100, length(SNP))
#' # SNP.gr <- GRanges("chr1", IRanges::IRanges(SNP, width=1, names=paste0("snp", SNP)),
#' #   value1=x, value2=100-x)
#' # SNP.gr$color <- rep(list(c("red", 'blue')), length(SNP))
#' # SNP.gr$border <- sample.int(7, length(SNP), replace=TRUE)
#' # features <- GRanges("chr1", IRanges::IRanges(c(1, 501, 1001),
#' #     width=c(120, 500, 405),
#' #     names=paste0("block", 1:3)),
#' #   color="black",
#' #   fill=c("#FF8833", "#51C6E6", "#DFA32D"),
#' #   height=c(0.1, 0.05, 0.08),
#' #   label.parameter.rot=45)
#' # lolliplot(SNP.gr, features, type="pie")
#' #
#'
#'
#' ############### handle legend ####################
#' ## set the legend as a list,
#' ## if all the legend for different tracks is same
#' ## set draw legend for last track later
#' handleLegend <- function(legend, len, dat){
#'     if(length(legend)>0){
#'         if(is.character(legend)){
#'             if(!missing(dat)){
#'                 if(is.list(dat)){
#'                     if(length(legend)<length(dat)){
#'                         legend <- rep(legend, length(dat))[seq_along(dat)]
#'                     }
#'                 }else{
#'                     dat <- list(dat)
#'                     legend <- legend[1]
#'                 }
#'                 para <- c("shape", "color", "border", "alpha")
#'                 preset <- list("circle", "white", "black", 1)
#'                 shapeMap <- c("circle"=21, "square"=22, "diamond"=23,
#'                               "triangle_point_up"=24, "triangle_point_down"=25)
#'                 names(preset) <- para
#'                 legend <- mapply(function(.legend, .dat){
#'                     coln <- colnames(mcols(.dat))
#'                     if(.legend %in% coln){
#'                         labels <- mcols(.dat)[, .legend]
#'                         gp <- lapply(para, function(.ele){
#'                             if(.ele %in% coln){
#'                                 mcols(.dat)[, .ele]
#'                             }else{
#'                                 rep(preset[[.ele]], length(.dat))
#'                             }
#'                         })
#'                         names(gp) <- para
#'                         names(gp)[names(gp)=="color"] <- "fill"
#'                         gp <- as.data.frame(gp, stringsAsFactors=FALSE)
#'                         gp <- cbind(labels=labels, gp)
#'                         gp[, "shape"] <- shapeMap[gp[, "shape"]]
#'                         names(gp)[names(gp)=="shape"] <- "pch"
#'                         gp <- gp[!duplicated(gp[, "labels"]), ]
#'                         gp <- gp[order(gp[, "labels"]), ]
#'                         gp <- as.list(gp)
#'                     }
#'                 }, legend, dat, SIMPLIFY = FALSE)
#'             }
#'         }
#'         if(!is.list(legend)){
#'             tmp <- legend
#'             legend <- vector(mode = "list", length = len)
#'             legend[[len]] <- tmp
#'             rm(tmp)
#'         }else{
#'             if(length(legend)==1){
#'                 tmp <- legend[[1]]
#'                 legend <- vector(mode = "list", length = len)
#'                 legend[[len]] <- tmp
#'                 rm(tmp)
#'             }else{
#'                 if("labels" %in% names(legend)){
#'                     tmp <- legend
#'                     legend <- vector(mode = "list", length = len)
#'                     legend[[len]] <- tmp
#'                     rm(tmp)
#'                 }else{
#'                     if(length(legend)<len){
#'                         length(legend) <- len
#'                     }
#'                 }
#'             }
#'         }
#'     }
#'     return(legend)
#' }
#'
#'
#' ################ handle ranges #####################
#' ## if !missing(ranges) set ranges as feature ranges
#' handleRanges <- function(ranges, SNP.gr, features, len){
#'     if(length(ranges)>0){
#'         stopifnot(inherits(ranges, c("GRanges", "GRangesList", "list")))
#'         if(is(ranges, "GRanges")){
#'             if(length(ranges)==1){
#'                 ranges <- split(rep(ranges, len)[seq.int(len)],
#'                                 seq.int(len))
#'             }else{
#'                 ranges <- split(rep(ranges, len),
#'                                 rep(seq.int(len), each=len))[seq.int(len)]
#'             }
#'         }else{## GRangesList
#'             if(length(ranges)!=len){
#'                 ranges <- rep(ranges, seq.int(len))[seq.int(len)]
#'             }
#'         }
#'         stopifnot(length(ranges)==len)
#'     }else{
#'         if(is(features, "GRanges")){
#'             ranges <- split(range(unname(features), ignore.strand=TRUE)[rep(1, len)],
#'                             seq.int(len))
#'         }else{
#'             if(length(features)!=len){
#'                 stop("if both SNP.gr and features is GRangesList,",
#'                      " the lengthes of them should be identical.")
#'             }
#'             ranges <- GRangesList(lapply(features, function(.ele){
#'                 range(unname(.ele), ignore.strand=TRUE)}))
#'         }
#'     }
#'     return(ranges)
#' }
#'
#'
#' ##cut all SNP.gr by the range
#' cutSNP <- function(SNP.gr, ranges, len){
#'     if(is(ranges, "GRanges")){
#'         for(i in seq.int(len)){
#'             range <- ranges[i]
#'             stopifnot(all(width(SNP.gr[[i]])==1))
#'             SNP.gr[[i]] <- subsetByOverlaps(SNP.gr[[i]], range, ignore.strand=FALSE)
#'         }
#'     }else{
#'         if(inherits(ranges, c("GRangesList", "list"))){
#'             for(i in seq.int(len)){
#'                 range <- ranges[[i]]
#'                 stopifnot(all(width(SNP.gr[[i]])==1))
#'                 SNP.gr[[i]] <- subsetByOverlaps(SNP.gr[[i]], range, ignore.strand=FALSE)
#'             }
#'         }
#'     }
#'     return(SNP.gr)
#' }
#'
#'
#' convertHeight2NPCnum <- function(.ele){
#'     if(is(.ele, "unit")){
#'         return(convertHeight(.ele, unitTo="npc", valueOnly=TRUE))
#'     }else{
#'         if(is.list(.ele)){
#'             .ele <- sapply(.ele, function(.e){
#'                 if(is(.e, "unit")){
#'                     .e <- convertHeight(.e, unitTo="npc", valueOnly=TRUE)
#'                 }
#'                 .e[1]
#'             })
#'             return(unlist(.ele))
#'         }else{
#'             if(is.numeric(.ele)){
#'                 return(.ele)
#'             }else{
#'                 if(is.integer(.ele)){
#'                     return(.ele)
#'                 }else{
#'                     return(.ele)
#'                 }
#'             }
#'         }
#'     }
#' }
#'
#'
#' ## multiple transcripts in one gene could be separated by featureLayerID
#' setFeatureLayerID <- function(feature, range){
#'     feature <- feature[end(feature)>=start(range) &
#'                            start(feature)<=end(range)]
#'     if(length(feature$featureLayerID)!=length(feature)){
#'         feature$featureLayerID <- rep("1", length(feature))
#'         feature$featureLayerID <- as.character(feature$featureLayerID)
#'         start(feature)[start(feature)<start(range)] <- start(range)
#'         end(feature)[end(feature)>end(range)] <- end(range)
#'     }
#'     return(feature)
#' }
#'
#'
#' ## bottomblank, the transcripts legend height
#' plotFeatureLegend <- function(feature, LINEH, range, xaxis, xaxis.gp, label_on_feature=FALSE){
#'     if(label_on_feature) return(0)
#'     if(length(xaxis)>1 || as.logical(xaxis[1])){
#'         xaxisSpace <- 2
#'         if(is.numeric(xaxis.gp$cex)) xaxisSpace <- 2*xaxis.gp$cex
#'     }else{
#'         xaxisSpace <- 0
#'     }
#'     if(length(names(feature))>0){ ## features legend
#'         feature.s <- feature[!duplicated(names(feature))]
#'         cex <- if(length(unlist(feature.s$cex))==length(feature.s))
#'             unlist(feature.s$cex) else 1
#'         ncol <- getColNum(names(feature.s), cex=cex)
#'         featureLegendSpace <- max(ceiling(length(names(feature.s)) / ncol) * cex + 1 )
#'         pushViewport(viewport(x=.5, y=featureLegendSpace*LINEH/2,
#'                               width=1,
#'                               height=featureLegendSpace*LINEH,
#'                               xscale=c(start(range), end(range))))
#'         color <- if(length(unlist(feature.s$color))==length(feature.s))
#'             unlist(feature.s$color) else "black"
#'         fill <- if(length(unlist(feature.s$fill))==length(feature.s))
#'             unlist(feature.s$fill) else "black"
#'         pch <- if(length(unlist(feature.s$pch))==length(feature.s))
#'             unlist(feature.s$pch) else 22
#'         grid.legend(label=names(feature.s), ncol=ncol,
#'                     byrow=TRUE, vgap=unit(.2, "lines"),
#'                     hgap=unit(.5, "lines"),
#'                     pch=pch,
#'                     gp=gpar(col=color, fill=fill, cex=cex))
#'         popViewport()
#'     }else{
#'         featureLegendSpace <- 0
#'     }
#'     bottomblank <- (xaxisSpace + featureLegendSpace) * LINEH
#'     return(bottomblank)
#' }
#'
#'
#' getColNum <- function(labels, spaces="WW", cex){
#'     ncol <- floor(as.numeric(convertX(unit(1, "npc"), "line")) /
#'                       maxStringWidth(labels, spaces=spaces, cex) /
#'                       as.numeric(convertX(stringWidth("W"), "line")))
#'     nrow <- ceiling(length(labels) / ncol)
#'     ncol <- ceiling(length(labels) / nrow)
#'     ncol
#' }
#'
#'
#' maxStringWidth <- function(labels, spaces="WW", cex){
#'     max(as.numeric(convertX(stringWidth(paste0(labels, spaces)), "line"))*cex)
#' }
#'
#'
#' plot.grid.xaxis <- function(xaxis, gp=gpar(col="black")){
#'     ## axis, should be in the bottom of transcripts
#'     if(length(xaxis)==1 && as.logical(xaxis)) {
#'         grid.xaxis(gp=gp)
#'     }
#'     if(length(xaxis)>1 && is.numeric(xaxis)){
#'         xaxisLabel <- names(xaxis)
#'         if(length(xaxisLabel)!=length(xaxis)) xaxisLabel <- TRUE
#'         grid.xaxis(at=xaxis, label=xaxisLabel, gp=gp)
#'     }
#' }
#'
#'
#' plotFeatures <- function(feature.splited, LINEH, bottomHeight,
#'                          label_on_feature=FALSE){
#'     feature.height <- 0
#'     for(n in seq_along(feature.splited)){
#'         this.feature.height <-
#'             max(c(feature.splited[[n]]$height/2,
#'                   .0001)) + 0.2 * LINEH
#'         feature.height <- feature.height + this.feature.height
#'         ##baseline
#'         grid.lines(x=c(0, 1), y=c(bottomHeight+feature.height,
#'                                   bottomHeight+feature.height))
#'         for(m in seq_along(feature.splited[[n]])){
#'             this.dat <- feature.splited[[n]][m]
#'             color <- if(is.list(this.dat$color)) this.dat$color[[1]] else
#'                 this.dat$color
#'             if(length(color)==0) color <- "black"
#'             fill <- if(is.list(this.dat$fill)) this.dat$fill[[1]] else
#'                 this.dat$fill
#'             if(length(fill)==0) fill <- "white"
#'             this.cex <- if(length(this.dat$cex)>0) this.dat$cex[[1]][1] else 1
#'             lwd <- if(length(this.dat$lwd)>0) this.dat$lwd[[1]][1] else 1
#'             this.feature.height.m <-
#'                 if(length(this.dat$height)>0)
#'                     this.dat$height[[1]][1] else
#'                         2*this.feature.height
#'             grid.rect(x=start(this.dat)-.1, y=bottomHeight+feature.height,
#'                       width=width(this.dat)-.8,
#'                       height=this.feature.height.m,
#'                       just="left", gp=gpar(col=color, fill=fill, lwd=lwd),
#'                       default.units = "native")
#'             if(label_on_feature & !is.na(names(this.dat)[1])){
#'                 grid.text(x=(start(this.dat)+end(this.dat))/2,
#'                           y=bottomHeight+feature.height,
#'                           just = "centre",
#'                           label = names(this.dat)[1],
#'                           gp= gpar(list(cex=this.cex *
#'                                             this.feature.height.m/
#'                                             this.feature.height,
#'                                         color=color)),
#'                           default.units = "native")
#'             }
#'         }
#'         feature.height <- feature.height + this.feature.height
#'     }
#'     feature.height
#' }
#'
#'
#' plotLollipops <- function(SNPs, feature.height, bottomHeight, baseline,
#'                           type, ranges, yaxis, yaxis.gp, scoreMax, scoreMax0, scoreType,
#'                           LINEW, cex, ratio.yx, GAP, pin, dashline.col,
#'                           side=c("top", "bottom"), jitter=c("node", "label"),
#'                           main=TRUE){
#'     side <- match.arg(side)
#'     jitter <- match.arg(jitter)
#'     if(side=="top"){
#'         pushViewport(viewport(y=bottomHeight,
#'                               height=1,
#'                               just="bottom",
#'                               xscale=c(start(ranges),
#'                                        end(ranges)),
#'                               clip="off"))
#'     }else{
#'         pushViewport(viewport(y=bottomHeight+feature.height,
#'                               height=1,
#'                               just="top",
#'                               xscale=c(start(ranges),
#'                                        end(ranges)),
#'                               yscale=c(1, 0),
#'                               clip="off"))
#'     }
#'     if(type=="pie.stack" && length(SNPs$stack.factor)>0){
#'         stopifnot(is.vector(SNPs$stack.factor, mode="character"))
#'         if(length(SNPs$stack.factor.order)>0 ||
#'            length(SNPs$stack.factor.first)>0){
#'             warning("stack.factor.order and stack.factor.first are used by this function!",
#'                     "The values in these column will be removed.")
#'         }
#'         ## condense the SNPs
#'         stack.factors <- unique(as.character(SNPs$stack.factor))
#'         stack.factors <- sort(stack.factors)
#'         stack.factors.order <- seq_along(stack.factors)
#'         names(stack.factors.order) <- stack.factors
#'         SNPs <- SNPs[order(as.character(seqnames(SNPs)), start(SNPs),
#'                            as.character(SNPs$stack.factor))]
#'         SNPs$stack.factor.order <- stack.factors.order[SNPs$stack.factor]
#'         SNPs$stack.factor.first <- !duplicated(SNPs)
#'         SNPs.condense <- SNPs
#'         SNPs.condense$oid <- seq_along(SNPs)
#'         SNPs.condense$factor <- paste(as.character(seqnames(SNPs)), start(SNPs), end(SNPs))
#'         SNPs.condense <- split(SNPs.condense, SNPs.condense$factor)
#'         SNPs.condense <- lapply(SNPs.condense, function(.ele){
#'             .oid <- .ele$oid
#'             .gr <- .ele[1]
#'             mcols(.gr) <- NULL
#'             .gr$oid <- NumericList(.oid)
#'             .gr
#'         })
#'         SNPs.condense <- unlist(GRangesList(SNPs.condense), use.names = FALSE)
#'         SNPs.condense <- sort(SNPs.condense)
#'         lab.pos.condense <- jitterLables(start(SNPs.condense),
#'                                          xscale=c(start(ranges), end(ranges)),
#'                                          lineW=LINEW*cex)
#'         lab.pos.condense <- reAdjustLabels(lab.pos.condense,
#'                                            lineW=LINEW*cex)
#'         condense.ids <- SNPs.condense$oid
#'         lab.pos <- rep(lab.pos.condense, elementNROWS(condense.ids))
#'         lab.pos <- lab.pos[order(unlist(condense.ids))]
#'     }else{
#'         lab.pos <- jitterLables(start(SNPs),
#'                                 xscale=c(start(ranges), end(ranges)),
#'                                 lineW=LINEW*cex)
#'         lab.pos <- reAdjustLabels(lab.pos,
#'                                   lineW=LINEW*cex)
#'     }
#'
#'     if(length(SNPs)>0){
#'         yaxisat <- NULL
#'         yaxisLabel <- TRUE
#'         if(length(yaxis)>1 && is.numeric(yaxis)){
#'             yaxisat <- yaxis
#'             if(length(names(yaxis))==length(yaxis)) yaxisLabel <- names(yaxis)
#'             yaxis <- TRUE
#'         }
#'         if(yaxis && scoreMax>1 && !type %in% c("pie", "pie.stack")){
#'             if(side=="top"){
#'                 grid.yaxis(at=yaxisat,
#'                            label=yaxisLabel,
#'                            main = main,
#'                            gp=yaxis.gp,
#'                            vp=viewport(x=.5+ifelse(main, -1, 1) *LINEW,
#'                                        y=feature.height+5.25*GAP*cex+
#'                                            scoreMax*LINEW*ratio.yx/2*cex,
#'                                        width=1,
#'                                        height=scoreMax*LINEW*ratio.yx*cex,
#'                                        yscale=c(0, scoreMax0+.5)))
#'             }else{
#'                 grid.yaxis(at=yaxisat,
#'                            label=yaxisLabel,
#'                            main = main,
#'                            gp=yaxis.gp,
#'                            vp=viewport(x=.5+ifelse(main, -1, 1) *LINEW,
#'                                        y=1-(feature.height+5.25*GAP*cex+
#'                                                 scoreMax*LINEW*ratio.yx/2*cex),
#'                                        width=1,
#'                                        height=scoreMax*LINEW*ratio.yx*cex,
#'                                        yscale=c(scoreMax0+.5, 0)))
#'             }
#'         }
#'         if(length(SNPs$alpha)==length(SNPs)){
#'             SNPs$alpha[is.na(SNPs$alpha)] <- 0
#'             if(all(is.numeric(SNPs$alpha))){
#'                 if(any(SNPs$alpha>1)){## convert to 0-1
#'                     SNPs$alpha <- SNPs$alpha/max(SNPs$alpha)
#'                 }
#'             }else{ ## not correct format.
#'                 SNPs$alpha <- as.numeric(factor(as.character(SNPs$alpha)))
#'                 SNPs$alpha <- (SNPs$alpha+max(SNPs$alpha))/max(SNPs$alpha)/2
#'             }
#'         }else{
#'             SNPs$alpha <- NULL
#'         }
#'         if(type=="circle"){
#'             if(length(SNPs$shape)==length(SNPs)){
#'                 ## shape could only be "circle", "square", "diamond", "triangle_point_up", "triangle_point_down"
#'                 if(!all(SNPs$shape %in% c("circle", "square", "diamond", "triangle_point_up", "triangle_point_down"))){
#'                     message('shape must be "circle", "square", "diamond", "triangle_point_up", or "triangle_point_down"')
#'                     SNPs$shape <- as.numeric(factor(SNPs$shape))
#'                     SNPs$shape <- rep(c("circle", "square", "diamond", "triangle_point_up", "triangle_point_down"),
#'                                       max(SNPs$shape))[SNPs$shape]
#'                 }
#'             }else{
#'                 SNPs$shape <- NULL
#'             }
#'         }
#'         for(m in seq_along(SNPs)){
#'             this.dat <- SNPs[m]
#'             color <- if(is.list(this.dat$color)) this.dat$color[[1]] else this.dat$color
#'             border <-
#'                 if(is.list(this.dat$border)) this.dat$border[[1]] else this.dat$border
#'             fill <- if(is.list(this.dat$fill)) this.dat$fill[[1]] else this.dat$fill
#'             alpha <- if(length(this.dat$alpha)>0) this.dat$alpha[[1]] else 1
#'             lwd <- if(is.list(this.dat$lwd)) this.dat$lwd[[1]] else this.dat$lwd
#'             id <- if(is.character(this.dat$label)) this.dat$label else NA
#'             id.col <- if(length(this.dat$label.col)>0) this.dat$label.col else "black"
#'             shape <- if(length(this.dat$shape)>0) this.dat$shape[[1]] else "circle"
#'             rot <- if(length(this.dat$label.rot)>0) this.dat$label.rot else 15
#'             this.cex <- if(length(this.dat$cex)>0) this.dat$cex[[1]][1] else 1
#'             this.dashline.col <-
#'                 if(length(this.dat$dashline.col)>0) this.dat$dashline.col[[1]][1] else dashline.col
#'             if(length(names(this.dat))<1) this.dashline.col <- NA
#'             this.dat.mcols <- mcols(this.dat)
#'             this.dat.mcols <- cleanDataMcols(this.dat.mcols, type)
#'
#'             grid.lollipop(x1=convertX(unit(start(this.dat), "native"), "npc",
#'                                       valueOnly=TRUE),
#'                           y1=baseline,
#'                           x2=convertX(unit(ifelse(jitter=="node",
#'                                                   lab.pos[m],
#'                                                   start(this.dat)),
#'                                            "native"), "npc", valueOnly=TRUE),
#'                           y2=feature.height,
#'                           y3=4*GAP*cex, y4=2.5*GAP*cex,
#'                           radius=LINEW*cex/2,
#'                           col=color,
#'                           border=border,
#'                           percent=this.dat.mcols,
#'                           edges=100,
#'                           type=type,
#'                           ratio.yx=ratio.yx,
#'                           pin=pin,
#'                           scoreMax=scoreMax * LINEW * cex,
#'                           scoreType=scoreType,
#'                           id=id, id.col=id.col,
#'                           cex=this.cex, lwd=lwd, dashline.col=this.dashline.col,
#'                           side=side, rot=rot, alpha=alpha, shape=shape)
#'
#'         }
#'         this.height <- getHeight(SNPs,
#'                                  ratio.yx, LINEW, GAP, cex, type,
#'                                  scoreMax=scoreMax,
#'                                  level="data")
#'         labels.rot <- 90
#'         if(length(names(SNPs))>0){
#'             if(type=="pie.stack"){
#'                 ## unique lab.pos and SNPs
#'                 idx <- !duplicated(names(SNPs))
#'                 lab.pos <- lab.pos[idx]
#'                 SNPs <- SNPs[idx]
#'             }
#'             labels.x <- lab.pos
#'             labels.text <- names(SNPs)
#'             labels.just <- ifelse(side=="top", "left", "right")
#'             labels.hjust <- NULL
#'             labels.vjust <- NULL
#'             labels.check.overlap <- FALSE
#'             labels.default.units <- "native"
#'             labels.gp <- gpar(cex=cex)
#'
#'             ## change the parameter by use definations.
#'             for(label.parameter in c("x", "y", "just", "hjust", "vjust",
#'                                      "rot", "check.overlap", "default.units",
#'                                      "gp")){
#'                 label.para <- paste0("label.parameter.", label.parameter)
#'                 if(label.para %in% colnames(mcols(SNPs))){
#'                     assign(paste0("labels.", label.parameter),
#'                            mcols(SNPs)[, label.para])
#'                 }
#'             }
#'             if(!"cex" %in% names(labels.gp)){
#'                 labels.gp <- c(labels.gp, cex=cex)
#'             }
#'             mergeList <- function(.ele){
#'                 .n <- unique(unlist(lapply(.ele, names)))
#'                 .out <- list()
#'                 if(length(.n)>0){
#'                     for(.name in .n){
#'                         .out[[.name]] <- sapply(.ele, function(.e){
#'                             if(.name %in% names(.e)){
#'                                 .e[[.name]][1]
#'                             }else{
#'                                 NA
#'                             }
#'                         })
#'                     }
#'                 }else{
#'                     .n <- unique(names(.ele))
#'                     for(.name in .n){
#'                         .out[[.name]] <- unlist(.ele[names(.ele) %in% .name])
#'                     }
#'                 }
#'                 .out
#'             }
#'             labels.gp <- mergeList(labels.gp)
#'             labels.gp[duplicated(names(labels.gp))] <- NULL
#'             labels.gp <- do.call(gpar, labels.gp)
#'             if(jitter=="label"){
#'                 ## add guide lines
#'                 rased.height <- 4*GAP*cex
#'                 guide.height <- 2.5*GAP*cex
#'                 for(i in seq_along(SNPs)){
#'                     this.dashline.col <-
#'                         if(length(SNPs[i]$dashline.col)>0)
#'                             SNPs[i]$dashline.col[[1]][1] else
#'                                 dashline.col
#'                     if(length(names(SNPs[i]))<1) this.dashline.col <- NA
#'                     grid.lines(x=c(start(SNPs[i]), labels.x[i]),
#'                                y=c(this.height+feature.height-cex*LINEW,
#'                                    this.height+feature.height+rased.height),
#'                                default.units = labels.default.units,
#'                                gp=gpar(col=this.dashline.col, lty=3))
#'                     grid.lines(x=c(labels.x[i], labels.x[i]),
#'                                y=c(this.height+rased.height+feature.height,
#'                                    this.height+rased.height+
#'                                        guide.height+feature.height),
#'                                default.units = labels.default.units,
#'                                gp=gpar(col=this.dashline.col, lty=3))
#'                 }
#'                 ## add this height
#'                 this.height <- this.height + rased.height + guide.height
#'             }
#'             grid.text(x=labels.x, y=this.height + feature.height,
#'                       label = labels.text,
#'                       just = labels.just,
#'                       hjust = labels.hjust,
#'                       vjust = labels.vjust,
#'                       rot=labels.rot,
#'                       check.overlap = labels.check.overlap,
#'                       default.units = labels.default.units,
#'                       gp=labels.gp)
#'         }
#'     }
#'     popViewport()
#' }
#'
#'
#' jitterLables <- function(coor, xscale, lineW, weight=1.2){
#'     if(weight==1.2) {
#'         stopifnot("Please sort your inputs by start position"=
#'                       order(coor)==1:length(coor))
#'     }
#'     if(weight<0.5) return(coor)
#'     stopifnot(length(xscale)==2)
#'     pos <- convertX(unit(coor, "native"), "npc", valueOnly=TRUE)
#'     pos.diff <- diff(c(0, pos, 1))
#'     idx <- which(pos.diff < weight*lineW)
#'     if(length(idx)<1){
#'         return(coor)
#'     }
#'     if(all(idx %in% c(1, length(pos)+1))){
#'         return(coor)
#'     }
#'     idx.diff <- diff(c(-1, idx))
#'     idx.grp <- rle(idx.diff)
#'     idx.grp$values[idx.grp$values==1] <- length(pos) + 1:sum(idx.grp$values==1)
#'     idx.grp <- inverse.rle(idx.grp)
#'     idx.grp.w <- which(idx.grp>length(pos))-1
#'     idx.grp.w <- idx.grp.w[idx.grp.w>0]
#'     idx.grp[idx.grp.w] <- idx.grp[idx.grp.w+1]
#'     idx.grp <- split(idx, idx.grp)
#'     flag <- as.numeric(names(idx.grp))>length(pos)
#'     idx.grp.mul <- lapply(idx.grp[flag], function(.ele){
#'         c(.ele[1]-1, .ele)
#'     })
#'     idx.grp.sin <- lapply(idx.grp[!flag], function(.ele){
#'         lapply(as.list(.ele), function(.ele){c(.ele-1, .ele)})
#'     })
#'     idx.grp.sin <- unlist(idx.grp.sin, recursive = FALSE)
#'     idx.grp <- c(idx.grp.mul, idx.grp.sin)
#'
#'     adj.pos <- lapply(idx.grp, function(.ele){
#'         .ele <- .ele[.ele>0 & .ele<=length(pos)]
#'         this.pos <- pos[.ele]
#'         names(this.pos) <- .ele
#'         if(length(this.pos)%%2==1){
#'             center <- ceiling(length(this.pos)/2)
#'         }else{
#'             center <- length(this.pos)/2 + .5
#'         }
#'         if(length(this.pos)>5){ ## too much, how to jitter?
#'             this.pos <- this.pos +
#'                 ((1:length(this.pos))-center) * (weight-.1) *
#'                 lineW/ceiling(log(length(this.pos), 5))
#'         }else{
#'             this.pos <- this.pos +
#'                 ((1:length(this.pos))-center) * (weight-.1) * lineW
#'         }
#'         this.pos
#'     })
#'     names(adj.pos) <- NULL
#'     adj.pos <- unlist(adj.pos)
#'     coor[as.numeric(names(adj.pos))] <- adj.pos*diff(xscale)+xscale[1]
#'
#'     Recall(coor, xscale=xscale, lineW=lineW, weight=weight-0.2)
#' }
#'
#'
#' reAdjustLabels <- function(coor, lineW){
#'     # resort
#'     coor <- sort(coor)
#'     bins <- ceiling(1/lineW)
#'     pos <- convertX(unit(coor, "native"), "npc", valueOnly=TRUE)
#'     pos.bin <- cut(pos, c(-Inf, (0:bins)*lineW, Inf), labels=0:(bins+1), right=FALSE)
#'
#'     ## split the coors by into clusters
#'     ## give the clusters with more idx more spaces if there are spaces between clusters
#'     tbl <- table(pos.bin)
#'     if(all(tbl<2)) return(coor)
#'     tbl.len <- length(tbl)
#'     if(tbl.len<3) return(coor)
#'     loops <- 1000
#'     loop <- 1
#'     while(any(tbl==0) && any(tbl>1) && loop < loops){
#'         tbl.bk <- tbl
#'         for(i in order(tbl.bk, decreasing=TRUE)){
#'             if(tbl[i]>1 && tbl.bk[i]==tbl[i]){
#'                 if(i==1){
#'                     if(tbl[2]<tbl[1]){
#'                         half <- sum(tbl[1:2])/2
#'                         tbl[2] <- ceiling(half)
#'                         tbl[1] <- floor(half)
#'                     }
#'                 }else{
#'                     if(i==tbl.len){
#'                         if(tbl[tbl.len]>tbl[tbl.len-1]){
#'                             half <- sum(tbl[(tbl.len-1):tbl.len])/2
#'                             tbl[tbl.len-1] <- ceiling(half)
#'                             tbl[tbl.len] <- floor(half)
#'                         }
#'                     }else{
#'                         if(tbl[i-1]<tbl[i+1]){
#'                             ## i-1 and i should be balanced
#'                             half <- sum(tbl[(i-1):i])/2
#'                             tbl[i-1] <- floor(half)
#'                             tbl[i] <- ceiling(half)
#'                         }else{
#'                             half <- sum(tbl[i:(i+1)])/2
#'                             tbl[i] <- floor(half)
#'                             tbl[i+1] <- ceiling(half)
#'                         }
#'                     }
#'                 }
#'             }
#'         }
#'         loop <- loop + 1
#'     }
#'     coef <- unlist(lapply(tbl, function(.ele){
#'         if(.ele==0) return(0)
#'         .ele <- seq(from=0, to=1, length.out=.ele+1)
#'         (.ele[-length(.ele)] + .ele[-1])/2
#'     }))
#'     coef <- coef[coef!=0]
#'     coor <- (rep(as.numeric(names(tbl)), tbl) - 1 + coef) * lineW
#'     coor <- convertX(unit(coor, "npc"), "native", valueOnly=TRUE)
#'     coor
#' }
#'
#'
#' cleanDataMcols <- function(this.dat.mcols, type){
#'     this.dat.mcols <-
#'         this.dat.mcols[,
#'                        !colnames(this.dat.mcols) %in%
#'                            c("color", "fill", "lwd", "id",
#'                              "cex", "dashline.col",
#'                              "id.col", "stack.factor", "SNPsideID",
#'                              "shape", "alpha"),
#'                        drop=FALSE]
#'     if(type!="pie.stack"){
#'         this.dat.mcols <-
#'             this.dat.mcols[, !colnames(this.dat.mcols) %in%
#'                                c("stack.factor.order",
#'                                  "stack.factor.first"),
#'                            drop=FALSE]
#'     }
#'     this.dat.mcols <-
#'         this.dat.mcols[, !grepl("^label.parameter",
#'                                 colnames(this.dat.mcols)),
#'                        drop=FALSE]
#'     return(this.dat.mcols)
#' }
#'
#'
#' grid.lollipop <- function (x1=.5, y1=.5,
#'                            x2=.5, y2=.75,
#'                            y3=.04, y4=.02,
#'                            radius=.8,
#'                            col=NULL,
#'                            border=NULL,
#'                            percent=NULL,
#'                            edges=100,
#'                            type=c("circle", "pie", "pin", "pie.stack", "flag"),
#'                            ratio.yx=1,
#'                            pin=NULL,
#'                            scoreMax,
#'                            scoreType,
#'                            id=NA, id.col="black",
#'                            cex=1, lwd=1,
#'                            dashline.col="gray80",
#'                            side=c("top", "bottom"),
#'                            rot=15,
#'                            alpha=NULL,
#'                            shape=shape){
#'     side <- match.arg(side)
#'     stopifnot(is.numeric(c(x1, x2, y1, y2, y3, y4, radius, edges)))
#'     type <- match.arg(type)
#'     side <- side!="top"
#'     if(!type %in% c("pie", "pie.stack")){
#'         this.score <- if(length(percent$score)>0) max(percent$score, 1) else 1
#'         if(type=="circle"){
#'             y0 <- c(y1, y2, y2+y3, y2+y3+y4+(this.score-1)*2*radius*ratio.yx+(1-cex)*radius*ratio.yx)
#'             if(scoreType) y0[4] <- y2+y3+y4
#'             if(side) y0 <- 1 - y0
#'             grid.lines(x=c(x1, x1, x2, x2), y=y0,
#'                        gp=gpar(col=border, lwd=lwd))
#'             y0 <- c(y2+y3+y4+this.score*2*radius*ratio.yx,
#'                     y2+y3+y4+scoreMax*ratio.yx)
#'             if(scoreType) y0[1] <- y2+y3+y4+this.score*2*radius*ratio.yx*cex
#'             if(side) y0 <- 1 - y0
#'             grid.lines(x=c(x2, x2),
#'                        y=y0,
#'                        gp=gpar(col=dashline.col, lty=3, lwd=lwd))
#'         }else{
#'             y0 <- c(y1, y2, y2+y3, y2+y3+y4+(this.score-.5)*2*radius*ratio.yx)
#'             if(side) y0 <- 1 - y0
#'             grid.lines(x=c(x1, x1, x2, x2), y=y0,
#'                        gp=gpar(col=border, lwd=lwd))
#'         }
#'
#'     }else{
#'         if(type=="pie.stack"){
#'             if(percent$stack.factor.first){
#'                 y0 <- c(y1, y2, y2+y3, y2+y3+y4)
#'                 if(side) y0 <- 1 - y0
#'                 grid.lines(x=c(x1, x1, x2, x2), y=y0,
#'                            gp=gpar(col=border, lwd=lwd))
#'                 y0 <- c(y2+y3+y4, y2+y3+y4+scoreMax*ratio.yx)
#'                 if(side) y0 <- 1 - y0
#'                 grid.lines(x=c(x2, x2),
#'                            y=y0,
#'                            gp=gpar(col=dashline.col, lty=3, lwd=lwd))
#'             }
#'         }else{
#'             y0 <- c(y1, y2, y2+y3, y2+y3+y4)
#'             if(side) y0 <- 1 - y0
#'             grid.lines(x=c(x1, x1, x2, x2), y=y0,
#'                        gp=gpar(col=border, lwd=lwd))
#'         }
#'     }
#'     if(length(pin)>0){
#'         if(length(border)>0) pin@paths[[2]]@rgb <- rgb2hex(col2rgb(border[1]))
#'         if(length(col)>0) pin@paths[[1]]@rgb <- rgb2hex(col2rgb(col[1]))
#'         if(length(col)>1) pin@paths[[3]]@rgb <- rgb2hex(col2rgb(col[2]))
#'     }
#'     switch(type,
#'            circle={
#'                if(length(border)==0) border <- "black"
#'                if(length(col)==0) col <- "white"
#'                if(scoreType){
#'                    for(i in 1:this.score){
#'                        y0 <- y2+y3+y4+2*radius*ratio.yx*(i-.5)*cex
#'                        if(side) y0 <- 1 - y0
#'                        switch(shape, #"circle", "square", "diamond", "triangle_point_up", "star", or "triangle_point_down"
#'                               circle=grid.circle1(x=x2, y=y0,
#'                                                   r=radius*ratio.yx*cex,
#'                                                   gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
#'                               square=grid.square(x=x2, y=y0,
#'                                                  r=radius*ratio.yx*cex,
#'                                                  gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
#'                               diamond=grid.diamond(x=x2, y=y0,
#'                                                    r=radius*ratio.yx*cex,
#'                                                    gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
#'                               triangle_point_up=grid.triangle_point_up(x=x2, y=y0,
#'                                                                        r=radius*ratio.yx*cex,
#'                                                                        gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
#'                               triangle_point_down=grid.triangle_point_down(x=x2, y=y0,
#'                                                                            r=radius*ratio.yx*cex,
#'                                                                            gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
#'                               star=grid.star(x=x2, y=y0,
#'                                              r=radius*ratio.yx*cex,
#'                                              gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
#'                               grid.circle1(x=x2, y=y0,
#'                                            r=radius*ratio.yx*cex,
#'                                            gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)))
#'
#'                    }
#'                }else{
#'                    y0 <- y2+y3+(this.score-.5)*2*radius*ratio.yx+y4
#'                    if(side) y0 <- 1 - y0
#'                    switch(shape,
#'                           circle=grid.circle1(x=x2, y=y0,
#'                                               r=radius*ratio.yx*cex,
#'                                               gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
#'                           square=grid.square(x=x2, y=y0,
#'                                              r=radius*ratio.yx*cex,
#'                                              gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
#'                           diamond=grid.diamond(x=x2, y=y0,
#'                                                r=radius*ratio.yx*cex,
#'                                                gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
#'                           triangle_point_up=grid.triangle_point_up(x=x2, y=y0,
#'                                                                    r=radius*ratio.yx*cex,
#'                                                                    gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
#'                           triangle_point_down=grid.triangle_point_down(x=x2, y=y0,
#'                                                                        r=radius*ratio.yx*cex,
#'                                                                        gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
#'                           star=grid.star(x=x2, y=y0,
#'                                          r=radius*ratio.yx*cex,
#'                                          gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)),
#'                           grid.circle1(x=x2, y=y0,
#'                                        r=radius*ratio.yx*cex,
#'                                        gp=gpar(col=border, fill=col, lwd=lwd, alpha=alpha)))
#'                    if(!is.na(id)){
#'                        y0 <- y2+y3+(this.score-.5)*2*radius*ratio.yx+y4
#'                        if(side) y0 <- 1 - y0
#'                        grid.text(label=id, x=x2,
#'                                  y=y0,
#'                                  just="centre", gp=gpar(col=id.col, cex=.75*cex))
#'                    }
#'                }
#'            },
#'            pie={
#'                y0 <- y2+y3+y4+radius*ratio.yx
#'                if(side) y0 <- 1 - y0
#'                grid.pie(x=x2, y=y0,
#'                         radius = radius*cex,
#'                         col = col,
#'                         border = border,
#'                         percent=percent,
#'                         edges=edges,
#'                         lwd=lwd, alpha=alpha)
#'            },
#'            pie.stack={
#'                y0 <- y2+y3+y4+(2*percent$stack.factor.order-1)*radius*ratio.yx
#'                if(side) y0 <- 1 - y0
#'                grid.pie(x=x2,
#'                         y=y0,
#'                         radius = radius*cex,
#'                         col = col,
#'                         border = border,
#'                         percent=percent[, !colnames(percent) %in%
#'                                             c("stack.factor.order",
#'                                               "stack.factor.first")],
#'                         edges=edges,
#'                         lwd=lwd, alpha=alpha)
#'            },
#'            pin={
#'                y0 <- y2+y3+(this.score-.5)*2*radius*ratio.yx+y4/2
#'                if(side) y0 <- 1 - y0
#'                grid.picture(picture=pin, x=x2,
#'                             y=y0,
#'                             width=2*radius*ratio.yx*cex,
#'                             height=3*radius*ratio.yx*cex+y4)
#'                if(!is.na(id)){
#'                    y0 <- y2+y3+(this.score-.25)*2*radius*ratio.yx+2*y4/3
#'                    grid.text(label=id, x=x2,
#'                              y=y0,
#'                              just="centre", gp=gpar(col=id.col, cex=.5*cex))
#'                }
#'            },
#'            flag={
#'                if(is.na(id)){
#'                    id <- " "
#'                }
#'                LINEH <- as.numeric(convertY(unit(1, "line"), "npc"))*cex
#'                y0 <- y2+y3+(this.score-.5)*2*radius*ratio.yx+y4/2
#'                if(side) y0 <- 1 - y0
#'                LINEW <- as.numeric(convertX(stringWidth(paste0("o", id, "u")), "npc"))*cex
#'                LINEW <- LINEW * sign(cos(pi*rot/180))
#'                LINEH0 <- LINEW*ratio.yx*tan(pi*rot/180)
#'                grid.polygon(x=c(x2, x2+LINEW, x2+LINEW, x2),
#'                             y=c(y0, y0+LINEH0, y0+LINEH0+LINEH*1.25, y0+LINEH*1.25),
#'                             gp=gpar(fill=col, col=border, alpha=alpha))
#'                grid.text(label=id, x=x2+LINEW*.5,
#'                          y=y0 + LINEH*.625+LINEH0*.5,
#'                          hjust=.5, vjust=.5,
#'                          gp=gpar(col=id.col, cex=cex),
#'                          rot=rot)
#'            },
#'            grid.pie(x=x2, y=y2+y3+y4+radius*ratio.yx,
#'                     radius = radius*cex,
#'                     col = col,
#'                     border = border,
#'                     percent=percent,
#'                     edges=edges,
#'                     lwd=lwd, alpha=alpha))
#' }
#'
#'
#' grid.pie <- function (x=.5, y=.5,
#'                       radius=.8,
#'                       col=NULL,
#'                       border=NULL,
#'                       percent=NULL,
#'                       edges=100,
#'                       lwd=1,
#'                       alpha=1) {
#'     if(length(percent)>0) percent <- unlist(percent[, sapply(percent, is.numeric)])
#'     if(length(percent)<1){
#'         percent <- 1
#'     }
#'     percent <- c(0, cumsum(percent)/sum(percent))
#'     if(any(is.na(percent))){
#'         warning("There are events with NA number after calculating the percentage.",
#'                 "Please make sure all the events must contain at least one values greater than 0")
#'         percent[is.na(percent)] <- 0
#'     }
#'     dx <- diff(percent)
#'     nx <- length(dx)
#'     if (is.null(col))
#'         col <- c("white", "lightblue", "mistyrose", "lightcyan",
#'                  "lavender", "cornsilk")
#'     if (!is.null(col))
#'         col <- rep_len(col, nx)
#'     if (!is.null(border))
#'         border <- rep_len(border, nx)
#'     twopi <- 2 * pi
#'     ratio.yx <- 1/as.numeric(convertX(unit(1, "snpc"), "npc"))
#'     t2xy <- function(t) {
#'         t2p <- twopi * t + pi/2
#'         list(x = radius * cos(t2p), y = radius * sin(t2p) * ratio.yx)
#'     }
#'     for (i in 1L:nx) {
#'         n <- max(2, floor(edges * dx[i]))
#'         P <- t2xy(seq.int(percent[i], percent[i + 1], length.out = n))
#'         grid.polygon(unit(c(P$x, 0)+x,"npc"), unit(c(P$y, 0)+y, "npc"), gp=gpar(col = border[i], fill = col[i], lwd=lwd, alpha=alpha))
#'     }
#'     invisible(NULL)
#' }
#'
#'
#' getHeight <- function(SNPs, ratio.yx, LINEW, GAP, cex, type, scoreMax,
#'                       level=c("data", "data&labels")){
#'     level=match.arg(level)
#'     stack.factors <- unique(as.character(SNPs$stack.factor))
#'     stack.factors <- sort(stack.factors)
#'     if(level=="data"){
#'         switch(type,
#'                circle={
#'                    labels.y <- LINEW + # add gaps for labels
#'                        6.5*GAP*cex +
#'                        scoreMax * LINEW * ratio.yx*cex
#'                },
#'                pin={
#'                    if(length(SNPs$score)>0) {
#'                        this.scores <- ceiling(SNPs$score)
#'                    }else {
#'                        this.scores <- .5
#'                    }
#'                    this.scores[is.na(this.scores)] <- .5
#'                    labels.y <- LINEW +
#'                        6.5*GAP*cex +
#'                        (this.scores-0.5) * LINEW * ratio.yx*cex
#'                },
#'                pie={
#'                    labels.y <- LINEW*max(ratio.yx, 1.2) +
#'                        6.5*GAP*cex + 0.5 * LINEW * ratio.yx * cex
#'                },
#'                pie.stack={
#'                    labels.y <- LINEW +
#'                        6.5*GAP*cex +
#'                        (scoreMax-0.5) * LINEW * ratio.yx*cex
#'                },
#'                flag={
#'                    labels.y <- LINEW +
#'                        6.5*GAP*cex +
#'                        scoreMax * LINEW * ratio.yx*cex
#'                })
#'         labels.y
#'     }else{
#'         if(length(SNPs$label.parameter.rot)>0) {
#'             labels.rot <- SNPs$label.parameter.rot
#'         }else{
#'             labels.rot <- 90
#'         }
#'         labels.cex <- 1
#'         if(length(SNPs$label.parameter.gp)>0){
#'             if(length(SNPs$label.parameter.gp$cex)>0)
#'                 labels.cex <- SNPs$label.parameter.gp$cex[[1]][1]
#'         }
#'         labels.length.rate <- labels.cex * max(cospi((labels.rot-90)/180), 0) * ratio.yx
#'         stringH <- as.numeric(convertY(stringHeight("W"), "npc"))
#'
#'         switch(type,
#'                circle={
#'                    if(length(names(SNPs))>0){
#'                        maxStrHeight <-
#'                            max(as.numeric(
#'                                convertX(stringWidth(names(SNPs)), "npc")
#'                            ))+stringH
#'                    }else{
#'                        maxStrHeight <- 0
#'                    }
#'                    maxStrHeight <- maxStrHeight * labels.length.rate
#'                    ypos <- LINEW + 6.5*GAP*cex +
#'                        scoreMax * LINEW * ratio.yx*cex + maxStrHeight
#'                },
#'                pin={
#'                    if(length(names(SNPs))>0){
#'                        thisStrHeight <- max(as.numeric(
#'                            convertX(stringWidth(names(SNPs)), "npc")) ) +
#'                            stringH
#'                    }else{
#'                        thisStrHeight <- 0
#'                    }
#'                    thisStrHeight <- thisStrHeight * labels.length.rate
#'                    if(length(SNPs$score)>0){
#'                        ypos <-
#'                            max(LINEW +
#'                                    6.5*GAP*cex +
#'                                    (SNPs$score-0.5) * LINEW * ratio.yx*cex +
#'                                    thisStrHeight)
#'                    }else{
#'                        ypos <- max(LINEW*max(ratio.yx, 1.2) +
#'                                        6.5*GAP*cex + thisStrHeight)
#'                    }
#'                },
#'                pie={
#'                    if(length(names(SNPs))>0){
#'                        maxStrHeight <-
#'                            max(as.numeric(
#'                                convertX(stringWidth(names(SNPs)), "npc")
#'                            ))+stringH
#'                    }else{
#'                        maxStrHeight <- 0
#'                    }
#'                    maxStrHeight <- maxStrHeight * labels.length.rate
#'                    ypos <- LINEW +
#'                        6.5*GAP*cex + maxStrHeight
#'                },
#'                pie.stack={
#'                    if(length(names(SNPs))>0){
#'                        maxStrHeight <-
#'                            max(as.numeric(
#'                                convertX(stringWidth(names(SNPs)), "npc")
#'                            ))+stringH
#'                    }else{
#'                        maxStrHeight <- 0
#'                    }
#'                    maxStrHeight <- maxStrHeight * labels.length.rate
#'                    ypos <- LINEW +
#'                        6.5*GAP*cex + maxStrHeight +
#'                        (scoreMax-0.5) * LINEW * ratio.yx*cex
#'                },
#'                flag={
#'                    if(length(names(SNPs))>0){
#'                        maxStrHeight <-
#'                            max(as.numeric(
#'                                convertX(stringWidth(names(SNPs)), "npc")
#'                            ))+stringH
#'                    }else{
#'                        maxStrHeight <- 0
#'                    }
#'                    maxStrHeight <- maxStrHeight * labels.length.rate
#'                    ypos <- LINEW + 6.5*GAP*cex +
#'                        scoreMax * LINEW * ratio.yx*cex + maxStrHeight
#'                }
#'         )
#'         ypos
#'     }
#' }
#'
#'
#' plotLegend <- function(legend, this.height, LINEH){
#'     ypos <- this.height
#'     pch <- 21
#'     if(length(legend)>0){
#'         if(is.list(legend)){
#'             thisLabels <- legend[["labels"]]
#'             if("pch" %in% names(legend)) pch <- legend[["pch"]]
#'             gp <- legend[!names(legend) %in% c("labels", "pch")]
#'             if(is.null(gp$cex)) gp$cex <- 1
#'             class(gp) <- "gpar"
#'         }else{
#'             thisLabels <- names(legend)
#'             gp <- gpar(fill=legend, cex=1)
#'         }
#'         if(length(thisLabels)>0){
#'             ncol <- getColNum(thisLabels, cex=gp$cex)
#'             topblank <- ceiling(length(thisLabels) / ncol) * gp$cex[1]
#'             pushViewport(viewport(x=.5,
#'                                   y=ypos+(topblank+.2*gp$cex[1])*LINEH/2,
#'                                   width=1,
#'                                   height=topblank*LINEH,
#'                                   just="bottom"))
#'             this.height <- ypos + (topblank+.2*gp$cex[1])*LINEH
#'             grid.legend(label=thisLabels, ncol=ncol,
#'                         byrow=TRUE, vgap=unit(.1*gp$cex[1], "lines"),
#'                         hgap=unit(.5*gp$cex[1], "lines"),
#'                         pch=pch,
#'                         gp=gp)
#'             popViewport()
#'         }
#'     }
#'     this.height + LINEH
#' }
#'
#'
#' #"circle", "square", "diamond", "triangle_point_up", "star", or "triangle point_down"
#' grid.circle1 <- function(x = 0.5, y = 0.5, r = 0.5,
#'                          default.units = "npc", name = NULL,
#'                          gp = gpar(), draw = TRUE, vp = NULL){
#'     fill <- gp$fill
#'     col <- gp$col
#'     lwd <- if(length(gp$lwd)>0) gp$lwd else 1
#'     alpha <- gp$alpha
#'     if(is.null(fill)) fill <- "white"
#'     twopi <- 2 * pi
#'     ratio.yx <- 1/as.numeric(convertX(unit(1, "snpc"), "npc"))
#'     t2xy <- function(t) {
#'         t2p <- twopi * t + pi/2
#'         list(x = r * cos(t2p)/ratio.yx, y = r * sin(t2p))
#'     }
#'     P <- t2xy(seq.int(0, 1, length.out = 100))
#'     invisible(grid.polygon(unit(P$x+x,"npc"), unit(P$y+y, "npc"),
#'                            gp=gpar(col = col, fill = fill, lwd=lwd, alpha=alpha)))
#' }
#'
#'
#' grid.diamond <- function(x = 0.5, y = 0.5, r = 0.5,
#'                          default.units = "npc", name = NULL,
#'                          gp = gpar(), draw = TRUE, vp = NULL){
#'     fill <- gp$fill
#'     col <- gp$col
#'     lwd <- if(length(gp$lwd)>0) gp$lwd else 1
#'     alpha <- gp$alpha
#'     if(is.null(fill)) fill <- "white"
#'     ratio.yx <- 1/as.numeric(convertX(unit(1, "snpc"), "npc"))
#'     P <-
#'         list(x = c(0, r/ratio.yx, 0, -r/ratio.yx),
#'              y = c(-r, 0, r, 0))
#'     invisible(grid.polygon(unit(P$x+x,"npc"), unit(P$y+y, "npc"),
#'                            gp=gpar(col = col, fill = fill, lwd=lwd, alpha=alpha)))
#' }
#'
#'
#' grid.triangle_point_up <- function(x = 0.5, y = 0.5, r = 0.5,
#'                                    default.units = "npc", name = NULL,
#'                                    gp = gpar(), draw = TRUE, vp = NULL){
#'     fill <- gp$fill
#'     col <- gp$col
#'     lwd <- if(length(gp$lwd)>0) gp$lwd else 1
#'     alpha <- gp$alpha
#'     if(is.null(fill)) fill <- "white"
#'     ratio.yx <- 1/as.numeric(convertX(unit(1, "snpc"), "npc"))
#'     P <-
#'         list(x = c(-r/ratio.yx, r/ratio.yx, 0, -r/ratio.yx),
#'              y = c(-r, -r, r, -r))
#'     invisible(grid.polygon(unit(P$x+x,"npc"), unit(P$y+y, "npc"),
#'                            gp=gpar(col = col, fill = fill, lwd=lwd, alpha=alpha)))
#' }
#'
#'
#' grid.triangle_point_down <- function(x = 0.5, y = 0.5, r = 0.5,
#'                                      default.units = "npc", name = NULL,
#'                                      gp = gpar(), draw = TRUE, vp = NULL){
#'     fill <- gp$fill
#'     col <- gp$col
#'     lwd <- if(length(gp$lwd)>0) gp$lwd else 1
#'     alpha <- gp$alpha
#'     if(is.null(fill)) fill <- "white"
#'     ratio.yx <- 1/as.numeric(convertX(unit(1, "snpc"), "npc"))
#'     P <-
#'         list(x = c(-r/ratio.yx, r/ratio.yx, 0, -r/ratio.yx),
#'              y = c(r, r, -r, r))
#'     invisible(grid.polygon(unit(P$x+x,"npc"), unit(P$y+y, "npc"),
#'                            gp=gpar(col = col, fill = fill, lwd=lwd, alpha=alpha)))
#' }
#'
#'
#' grid.star <- function(x = 0.5, y = 0.5, r = 0.5,
#'                       default.units = "npc", name = NULL,
#'                       gp = gpar(), draw = TRUE, vp = NULL){
#'     fill <- gp$fill
#'     col <- gp$col
#'     lwd <- if(length(gp$lwd)>0) gp$lwd else 1
#'     alpha <- gp$alpha
#'     if(is.null(fill)) fill <- "white"
#'     ratio.yx <- 1/as.numeric(convertX(unit(1, "snpc"), "npc"))
#'     i <- 1:11
#'     angle <- 180
#'     alpha <- 2*pi / 10
#'     r <- r * (i %% 2 + 1)/2
#'     omega <- alpha * i + angle * pi /180
#'     invisible(grid.polygon(unit(r*sin(omega)/ratio.yx+x,"npc"),
#'                            unit(r*cos(omega)+y, "npc"),
#'                            gp=gpar(col = col, fill = fill, lwd=lwd, alpha=alpha)))
#' }
#'
#'
#' rgb2hex <- function(x){
#'     hex <- function(a) format(as.hexmode(a), width=2, upper.case=TRUE)
#'     if(length(x==3))
#'         paste0("#",hex(x[1]),hex(x[2]),hex(x[3]))
#'     else
#'         paste0("#",hex(x[1]),hex(x[2]),hex(x[3]),hex(x[4]))
#' }
#'
#'
#' grid.square <- function(x = 0.5, y = 0.5, r = 0.5,
#'                         default.units = "npc", name = NULL,
#'                         gp = gpar(), draw = TRUE, vp = NULL){
#'     fill <- gp$fill
#'     col <- gp$col
#'     lwd <- if(length(gp$lwd)>0) gp$lwd else 1
#'     alpha <- gp$alpha
#'     if(is.null(fill)) fill <- "white"
#'     ratio.yx <- 1/as.numeric(convertX(unit(1, "snpc"), "npc"))
#'     invisible(grid.rect(unit(x,"npc"), unit(y, "npc"),
#'                         width = unit(r*2/ratio.yx, "npc"),
#'                         height = unit(r*2, "npc"),
#'                         gp=gpar(col = col, fill = fill, lwd=lwd, alpha=alpha)))
#' }
