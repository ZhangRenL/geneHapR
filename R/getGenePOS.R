#' @title Get Gene Position
#' @name getGenePOS
#' @export
#' @usage
#' getGenePOS(gff= gff,
#'            geneID = geneID,
#'            type = type,
#'            gffTermContaingeneID = "Parent")
#' @param gff imported gff
#' @param geneID target geneID
#' @param type vector consist with one or more types in gff
#' @param gffTermContaingeneID which term contains the geneID in your gff,
#' defalt is Parent
#' @examples
#' data("geneHapR_test")
#' genePOS <- getGenePOS(gff = gff,
#'            geneID = "test1G0387",
#'            type = "CDS",
#'            gffTermContaingeneID = "Parent")
#'
#' @return named vectors contains start, end and strand
#' @importFrom IRanges start
#' @importFrom IRanges end
#' @importFrom GenomicRanges strand
getGenePOS <- function(gff = gff,
                       geneID = geneID,
                       type = type,
                       gffTermContaingeneID = "Parent") {
    # get target gene range
    if (missing(type)) {
        gffGene <- getGeneRanges(
            gff = gff,
            geneID = geneID,
            gffTermContaingeneID = "Parent"
        )
    } else {
        gffGene <- getGeneRanges(
            gff = gff,
            geneID = geneID,
            type = type,
            gffTermContaingeneID = "Parent"
        )

    }

    # Get start, end and strand
    strd <- unique(GenomicRanges::strand(range(gffGene)))
    strd <- as.vector(strd)
    strts <- IRanges::start(gffGene)
    ends <- IRanges::end(gffGene)
    strt <- min(strts, ends)
    end <- max(strts, ends)
    genePOSInfo <- c(strt, end, strd)
    names(genePOSInfo) <- c("start", "end", "strand")
    return(genePOSInfo)
}


#' @title Get Gene Ranges
#' @name getGeneRanges
#' @param gff imported gff
#' @param geneID target geneID
#' @param type vector consist with one or more types in gff
#' @param gffTermContaingeneID which term contains the geneID in your gff,
#' defalt is Parent
#' @usage
#' getGeneRanges(gff= gff,
#'              geneID = geneID,
#'              type = type,
#'              gffTermContaingeneID = "Parent")
#' @examples
#' data("geneHapR_test")
#' geneRanges <- getGeneRanges(gff = gff,
#'                             geneID = "test1G0387",
#'                             type = "CDS",
#'                             gffTermContaingeneID = "Parent")
#' @importFrom stringr str_detect
#' @importFrom GenomicRanges as.data.frame
#' @return GRanges
#' @export
getGeneRanges <- function(gff = gff,
                          geneID = geneID,
                          type = type,
                          gffTermContaingeneID = "Parent") {
    # set probes for target gene
    probeCol <- GenomicRanges::as.data.frame(gff)
    probeCol <- probeCol[, gffTermContaingeneID]
    probeCol <- unlist(probeCol)
    probe <- stringr::str_detect(probeCol, geneID)
    probe[is.na(probe)] <- FALSE

    # get ranges of target gene
    if (!TRUE %in% probe)
        stop("Can't found geneID in current gff")
    gffGene <- gff[probe]

    # filter range types
    if (!missing(type)) {
        probe <- gffGene$type %in% type
        probe[is.na(probe)] <- FALSE
        gffGene <- gffGene[probe]
    }
    return(gffGene)
}
