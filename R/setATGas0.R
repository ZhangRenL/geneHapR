setATGas0 <- function(gff = gff, hap = hap,
                      geneID = geneID,
                      Chr = Chr, POS = c(start, end)){
    require(IRanges)
    # filter gff by postion
    if(missing(Chr)) {
        warning("Chr is missing, using CHROM in hap")
        Chr <- hap[1,2]
    }
    if(missing(POS)){
        warning("POS is missing, using postion info in hap")
        POS <- hap[2,]
        POS <- suppressWarnings(as.numeric(POS))
        POS <- na.omit(POS)
        POS <- c(min(POS), max(POS))
    } else {
        POS <- c(min(POS), max(POS))
    }
    POSRange <- GenomicRanges::GRanges(
        seqnames = Chr,
        IRanges::IRanges(start = min(POS),
                         end = max(POS)))
    gffOver <- gff[gff %over% POSRange]

    # filter gff by geneID
    if(missing(geneID)){
        stop("geneID is missing")
    } else {
        probe <- stringr::str_detect(gffOver$Parent, geneID)
        probe[is.na(probe)] <- FALSE
        gffOver <- gffOver[probe]
    }

    # get position of ATG
    # filter gff by type: CDS
    gffCDS <- gffOver[gffOver$type == "CDS"]

    strd = GenomicRanges::strand(gffCDS)
    strd = unique(strd)
    if("-" %in% strd & "+" %in% strd)
        stop("Can't decide Gene strand, please check your input")
    if(strd == "+"){
        POSatg <- min(IRanges::start(gffOver))
    } else if(strd == "-"){
        POSatg <- max(IRanges::end(gffOver))
    } else {
        stop("Can't decide Gene strand, please check your input")
    }

    # get postions for hap and gffOver
    newPOS <- suppressWarnings(as.numeric(hap[2,]) - POSatg)
    starts <- IRanges::start(gffOver) - POSatg
    ends <- IRanges::end(gffOver) - POSatg

    # substitute postions in hap and gff
    if(strd == "-") {
        hap[2, -1] <- -newPOS[-1]
        # rename colnames of hap
        probe <- intersect(c("Accession","freq"), colnames(hap))
        newColNam <- hap[2,c(2:(ncol(hap)-length(probe)))]
        colnames(hap)[c(2:(ncol(hap)-length(probe)))] <- newColNam
        # reorder cols of hap
        newOrd <- order(as.numeric(newColNam))
        newOrd <- newColNam[newOrd]
        firColNam <- colnames(hap)[1]
        lastColNam <- colnames(hap)[(ncol(hap) - length(probe)+1):ncol(hap)]
        hap <- hap[,c(firColNam, names(newOrd), lastColNam)]

        IRanges::start(gffOver) <- -ends
        IRanges::end(gffOver) <- -starts
    } else {
        hap[2, -1] <- newPOS[-1]
        IRanges::start(gffOver) <- starts
        IRanges::end(gffOver) <- ends
    }
    GenomicRanges::strand(gffOver) <- "+"

    return(list(gff = gffOver, hap = hap))
}

#' @title gffSetATGas0
#' @description set ATG posision as 0.
#' @importFrom GenomicRanges strand
#' @usage
#' gffSetATGas0(gff = gff, hap = hap,
#'              geneID = geneID,
#'              Chr = Chr, POS = c(start, end))
#' @param gff imported gff
#' @param hap hap result with meta infos
#' @param geneID geneID
#' @param Chr Chrom name
#' @param POS vector defined by start and end position
#' @examples
#' data("quickHap_test")
#' hap <- get_hap(vcf)
#' gffn <- gffSetATGas0(gff = gff, hap = hap,
#'                   geneID = "test1G0387",
#'                   Chr = "scaffold_1",
#'                   POS = c(4300, 7910))
#' @return filtered gff with set ATG position as 0
#' @export
gffSetATGas0 <- function(gff = gff, hap = hap,
                         geneID = geneID,
                         Chr = Chr, POS = c(start, end)){
    tmp <- setATGas0(gff = gff, hap = hap,
                     geneID = geneID,
                     Chr = Chr, POS = POS)
    return(tmp$gff)
}


#' @title hapSetATGas0
#' @description set ATG posision as 0
#' @usage
#' hapSetATGas0(gff = gff, hap = hap,
#'              geneID = geneID,
#'              Chr = Chr, POS = c(start, end))
#' @param gff imported gff
#' @param hap hap result with meta infos
#' @param geneID geneID
#' @param Chr Chrom name
#' @param POS vector defined by start and end position
#' @examples
#' data("quickHap_test")
#' hap <- get_hap(vcf)
#' hapn <- hapSetATGas0(gff = gff, hap = hap,
#'                      geneID = "test1G0387",
#'                      Chr = "scaffold_1",
#'                      POS = c(4300, 7910))
#' @return hap results with set ATG position as 0
#' @export
hapSetATGas0 <- function(gff = gff, hap = hap,
                         geneID = geneID,
                         Chr = Chr, POS = c(start, end)){
    tmp <- setATGas0(gff = gff, hap = hap,
                     geneID = geneID,
                     Chr = Chr, POS = POS)
    hap <- tmp$hap
    return(hap)
}
