# checked
#' @name filter_vcf_by_gff
#' @title filter vcf by gff
#' @description filter vcf by gff annotation or by position or both
#' @usage
#' filter_vcf(vcf, gff = gff,
#'            mode = c("POS", "type", "both"),
#'            Chr = Chr, start = start, end = end,
#'            type = c("CDS", "exon", "gene", "genome", "custom"),
#'            cusTyp = c("CDS", "five_prime_UTR", "three_prime_UTR"))
#' @param vcf object of `vcfR` class, vcf file imported by `import_vcf()`
#' @param gff object of `GRanges` class, genome annotations imported by
#' `import_gff()`
#' @param mode filter mode, one of "POS", "type", "both"
#' @param Chr CHROM name, needed if mode set to "POS" or "both"
#' @param start start position, needed if mode set to "POS" or "both"
#' @param end end position, needed if mode set to "POS" or "both"
#' @param type filter type, needed if mode set to "type" or "both",
#' one of "CDS", "exon", "gene", "genome", "custom",
#' if `type` was set to "custom", then `custom_type` is needed.
#' @param cusTyp character vector, custom filter type,
#' needed if `type` set to "custom"
#' @importFrom IRanges start
#' @importFrom IRanges `%over%`
#' @examples
#' # filtet hap
#' data("geneHap_test")
#' vcf_f1 <- filter_vcf(vcf, mode = "POS",
#'                     Chr = "scaffold_1",
#'                     start = 4300, end = 5890)
#'
#' vcf_f2 <- filter_vcf(vcf, mode = "type",
#'                     gff = gff,
#'                     type = "CDS")
#' vcf_f3 <- filter_vcf(vcf, mode = "both",
#'                     Chr = "scaffold_1",
#'                     start = 4300, end = 5890,
#'                     gff = gff,
#'                     type = "CDS")
#'
#' @export
filter_vcf <- function(vcf,
                       gff = gff,
                       mode = c("POS", "type", "both"),
                       Chr = Chr,
                       start = start,
                       end = end,
                       type = c("CDS", "exon", "gene", "genome", "custom"),
                       cusTyp = c("CDS", "five_prime_UTR", "three_prime_UTR")) {
    if (mode == "POS" | mode == "both") {
        if (missing(Chr))
            stop("Chr is missing!")
        if (missing(start))
            stop("start is missing!")
        if (missing(end))
            stop("end is missing!")
        POS <- vcfR::getPOS(vcf)
        POS <- as.numeric(POS)
        Chrs <- vcfR::getCHROM(vcf)
        probe <- POS >= min(start, end) & POS <= max(start, end)
        probe <- probe & Chrs == Chr
        vcf@fix <- vcf@fix[probe,]
        vcf@gt <- vcf@gt[probe,]
    }

    if (mode == "type" | mode == "both") {
        if (missing(gff))
            stop("gff is missing!")
        if (type == "custom")
            type <- cusTyp
        if (length(type) != 1)
            stop('Type must be one of c("CDS", "exon", "gene", "genome", "custom")')
        if (type == "genome") {
            gff <- gff
        } else {
            gff <- gff[gff$type %in% type]
        }
        if (missing(Chr))
            Chr <- vcfR::getCHROM(vcf)[1]
        POS <- vcfR::getPOS(vcf)
        POS <- as.numeric(POS)
        POSRange <- POS2GRanges(Chr = Chr, POS = POS)
        POSRange_rm <- POSRange[!(POSRange %over% gff)]

        POS_rm <- IRanges::start(POSRange_rm)
        probe <- !(POS %in% POS_rm)
        vcf@fix <- vcf@fix[probe,]
        vcf@gt <- vcf@gt[probe,]
    }

    return(vcf)
}




POS2GRanges <- function(Chr, POS) {
    POSRange <- GenomicRanges::GRanges(seqnames = Chr,
                                       IRanges::IRanges(start = POS,
                                                        width = rep(1, length(POS))))
    return(POSRange)
}
