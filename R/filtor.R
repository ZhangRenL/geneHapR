#' @name filter_vcf
#' @title Filter VCF
#' @description filter VCF by GFF annotation or by position or both
#' @usage
#' filter_vcf(vcf, gff = gff,
#'            mode = c("POS", "type", "both"),
#'            Chr = Chr, start = start, end = end,
#'            type = c("CDS", "exon", "gene", "genome", "custom"),
#'            cusTyp = c("CDS", "five_prime_UTR", "three_prime_UTR"))
#' @param vcf object of vcfR class, VCF file imported by `import_vcf()`
#' @param gff object of GRanges class, genome annotations imported by
#' `import_gff()`
#' @param mode filter mode, one of "POS", "type", "both"
#' @param Chr chromosome name, needed if mode set to "POS" or "both"
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
#'
#' # filtet vcf
#' data("geneHapR_test")
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
#' @return vcfR
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
        p <- type %in% unique(gff$type)
        m <- paste(unique(gff$type), collapse = "','")
        if (FALSE %in% p)
            stop("type should in c('",m,"')")
        if (type == "custom")
            type <- cusTyp
        if ("genome" %in% type) {
            gff <- gff
        } else {
            gff <- gff[gff$type %in% type]
        }

        POS <- vcfR::getPOS(vcf)
        POS <- as.numeric(POS)
        POSRange <- POS2GRanges(Chr = vcfR::getCHROM(vcf),
                                POS = POS)
        POSRange_rm <- POSRange[!(POSRange %over% gff)]

        POS_rm <- IRanges::start(POSRange_rm)
        probe <- !(POS %in% POS_rm)
        vcf@fix <- vcf@fix[probe,]
        vcf@gt <- vcf@gt[probe,]
    }

    return(vcf)
}


#' @name filter_hap
#' @title Filter hap
#' @description filter hapResult or hapSummary by remove
#'   positions or accessions or haplotypes
#' @usage
#' filter_hap(hap, rm.mode = c("position", "accession", "haplotype"),
#'            position.rm = position.rm,
#'            accession.rm = accession.rm,
#'            haplotype.rm = haplotype.rm)
#' @param hap object of hapSummary or hapResult class
#' @param rm.mode filter mode, one of "position", "accession", "haplotype"
#' @param position.rm numeric vector contains positions need to be removed
#' @param accession.rm character vector contains accessions need to be removed,
#'   only hapResult can be filtered by accessions
#' @param haplotype.rm character vector contains haplotypes need to be removed
#' @param freq.min numeric, hapltypes with accessions number less than freq.min
#' will be removed
#' @examples
#' data("geneHapR_test")
#' hap <- filter_hap(hap,
#'                   rm.mode = c("position", "accession", "haplotype", "freq"),
#'                   position.rm = c(4879, 4950),
#'                   accession.rm = c("C1", "C9"),
#'                   haplotype.rm = c("H009", "H008"),
#'                   freq.min = 5)
#' @return hapSummary or hapResult depend input
#' @export
filter_hap <- function(hap,
                       rm.mode = c("position", "accession", "haplotype", "freq"),
                       position.rm = position.rm,
                       accession.rm = accession.rm,
                       haplotype.rm = haplotype.rm,
                       freq.min = 5){
    if(! (inherits(hap, "hapSummary") | inherits(hap, "hapResult")))
        stop("hap should be class of hapSummary or hapResult")
    cls <- class(hap)
    # attributes
    options <- attr(hap, "options")
    AccAll <- attr(hap, "AccAll")
    AccRemoved <- attr(hap, "AccRemoved")
    hap2acc <- attr(hap, "hap2acc")

    if("accession" %in% rm.mode){
        if(inherits(hap, "hapSummary")){
            warning("Only hapResult class can be filtered by accession")
        } else {
            if(missing(accession.rm))
                warnning("accession.rm is missing")
            rm.probe <- hap$Accession %in% accession.rm
            rm.probe <- sapply(function(x) isTRUE(x))
            hap <- hap[which(!rm.probe),]
            AccRemoved <- c(AccRemoved, accession.rm)
        }
    }

    if("haplotype" %in% rm.mode){
        if(missing(haplotype.rm))
            warnning("haplotype is missing")
        rm.probe <- hap$Hap %in% haplotype.rm
        rm.probe <- sapply(function(x) isTRUE(x))
        hap <- hap[which(! rm.probe),]
        AccRemoved <- c(AccRemoved, hap2acc[names(hap2acc) %in% haplotype.rm])
    }

    if("freq" %in% rm.mode){
        if(missing(freq.min))
            warnning("freq.min is missing")
        rm.probe <- hap$freq < 5
        rm.probe <- sapply(rm.probe, function(x) isTRUE(x))
        hap <- hap[which(! rm.probe),]
        AccRemoved <- c(AccRemoved, hap2acc[! names(hap2acc) %in% hap$Hap])
    }


    if("position" %in% rm.mode){
        if(missing(position.rm))
            stop("position.rm is missing")

        rm.probe <- hap[hap$Hap == "POS",] %>%
            t() %>%
            as.vector()
        rm.probe <- sapply(function(x) isTRUE(x))
        probe <- probe %in% as.character(position.rm)
        hap <- hap[, which(! probe)]
        AccRemoved <- c(AccRemoved, hap2acc[! names(hap2acc) %in% hap$Hap])
    }

    hap2acc <- hap2acc[names(hap2acc) %in% hap$Hap]
    hap <- remove_redundancy_col(hap)
    attr(hap, "options") <- options
    attr(hap, "AccAll") <- AccAll
    attr(hap, "AccRemoved") <- AccRemoved
    attr(hap, "AccRemain") <- AccAll[! AccAll %in% AccRemoved]
    attr(hap, "hap2acc") <- hap2acc
    attr(hap, "freq") <- table(names(hap2acc))
    class(hap) <- cls
    return(hap)
}


POS2GRanges <- function(Chr, POS) {
    POSRange <- GenomicRanges::GRanges(seqnames = Chr,
                                       IRanges::IRanges(start = POS,
                                                        width = rep(1, length(POS))))
    return(POSRange)
}
