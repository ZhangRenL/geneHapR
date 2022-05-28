POS2GRanges <- function(Chr, POS){
    POSRange <- GenomicRanges::GRanges(
        seqnames = Chr,
        IRanges::IRanges(
            start = POS,
            width = rep(1, length(POS))))
    return(POSRange)
}


#' @description 根据变异位点所在区域的类型进行筛选
#' @title filter_hap_by_gff
#' @usage
#' filter_hap(gff = gff,
#'            hap = hap,
#'            hapPrefix_new = "C_H",
#'            type = c("CDS", "exon", "gene","genome", "custom"),
#'            custom_type = c("CDS", "exon", "gene"))
#' @param gff gff文件
#' @param hap hap或hapResult
#' @param hapPrefix_new new prefix for filtered haps
#' @param type 筛选模式，CDS/exon/gene/genome/custom之一，
#' 如果是custom则必须设定custom_type
#' @param custom_type 自定义筛选模式
#' @importFrom IRanges start
#' @importFrom IRanges `%over%`
#' @examples
#' # filtet hap
#' data("quickHap_test")
#' hap <- get_hap(vcf)
#' hap <- filter_hap(gff = gff, hap = hap, type = "CDS")
#'
#' @export
filter_hap <- function(gff = gff,
                       hap = hap,
                       hapPrefix_new = "C_H",
                       type = c("CDS", "exon", "gene", "genome", "custom"),
                       custom_type = c("CDS", "exon", "gene")){
    if(!type %in% c("CDS", "exon", "gene","genome", "custom"))
        stop('type must be one of  c("CDS", "exon", "gene","genome", "custom")')
    if(type == "custom") type <- custom_type
    # filter gff file by type type
    if(type == "genome") {
        gff <- gff
    } else {
        gff <- gff[gff$type %in% type]
    }
    meta <- hap[1:4,]
    POS <- suppressWarnings(as.numeric(meta[2,]))
    POS <- na.omit(POS)
    Chr <- meta[1,2]
    hapRange <- POS2GRanges(Chr = Chr, POS = POS)
    hapRange <- hapRange[!(hapRange %over% gff)]
    POS_rm <- IRanges::start(hapRange)
    probe <- !(colnames(hap) %in% POS_rm)
    hap <- hap[,probe]
    meta <- hap[1:4,]
    hap <- hap[5:nrow(hap), -1]
    hap <- assign_hapID(hap, hapPrefix = hapPrefix_new)
    hap <- rbind(meta, hap)
    return(hap)
}



#' @description filter vcf by gff file or by position or both
#' @title filter_vcf_by_gff
#' @usage
#' filter_vcf(vcf, gff = gff,
#'            mode = c("POS", "type", "both"),
#'            Chr = Chr, start = start, end = end,
#'            type = c("CDS", "exon", "gene", "genome", "custom"),
#'            cusTyp = c("CDS", "five_prime_UTR", "three_prime_UTR"))
#' @param vcf imported vcf
#' @param gff imported gff
#' @param mode filter mode, one of POS/type/both
#' @param Chr CHROM name, needed if mode set to 'POS' or 'both'
#' @param start Start position, needed if mode set to 'POS' or 'both'
#' @param end End position, needed if mode set to 'POS' or 'both'
#' @param type filter type, needed if mode set to 'type' or 'both',
#' one of CDS/exon/gene/genome/custom,
#' if type set to custom, the custom_type is needed.
#' @param cusTyp vector, custom filter type, needed if type set to custom
#' @importFrom IRanges start
#' @importFrom IRanges `%over%`
#' @examples
#' # filtet hap
#' data("quickHap_test")
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
filter_vcf <- function(vcf, gff = gff,
                       mode = c("POS", "type", "both"),
                       Chr = Chr, start = start, end = end,
                       type = c("CDS", "exon", "gene", "genome", "custom"),
                       cusTyp = c("CDS", "five_prime_UTR", "three_prime_UTR")){

    if(mode == "POS" | mode == "both"){
        if(missing(Chr)) stop("Chr is missing!")
        if(missing(start)) stop("start is missing!")
        if(missing(end)) stop("end is missing!")
        POS <- vcfR::getPOS(vcf)
        POS <- as.numeric(POS)
        Chrs <- vcfR::getCHROM(vcf)
        probe <- POS >= min(start, end) & POS <= max(start, end)
        probe <- probe & Chrs == Chr
        vcf@fix <- vcf@fix[probe,]
        vcf@gt <- vcf@gt[probe,]
        }

    if(mode == "type" | mode == "both"){
        if(missing(gff)) stop("gff is missing!")
        if(type == "custom") type <- cusTyp
        if(length(type) != 1 )
            stop('Type must be one of c("CDS","exon","gene","genome","custom")')
        if(type == "genome") {
            gff <- gff
        } else {
            gff <- gff[gff$type %in% type]
        }
        if(missing(Chr)) Chr <- vcfR::getCHROM(vcf)[1]
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
