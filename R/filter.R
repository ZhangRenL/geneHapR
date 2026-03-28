#' @name filter_vcf
#' @title Filter VCF
#' @description filter VCF by GFF annotation or by position or both
#' @usage
#' filter_vcf(vcf, gff = gff,
#'            mode = c("POS", "type", "both"),
#'            Chr = Chr, start = start, end = end,
#'            type = c("CDS", "exon", "gene", "genome", "custom"),
#'            cusTyp = c("CDS", "five_prime_UTR", "three_prime_UTR"),
#'            geneID = geneID)
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
#' @param geneID gene ID
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
                       cusTyp = c("CDS", "five_prime_UTR", "three_prime_UTR"),
                       geneID = geneID) {
    if (mode == "POS" | mode == "both") {
        POS <- vcfR::getPOS(vcf) %>% as.numeric()
        Chrs <- vcfR::getCHROM(vcf)
        probe <- rep(TRUE, length(POS))
        if(!missing(Chr))
            probe <- probe & Chrs == Chr
        if(!missing(start))
            probe <- probe & (POS >= min(start, end))
        if(!missing(end))
            probe <- probe & (POS <= max(start, end))

        nr <- sum(probe)
        if(nr > 1){
            vcf@fix <- vcf@fix[probe,]
            vcf@gt <- vcf@gt[probe,]
        } else if(nr == 1){
            nf <- colnames(vcf@fix)
            ng <- colnames(vcf@gt)
            vcf@fix <- matrix(vcf@fix[probe,], byrow = TRUE, nrow = nr)
            colnames(vcf@fix) <- nf
            vcf@gt <- matrix(vcf@gt[probe,], byrow = TRUE, nrow = nr)
            colnames(vcf@gt) <- ng
        } else stop(" No variants after filter by position")
    }

    if (mode == "type" | mode == "both") {
        if (missing(gff))
            stop("gff is missing!")
        if (type == "custom")
            type <- cusTyp
        p <- type %in% unique(gff$type)
        m <- paste(unique(gff$type), collapse = "','")
        if (FALSE %in% p)
            stop("type should in c('",m,"')")
        if ("genome" %in% type) {
            gff <- gff
        } else {
            if(!missing(geneID)){
                ids <- tolower(gff$ID)
                nms <- tolower(gff$Name)
                id <- tolower(geneID)
                p <- (stringr::str_detect(ids,id)) | (stringr::str_detect(nms,id))
            }
            gff <- gff[(gff$type %in% type) & p]
        }


        POS <- vcfR::getPOS(vcf)
        POS <- as.numeric(POS)
        POSRange <- POS2GRanges(Chr = vcfR::getCHROM(vcf),
                                POS = POS)
        POSRange_rm <- POSRange[!(POSRange %over% gff)]

        POS_rm <- IRanges::start(POSRange_rm)
        probe <- !(POS %in% POS_rm)
        nr <- sum(probe)
        if(nr > 1){
            vcf@fix <- vcf@fix[probe,]
            vcf@gt <- vcf@gt[probe,]
        } else if(nr == 1){
            nf <- colnames(vcf@fix)
            ng <- colnames(vcf@gt)
            vcf@fix <- matrix(vcf@fix[probe,], byrow = TRUE, nrow = nr)
            colnames(vcf@fix) <- nf
            vcf@gt <- matrix(vcf@gt[probe,], byrow = TRUE, nrow = nr)
            colnames(vcf@gt) <- ng
        } else stop(" No variants after filter by ", type)
    }

    return(vcf)
}


#' @name filter_table
#' @title filter variants stored in table
#' @param x genotype dataset in hapmap format, object of data.frame class
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
#' @param geneID gene ID
#' @examples
#' # example
#' data("geneHapR_test")
#' table <- filter_table(gt.geno, mode = "POS",
#'                       Chr = "scaffold_1", start = 4100, end = 5000)
#' @export
filter_table <- function(x,
                         mode = c("POS", "type", "both"),
                         Chr = Chr, start = start, end = end,
                         gff = gff, type = type, cusTyp = cusTyp,
                         geneID = geneID){
    probe <- rep(TRUE, nrow(x))

    if(mode %in% c("POS", "both")){
        if(missing(start) | missing(end)) stop("start or end is missing")
        POS <- x[,2] %>% as.numeric()
        if(!missing(Chr))
            probe <- probe & as.character(x[, 1]) %in% as.character(Chr)
        if(!missing(start))
            probe <- probe & x[, 2] > min(start, end)
        if(!missing(end))
            probe <- probe & x[, 2] < max(start, end)
    }


    if(mode %in% c("type", "both")){
        if (missing(gff))
            stop("gff is missing!")
        p <- type %in% unique(gff$type)
        m <- paste(unique(gff$type), collapse = "','")
        if (!all(p))
            stop("type should in c('",m,"')")
        if (type == "custom")
            type <- cusTyp
        if ("genome" %in% type) {
            gff <- gff
        } else {
            if(!missing(geneID)){
                ids <- tolower(gff$ID)
                nms <- tolower(gff$Name)
                id <- tolower(geneID)
                p <- (stringr::str_detect(ids,id)) | (stringr::str_detect(nms,id))
            }
            gff <- gff[(gff$type %in% type) & p]
        }
        POSRange <- POS2GRanges(Chr = "scaffold_1", POS = as.numeric(x[, 2]))
        probe <- probe & POSRange %over% gff
    }
    if(! any(probe))
        warning("There is no overlaps")

    return(x[which(probe), ])
}


#' @name filter_hmp
#' @title filter variants in hapmap format
#' @param x genotype dataset in hapmap format, object of data.frame class
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
#' @param geneID gene ID
#' @examples
#' # create a dataset of genotype in hapmap format
#' hmp <- hap2hmp(hapResult);
#'
#' # example
#' hmp <- filter_hmp(hmp, mode = "POS",
#'                       Chr = "scaffold_1", start = 4100, end = 5000)
#'
#' @export
filter_hmp <- function(x,
                       mode = c("POS", "type", "both"),
                       Chr = Chr, start = start, end = end,
                       gff = gff, type = type, cusTyp = cusTyp,
                       geneID = geneID){
    probe <- rep(TRUE, nrow(x))

    if(mode %in% c("POS", "both")){
        if(!missing(Chr))
            probe <- probe & as.character(x[, 3]) %in% as.character(Chr)
        if(!missing(start))
            probe <- probe & x[, 4] > min(start, end)
        if(!missing(end))
            probe <- probe & x[, 4] < max(start, end)
    }


    if(mode %in% c("type", "both")){
        if (missing(gff))
            stop("gff is missing!")
        p <- type %in% unique(gff$type)
        m <- paste(unique(gff$type), collapse = "','")
        if (!all(p))
            stop("type should in c('",m,"')")
        if (type == "custom")
            type <- cusTyp
        if ("genome" %in% type) {
            gff <- gff
        } else {
            if(!missing(geneID)){
                ids <- tolower(gff$ID)
                nms <- tolower(gff$Name)
                id <- tolower(geneID)
                p <- (stringr::str_detect(ids,id)) | (stringr::str_detect(nms,id))
            }
            gff <- gff[(gff$type %in% type) & p]
        }
        POSRange <- POS2GRanges(Chr = "scaffold_1", POS = as.numeric(x[, 4]))
        probe <- probe & POSRange %over% gff
    }
    if(! any(probe))
        warning("There is no overlaps")

    return(x[which(probe), ])
}


#' @name filter_plink.pedmap
#' @title filter_plink.pedmap
#' @description used for filtration of p.link
#' @usage
#' filter_plink.pedmap(x,
#'                     mode = c("POS", "type", "both"),
#'                     Chr = Chr, start = start, end = end,
#'                     gff = gff, type = type, cusTyp = cusTyp,
#'                     geneID = geneID)
#' @param x a list stored the p.link information
#' @param mode filtration mode, one of c("POS", "type", "both")
#' @param Chr the chromosome name, need if mode set as POS or both
#' @param start,end numeric, the range of filtration, and the start should smaller than end
#' @param gff the imported gff object
#' @param type should be in `unique(gff$type)`, usually as "CDS", "genome".
#' @param cusTyp if `type` set as custom, then `cusTyp` is needed
#' @param geneID geneID
#' @inherit plink.pedmap2hap examples
#' @return list, similar with `x`, but filtered
#' @export
filter_plink.pedmap <- function(x,
                                mode = c("POS", "type", "both"),
                                Chr = Chr, start = start, end = end,
                                gff = gff, type = type, cusTyp = cusTyp,
                                geneID = geneID){
    on.exit(gc(verbose = FALSE))
    map <- x$map
    probe <- rep(TRUE, nrow(map))


    if(mode %in% c("POS", "both")){
        if(missing(start) | missing(end)) stop("start or end is missing")
        POS <- c(start, end)
        if(!missing(Chr))
            probe <- probe & as.character(map[, 1]) %in% as.character(Chr)
        if(!missing(start))
            probe <- probe & map[, 4] > min(POS)
        if(!missing(end))
            probe <- probe & map[, 4] < max(POS)
    }


    if(mode %in% c("type", "both")){
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
            if(!missing(geneID)){
                ids <- tolower(gff$ID)
                nms <- tolower(gff$Name)
                id <- tolower(geneID)
                p <- (stringr::str_detect(ids,id)) | (stringr::str_detect(nms,id))
            }
            gff <- gff[(gff$type %in% type) & p]
        }
        POSRange <- POS2GRanges(Chr = "scaffold_1", POS = map[, 4])
        probe <- probe & POSRange %over% gff
    }
    if(! any(probe))
        warning("There is no overlaps")

    probe_ped <- c(1:6, which(probe) * 2 + 6, which(probe) * 2 + 5)
    probe_ped <- probe_ped[order(probe_ped)]
    return(list(map = x$map[which(probe), ],
                ped = x$ped[, probe_ped]))
}


#' @name filter_hap
#' @title Filter hap
#' @description filter hapResult or hapSummary by remove
#'   positions or accessions or haplotypes
#' @usage
#' filter_hap(hap,
#'            rm.mode = c("position", "accession", "haplotype", "freq"),
#'            position.rm = position.rm,
#'            accession.rm = accession.rm,
#'            haplotype.rm = haplotype.rm,
#'            freq.min = 5)
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
#' hap <- filter_hap(hapResult,
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
                warning("accession.rm is missing")
            rm.probe <- hap$Accession %in% accession.rm
            rm.probe <- sapply(rm.probe, function(x) isTRUE(x))
            hap <- hap[which(!rm.probe),]
            AccRemoved <- c(AccRemoved, accession.rm)
        }
    }

    if("haplotype" %in% rm.mode){
        if(missing(haplotype.rm))
            warning("haplotype is missing")
        rm.probe <- hap$Hap %in% haplotype.rm
        rm.probe <- sapply(rm.probe, function(x) isTRUE(x))
        hap <- hap[which(! rm.probe),]
        AccRemoved <- c(AccRemoved, hap2acc[names(hap2acc) %in% haplotype.rm])
    }

    if("freq" %in% rm.mode){
        if(inherits(hap, "hapSummary")) {
            if(missing(freq.min))
                warning("freq.min is missing")
            rm.probe <- as.numeric(hap$freq) < freq.min
            rm.probe <- sapply(rm.probe, function(x) isTRUE(x))
            hap <- hap[which(! rm.probe),]
            AccRemoved <- c(AccRemoved, hap2acc[! names(hap2acc) %in% hap$Hap])
        } else {
            freq <- table(hap$Hap[-c(1:4)])
            rm.probe <- freq < freq.min
            rm.probe <- hap$Hap %in% names(freq)[which(rm.probe)]
            hap <- hap[which(!rm.probe),]
            # warning("only 'hapSummary' class surrport filtered by 'freq'")
            AccRemoved <- c(AccRemoved, hap$Accession[which(rm.probe)])
        }
    }


    if("position" %in% rm.mode){
        if(missing(position.rm))
            stop("position.rm is missing")

        rm.probe <- hap[hap$Hap == "POS",] %>%
            t() %>%
            as.vector()
        rm.probe <- rm.probe %in% as.character(position.rm)
        rm.probe <- sapply(rm.probe, function(x) isTRUE(x))
        hap <- hap[, which(! rm.probe)]
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
