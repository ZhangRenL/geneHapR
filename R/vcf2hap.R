#' @name vcf2hap
#' @title Generat Haps from VCF
#' @description
#' Generate hapResult from vcfR object
#' A simple filter by position was provided in this function,
#' however it's prefer to filter VCF (vcfR object) through
#' \code{\link[geneHapR:filter_vcf]{filter_vcf()}}.
#' @usage
#' vcf2hap(vcf,
#'         hapPrefix = "H",
#'         filter_Chr = FALSE,
#'         filter_POS = FALSE,
#'         hyb_remove = TRUE,
#'         na.drop = TRUE, ...)
#' @author Zhangrenl
#' @examples
#' data("geneHapR_test")
#' hapResult <- vcf2hap(vcf)
#' @param vcf vcfR object imported by `import_vcf()`
#' @param filter_Chr logical, whether filter vcf by chromosome or not. Default
#' as `FALSE`. If set as `TRUE`, `Chr` is needed
#' @param filter_POS logical, whether filter vcf by position or not. Default
#' as `FALSE`. If set as `TRUE`, `startPOS` and `endPOS` are needed
#' @param hapPrefix prefix of hap names, default as "H"
#' @param hyb_remove whether remove accessions contains hybrid site or not.
#' Default as `TRUE`
#' @param na.drop whether remove accessions contains unknown allele site or not
#' Default as `TRUE`.
#' @param Chr Chromosome name, needed when `filter_Chr` was set as `TRUE`
#' @param startPOS,endPOS start and end position, needed when `filter_POS` was
#' set as `TRUE`. In addition, `startPOS` must less than `endPOS`
#' @param ... Parameters not used
#' @seealso
#' extract genotype from vcf:
#' \code{\link[vcfR:extract_gt_tidy]{vcfR::extract_gt_tidy()}},
#' import vcf files:
#' \code{\link[geneHapR:import_vcf]{import_vcf()}} (preferred) and
#' \code{\link[vcfR:read.vcfR]{vcfR::read.vcfR()}},
#' filter vcf according **position** and **annotations**:
#' \code{\link[geneHapR:filter_vcf]{filter_vcf()}}
#' @return
#' object of hapResult class
#' @import tidyr
#' @import vcfR
#' @importFrom rlang .data
#' @importFrom stats na.omit
#' @export
vcf2hap <- function(vcf,
                    hapPrefix = "H",
                    filter_Chr = FALSE,
                    # Chr = Chr,
                    filter_POS = FALSE,
                    # startPOS = startPOS,
                    # endPOS = endPOS,
                    hyb_remove = TRUE,
                    na.drop = TRUE, ...) {
    requireNamespace('tidyr')
    allS_new <- allS
    options <- c(hapPrefix = hapPrefix)
    if (filter_Chr) {
        message("Filter VCF by chromosome has been detached from vcf2hap")
    }
    # else {
    #         if (length(Chr) > 1)
    #             stop("
    # only one CHROM should be in vcf, consider set 'filter_Chr' as 'TRUE'
    #                  ")
    # }

    # filter Postion according given range
    if (filter_POS) {
        message("Filter VCF by position has been detached from vcf2hap")
    }

    # extract information from vcf
    CHR <- vcfR::getCHROM(vcf)
    POS <- vcfR::getPOS(vcf)
    options <- c(options, CHROM = CHR)
    options <- c(options, POS = paste0(min(POS), "-", max(POS)))
    REF <- vcfR::getREF(vcf)
    ALT <- vcfR::getALT(vcf)
    ALLELE <- paste0(REF, "/", ALT)
    INFO <- vcfR::getINFO(vcf)

    # generate hap data from vcf
    hapData <- vcf2hap_data(
        vcf,
        allS_new = allS_new,
        REF = REF,
        ALT = ALT,
        ALLELE = ALLELE,
        POS = POS
    )
    allS_new <- hapData$allS_new
    hap <- hapData$hap
    # Drop hyb or N
    if (hyb_remove) {
        hap[grepl("|", hap, fixed = T)] <- NA
        hap[grepl("/", hap, fixed = T)] <- NA
        hap <- na.omit(hap)
        options <- c(options, hyb_remove = "YES")
    } else
        options <- c(options, hyb_remove = "NO")

    if (na.drop) {
        hap[hap == "N"] <- NA
        hap <- na.omit(hap)
        options <- c(options, NA_remove = "YES")
    } else
        options <- c(options, NA_remove = "NO")

    hap <- assign_hapID(hap, hapPrefix)

    # add infos
    meta <- rbind(c("CHR", CHR, ""),
                  c("POS", POS, ""),
                  c("INFO", INFO, ""),
                  c("ALLELE", ALLELE, ""))
    colnames(meta) <- colnames(hap)
    hap <- rbind(meta, hap)
    rownames(hap) <- seq_len(nrow(hap))

    # set attributes
    hap <- remove_redundancy_col(hap)
    class(hap) <- unique(c('hapResult', "data.frame"))
    accAll <- colnames(vcf@gt)[-1]
    attr(hap, "AccAll") <- accAll
    accRemain <- hap$Accession[hap$Accession != ""]
    attr(hap, "AccRemain") <- accRemain
    attr(hap, "AccRemoved") <- accAll[!accAll %in% accRemain]
    attr(hap, "options") <- options
    hap2acc <- hap$Accession[-c(1:4)]
    names(hap2acc) <- hap$Hap[-c(1:4)]
    attr(hap, "hap2acc") <- hap2acc
    attr(hap, "freq") <- table(hap$Hap[-c(1:4)])

    return(hap)
}



# return: data.frame, individuals in rows and positions in cols
# with additional attrs
# eg. # $hap
#     #    41  136
#     # a  A   A|T
#     # b  G   G
#     # $allS_new
vcf2hap_data <- function(vcf,
                         allS_new = allS_new,
                         REF = REF,
                         ALT = ALT,
                         ALLELE = ALLELE,
                         POS = POS) {
    # vcf2data.frame for analysis
    gt <- vcfR::extract_gt_tidy(vcf)
    hap <- tidyr::pivot_wider(
        data = gt,
        id_cols = .data$Key,
        names_from = .data$Indiv,
        values_from = .data$gt_GT_alleles
    )
    hap <- hap[, colnames(hap) != "Key"]
    hap <- as.matrix(hap)
    rownames(hap) <- POS

    # convert "." into "N/N"
    hap[hap == "."] <- "N/N"
    hap[is.na(hap)] <- "N/N"

    # convert Indel(biallelic site) into +/-
    nr = nrow(hap)
    indel_probe <- vcfR::is.indel(vcf)
    biallelic_probe <- vcfR::is.biallelic(vcf)
    for (l in seq_len(nr)) {
        if (biallelic_probe[l]) {
            if (indel_probe[l])
                allS_new <-
                    update_allS(allS_new, REF = REF[l], ALT = ALT[l])
        } else
            allS_new <-
                update_allS(allS_new, REF = REF[l], ALT = ALT[l])
    }

    hap <- gsub(pattern = "|", replacement = "/", x = hap, fixed = TRUE) %>%
        t() %>%
        toupper()

    # reform the genotypes
    # homo site convert into single
    # ++ -> +; -- -> -
    # A/A -> A; T/T ->T; C/C -> C; G/G ->G
    # N/N -> N
    # hetero convert "/" to "|"
    for (i in seq_len(ncol(hap))) {
        hap[, i] <- allS_new$all[hap[, i]]
    }
    return(list(hap = hap, allS_new = allS_new))
}





assign_hapID <- function(hap, hapPrefix) {
    # name haps
    hap <- data.frame(hap, check.rows = FALSE, check.names = FALSE)

    HapID <- tidyr::unite(hap,
                          tidyr::matches("[0-9]{1, }"),
                          col = "IDs",
                          sep = "")
    HapID <- HapID$IDs
    hap <- cbind(Hap = HapID, hap)
    if (!"Accession" %in% colnames(hap))
        hap$Accession <- row.names(hap)
    haps <- table(hap$Hap)
    haps <- haps[order(haps, decreasing = TRUE)]
    n_hs <- length(haps)
    hapnms <- stringr::str_pad(seq_len(n_hs), 3, "left", "0")
    hapnms <- paste0(hapPrefix, hapnms)
    names(hapnms) <- names(haps)
    hap$Hap <- hapnms[hap$Hap]
    hap <- hap[order(hap$Hap),]
    return(hap)
}


remove_redundancy_col <- function(hap) {
    # removed Redundancy cols
    removecols <- c()
    nc_hap <- ncol(hap)
    if (hap[1, 1] == "CHR")
        gth <- hap[-c(1, 2, 3, 4),]
    else
        gth <- hap
    for (c in seq_len(nc_hap)) {
        namec <- colnames(hap)[c]
        if (!(namec %in% c("Hap", "Accession", "freq"))) {
            gtc <- unique(gth[, c])
            if (length(gtc) == 1) {
                removecols <- c(removecols, c)
            }
        }
    }
    if (!is.null(removecols))
        hap <- hap[,-removecols]
    return(hap)
}

