#' @name table2hap
#' @title table2hap
#' @description convert variants stored in table format into hapResult
#' @param x a data.frame contains variants information.
#' The first file column are fix as Chrome name, position, reference nuclieotide,
#' alter nuclieotide and INFO. Accession genotype should be in followed columns.
#' "-" will be treated as Indel. "." and "N" will be treated as missing data.
#' Heterozygotes should be "A/T", "AAA/A"
#' @param hapPrefix prefix of haplotype names
#' @param hetero_remove whether remove accessions contains hyb-sites, Character not A T C G
#' @param na_drop whether drop accessions contains missing data ("N", "NA", ".")
#' @param pad The number length in haplotype names should be extend to.
#' @examples
#' \donttest{
#'    data("geneHapR_test")
#'    hapResult <- table2hap(gt.geno, hapPrefix = "H",
#'                           hetero_remove = TRUE,
#'                           na_drop = TRUE)
#' }
#' @return object of hapSummary class
#' @export
table2hap <- function(x,
                      hapPrefix = "H",
                      pad = 3,
                      hetero_remove = TRUE,
                      na_drop = TRUE) {

    if(! inherits(x, "data.frame"))
        x <- try(as.data.frame(x))
    if(! inherits(x, "data.frame"))
        stop("x is not in data.frame and can't convert into a data.fram")
    if(! is.numeric(x[,2]))
        stop("second column (position) should be numeric")

    # set options
    options <- c(hapPrefix = hapPrefix)
    options <- c(options, CHROM = unique(x[,1]))
    options <- c(options, POS = paste0(range(x[,2]), collapse = "-"))


    CHR <- x[,1]
    POS <- x[,2]
    REF <- x[,3]
    ALT <- x[,4]
    ALLELE <- paste(REF, ALT, sep = "/")
    INFO <- x[,5]

    hap <- table2hapdata(x, POS = POS)
    allS_new <- unique(ALLELE)

    # Drop hyb or N
    if (hetero_remove) {
        hap[grepl("|", hap, fixed = T)] <- NA
        hap[grepl("/", hap, fixed = T)] <- NA
        hap <- na.omit(hap)
        options <- c(options, hetero_remove = "YES")
    } else
        options <- c(options, hetero_remove = "NO")

    if (na_drop) {
        hap[hap == "N"] <- NA
        hap <- na.omit(hap)
        options <- c(options, NA_remove = "YES")
    } else
        options <- c(options, NA_remove = "NO")

    hap <- assign_hapID(hap, hapPrefix, pad)

    # add infos
    meta <- rbind(c("CHR", CHR, ""),
                  c("POS", POS, ""),
                  c("INFO", INFO, ""),
                  c("ALLELE", ALLELE, "")) %>% data.frame()
    colnames(meta) <- colnames(hap)
    hap <- rbind(meta, hap)
    rownames(hap) <- seq_len(nrow(hap))


    # set attributes
    hap <- remove_redundancy_col(hap)
    class(hap) <- unique(c('hapResult', "data.frame"))
    accAll <- colnames(x)[-c(1:4)]
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


table2hapdata <- function(x, POS = POS){
    rownames(x) <- POS
    hap <-  x %>%
        t() %>%
        toupper()
    hap <- hap[-c(1 : 5),]

    # deal with missing data
    hap[hap == "."] <- "N"
    hap[hap == "NA"] <- "N"
    hap[is.na(hap)] <- "N"

    p <- c("A|G", "C|T", "A|C", "G|T", "G|C", "A|T",
           "A|T|C", "G|T|C", "G|A|C", "G|A|T")
    names(p) <- c("R", "Y", "M", "K", "S", "W",
                  "H", "B", "V", "D")
    hap[hap %in% names(p)] <- p[hap[hap %in% names(p)]]
    hap <- gsub(pattern = "|", replacement = "/", x = hap, fixed = TRUE)
    p1 <- c("A", "T", "C", "G")
    names(p1) <- c("A|A","T|T","C|C","G|G")
    hap[hap %in% names(p1)] <- p1[hap[hap %in% names(p1)]]

    # prepare hapdata
    hap <- data.frame(hap, check.names = FALSE)
    hap$Accession <- rownames(hap)
    return(as.matrix(hap))
}

