#' @name plink.pedmap2hap
#' @title plink.pedmap2hap
#' @description convert p.link format data into hapResult
#' @param p.link list contains p.link information
#' @param hapPrefix prefix of haplotype names
#' @param hyb_remove whether remove accessions contains hyb-sites
#' @param na.drop whether drop accessions contains missing data ("N", NA)
#' @usage
#'   plink.pedmap2hap(p.link,
#'                    hapPrefix = "H",
#'                    hyb_remove = TRUE,
#'                    na.drop = TRUE)
#' @examples
#' \donttest{
#'    pedfile <- system.file("extdata",
#'                           "snp3kvars-CHR8-25947258-25951166-plink.ped",
#'                           package = "geneHapR")
#'    mapfile <- system.file("extdata",
#'                           "snp3kvars-CHR8-25947258-25951166-plink.map",
#'                           package = "geneHapR")
#'    p.link <- import_plink.pedmap(pedfile = pedfile, mapfile = mapfile,
#'                                  sep_map = "\t", sep_ped = "\t")
#'    p.link <- filter_plink.pedmap(p.link, mode = "POS",
#'                                  Chr = "chr08", start = 25948004,
#'                                  end = 25949944)
#'    hapResult <- plink.pedmap2hap(p.link, hapPrefix = "H",
#'                                  hyb_remove = TRUE,
#'                                  na.drop = TRUE)
#' }
#' @return object of hapSummary class
#' @export
plink.pedmap2hap <- function(p.link,
                             hapPrefix = "H",
                             hyb_remove = TRUE,
                             na.drop = TRUE) {
    map <- p.link$map

    # set options
    options <- c(hapPrefix = hapPrefix)
    options <- c(options, CHROM = unique(map[,1]))
    options <- c(options, POS = paste0(range(map[,4]), collapse = "-"))

    POS <- map[,4]
    CHR <- map[,1]
    INFO <- map[,2]
    hapData <- plink.pedmap2hapdata(p.link, POS = POS)
    hap <- hapData$hap
    allS_new <- hapData$allS_new
    ALLELE <- hapData$alleles

    # Drop hyb or N
    if (hyb_remove) {
        hap[grepl("|", hap, fixed = T)] <- NA
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
                  c("ALLELE", ALLELE, "")) %>% data.frame()
    colnames(meta) <- colnames(hap)
    hap <- rbind(meta, hap)
    rownames(hap) <- seq_len(nrow(hap))



    # set attributes
    hap <- remove_redundancy_col(hap)
    class(hap) <- unique(c('hapResult', "data.frame"))
    accAll <- p.link$ped[, 1]
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


plink.pedmap2hapdata <- function(p.link, POS = POS){
    allS_new <- allS

    # prepare for hapdata
    map <- p.link$map
    ped <- as.matrix(p.link$ped)
    ped.cls <- ncol(ped)
    ped[ped == 0] <- "N"
    hap <- data.frame(row.names = ped[,1])
    alleles <- c()
    for(i in seq_len(ped.cls / 2 - 3)){
        al.i <- ped[, 5 + 2 * i]
        al.i1 <- ped[, 6 + 2 * i]

        # prepare genotype
        col.i <- p <- al.i == al.i1
        col.i[p] <- al.i[p]
        col.i[!p] <- paste(al.i[!p], al.i1[!p], sep = "|")

        # hapdata column
        hap <- cbind(hap, col.i)

        # prepare alleles
        alele.i <- c(al.i, al.i1)
        alele.i[alele.i == "N"] <- NA
        alele.i <- na.omit(alele.i) %>% unique()
        if(length(alele.i) > 2){
            alele.i[2] <- paste0(alele.i[1:2], collapse = "/")
            alele.i <- paste(alele.i[2:length(alele.i)], collapse = ",")
        } else if(length(alele.i) == 2){
            alele.i <- paste(alele.i, collapse = "/")
        } else if(length(alele.i) == 1){
            alele.i <- paste(alele.i, alele.i, sep = "/")
        }
        alleles <- c(alleles, alele.i)
    }
    colnames(hap) <- POS
    hap <- as.matrix(hap) %>% toupper()
    allS_new <- table(hap) %>% names()
    allS_new
    return(list(hap = hap, allS_new = allS_new, alleles = alleles))
}

