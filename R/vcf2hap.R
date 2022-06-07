#' @name hap_result
#' @title generate hap results
#' @description summarize hap result and output a txt file
#' @usage hap_result(hap, hapPrefix = "H", out = FALSE, file = "hapResult.txt")
#' @examples
#' data("quickHap_test")
#' # Your vcf data should imported by `import_vcf()` here
#' hap <- get_hap(vcf, hyb_remove = TRUE, na.drop =TRUE)
#' hapResult <- hap_result(hap)
#' plotHapTable(hapResult)
#' @import tidyr
#' @importFrom utils write.table
#' @param hap hap
#' @param hapPrefix prefix of hap names
#' @param out write hap results to a txt file
#' @param file file path
#' @export
#' @return data.frame, first four rows are fixed to meta information: CHR, POS, INFO, ALLELE
#' Hap names were placed in col1, Accessions and freqs were placed at the last two cols.
hap_result <- function(hap, hapPrefix = "H", out = FALSE, file = "hapResult.txt"){
    if(!inherits(hap,"haptypes"))
        stop("
'hap' must be of class 'haptypes',
you should generate haps with function 'get_hap' first")
    requireNamespace('tidyr')
    hapResults <- hap %>% data.frame(check.names = FALSE)
    hapfre <- table(hapResults[,1])
    hapfre <- hapfre[stringr::str_starts(names(hapfre), hapPrefix)]
    hapResults <- hapResults %>% tidyr::chop(cols = "Accession")
    prob = hapResults[5:nrow(hapResults),1]
    hapResults$freq[5:nrow(hapResults)] <- hapfre[prob]
    Acc <- c()
    nAcc = length(hapResults$Accession)
    for(i in seq_len(nAcc)) {
        Acc[i] <- paste(hapResults$Accession[[i]],collapse = ";")
    }
    hapResults$Accession <- Acc

    if(out)  {
        utils::write.table(
            hapResults,
            file = file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
    }
    class(hapResults) <- c("hapSummary", "data.frame")
    attr(hapResults, "options") <- attr(hap, "options")
    attr(hapResults, "hap2acc") <- attr(hap, "hap2acc")

    return(hapResults)
}


#' @name get_hap
#' @title generat haps from vcf
#' @description  generate hap format from vcf.
#' @usage
#' get_hap(vcf,
#'         hapPrefix = "H",
#'         filter_Chr = FALSE,
#'         Chr = "scaffold_1",
#'         filter_POS = FALSE,
#'         startPOS = 104,
#'         endPOS = 9889,
#'         hyb_remove = TRUE,
#'         na.drop = TRUE)
#' @author Zhangrenl
#' @examples
#'
#' data("quickHap_test")
#' # You can import your vcf data by `import_vcf()` here
#' hap <- get_hap(vcf,hyb_remove = TRUE, na.drop =TRUE)
#' hapResult <- hap_result(hap)
#' plotHapTable(hapResult)
#' @import tidyr
#' @import vcfR
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom  stats na.omit
#' @param vcf vcf
#' @param filter_Chr filter vcf by Chrome or not, defalt is FALSE.
#'  If TRUE, the Chr is needed.
#' @param filter_POS filter vcf by Position or not, defalt is FALSE.
#'  If TRUE, startPOS and endPOS are needed.
#' @param hapPrefix hapPrefix, defalt is "H"
#' @param hyb_remove Remove accessions contains hybrid site or not.
#' Defalt is TRUE.
#' @param na.drop Drop Accessions contains unknown allele site or note.
#' Defalt is TRUE
#' @param Chr Needed for filter vcf by Chrom
#' @param startPOS,endPOS Needed for filter vcf by position.
#' startPOS must less than endPOS
#' @export
#' @return data.frame, first four rows are fixed to meta information: CHR, POS, INFO, ALLELE
#' Hap names were placed in col1, Accessions were placed at the last col.
get_hap <- function(
    vcf, hapPrefix = "H",
    filter_Chr = FALSE, Chr = "scaffold_1",
    filter_POS = FALSE, startPOS = 104, endPOS = 9889,
    hyb_remove = TRUE,
    na.drop = TRUE) {
    requireNamespace('tidyr')
    requireNamespace('dplyr')
    allS_new <- allS
    options <- c(hapPrefix = hapPrefix)
    if(filter_Chr){
        if(missing(Chr)) stop("Chr must be character")
        vcfChr <- vcf@fix[,1]
        probe <- vcfChr %in% Chr
        vcf@fix <- vcf@fix[probe,]
        vcf@gt <- vcf@gt[probe,]
        options <- c(options, CHROM = Chr)
    } else {
        Chr <- unique(vcfR::getCHROM(vcf))
        if(length(Chr) > 1) stop("only one CHROM should be in vcf, consider set 'filter_Chr' as 'TRUE'")
        options <- c(options, CHROM = Chr)
    }

    # filter Postion according given range
    if(filter_POS){
        if(missing(startPOS)) stop("startPOS must be numeric")
        if(missing(endPOS)) stop("endPOS must be numeric")
        if(startPOS >= endPOS) stop("startPOS must less tan endPOS")
        vcfPOS <- as.numeric(vcf@fix[,2])
        probe <- c(vcfPOS >= startPOS & vcfPOS <= endPOS)
        if(!(TRUE %in% probe)) {
            e = paste0(
                "There is no variant in selected range.
                Please check vcf file between ",
                startPOS," and ", endPOS, ".")
            return(e)
        }
        vcf@fix <- vcf@fix[probe,]
        vcf@gt <- vcf@gt[probe,]
        options <- c(options, POS = paste0(startPOS,"-",endPOS))
    } else {
        POS <- vcfR::getPOS(vcf)
        options <- c(options, POS = paste0(min(POS),"-",max(POS)))
    }

    # extract information from vcf
    vcf <- order_vcf(vcf)
    CHR <- vcfR::getCHROM(vcf)
    POS <- vcfR::getPOS(vcf)
    REF <- vcfR::getREF(vcf)
    ALT <- vcfR::getALT(vcf)
    ALLELE <- paste0(REF,"/",ALT)
    INFO <- vcfR::getINFO(vcf)

    # generate hap data from vcf
    hapData <- vcf2hap_data(vcf, allS_new = allS_new,
                        REF = REF, ALT= ALT, ALLELE = ALLELE, POS = POS)
    allS_new <- hapData$allS_new
    hap <- hapData$hap
    # Drop hyb or N
    if(hyb_remove) {
        hap[!hap %in% allS_new$homo] <- NA
        hap <- na.omit(hap)
        options <- c(options, hyb_remove = "YES")
    } else options <- c(options, hyb_remove = "NO")

    if(na.drop) {
        hap[hap == "N"] <- NA
        hap <- na.omit(hap)
        options <- c(options, NA_remove = "YES")
    } else options <- c(options, NA_remove = "NO")

    hap <- assign_hapID(hap, hapPrefix)

    # add infos
    meta <- rbind(
        c("CHR", CHR, ""),
        c("POS", POS, ""),
        c("INFO", INFO, ""),
        c("ALLELE", ALLELE, ""))
    colnames(meta) <- colnames(hap)
    hap <- rbind(meta, hap)
    rownames(hap) <- seq_len(nrow(hap))

    # set attributes
    hap <- remove_redundancy_col(hap)
    class(hap) <- unique(c("haptypes", "data.frame"))
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


order_vcf <- function(vcf){
    POS <- vcfR::getPOS(vcf)
    probe <- order(POS)
    vcf@fix <- vcf@fix[probe,]
    vcf@gt <- vcf@gt[probe,]
    return(vcf)
}


#' @importFrom rlang .data
vcf2hap_data <- function(vcf, allS_new = allS_new,
                         REF = REF, ALT= ALT, ALLELE = ALLELE, POS = POS){
    # vcf2data.frame for analysis
    gt <- vcfR::extract_gt_tidy(vcf)
    hap <- tidyr::pivot_wider(
        data = gt,
        id_cols = .data$Key,
        names_from = .data$Indiv,
        values_from = .data$gt_GT_alleles)
    hap <- dplyr::select(hap, -c(.data$Key))
    hap <- as.matrix(hap)
    rownames(hap) <- POS

    # convert "." into "N/N"
    hap[hap == "."] <- "N/N"

    # convert Indel(biallelic site) into +/-
    nr = nrow(hap)
    indel_probe <- vcfR::is.indel(vcf)
    biallelic_probe <- vcfR::is.biallelic(vcf)
    for(l in seq_len(nr)){
        if(biallelic_probe[l]){
            if(indel_probe[l])
                allS_new <- update_allS(allS_new, REF = REF[l], ALT = ALT[l])
        } else
            allS_new <- update_allS(allS_new, REF = REF[l], ALT = ALT[l])
    }

    hap <- t(hap)

    # reform the genotypes
    # homo site convert into single
    # ++ -> +; -- -> -
    # A/A -> A; T/T ->T; C/C -> C; G/G ->G
    # N/N -> N
    # hetero convert "/" to "|"
    for(i in seq_len(ncol(hap))) {
        hap[,i] <- allS_new$all[hap[,i]]
    }
    return(list(hap = hap, allS_new = allS_new))
}





assign_hapID <- function(hap, hapPrefix){
    # name haps
    hap <- data.frame(hap, check.rows = FALSE, check.names = FALSE)

    HapID <- tidyr::unite(
        hap,
        dplyr::matches("[0-9]{1,}"),
        col = "IDs",
        sep = "")
    HapID <- HapID$IDs
    hap <- cbind(Hap = HapID, hap)
    if(!"Accession" %in% colnames(hap))
        hap$Accession <- row.names(hap)
    haps <- table(hap$Hap)
    haps <- haps[order(haps,decreasing = TRUE)]
    n_hs <- length(haps)
    hapnms <- stringr::str_pad(seq_len(n_hs),3,"left","0")
    hapnms <- paste0(hapPrefix, hapnms)
    names(hapnms) <- names(haps)
    hap$Hap <- hapnms[hap$Hap]
    hap <- hap[order(hap$Hap),]
    return(hap)
}


remove_redundancy_col <- function(hap){
    # removed Redundancy cols
    removecols <- c()
    nc_hap <- ncol(hap)
    if(hap[1,1] == "CHR")
        gth <- hap[-c(1, 2, 3, 4),] else
            gth <- hap
    for(c in seq_len(nc_hap)){
        namec <- colnames(hap)[c]
        if(!(namec %in% c("Hap", "Accession", "freq"))){
            gtc <- unique(gth[,c])
            if(length(gtc) == 1) {
                removecols <- c(removecols, c)
            }
        }
    }
    if(!is.null(removecols)) hap <- hap[, -removecols]
    return(hap)
}


