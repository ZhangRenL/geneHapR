#' @name get_hap
#' @title generat haps from vcf
#' @description  generate hap format from vcf.
#' @usage get_hap(vcf,
#' hap_prefix = "H",
#' filter_Chr = FALSE,
#' Chr = "scaffold_1",
#' filter_POS = FALSE,
#' startPOS = 104,
#' endPOS = 9889,
#' hyb_remove = TRUE,
#' na.drop = TRUE)
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
#' @param hap_prefix hap_prefix, defalt is "H"
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
    vcf, hap_prefix = "H",
    filter_Chr = FALSE, Chr = "scaffold_1",
    filter_POS = FALSE, startPOS = 104, endPOS = 9889,
    hyb_remove = TRUE,
    na.drop = TRUE) {
    requireNamespace('tidyr')
    requireNamespace('dplyr')
    if(filter_Chr){
        if(missing(Chr)) stop("Chr must be character")
        vcfChr <- vcf@fix[,1]
        probe <- vcfChr %in% Chr
        vcf@fix <- vcf@fix[probe,]
        vcf@gt <- vcf@gt[probe,]
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
    }

    vcf <- order_vcf(vcf)
    CHR <- vcfR::getCHROM(vcf)
    POS <- vcfR::getPOS(vcf)
    REF <- vcfR::getREF(vcf)
    ALT <- vcfR::getALT(vcf)
    ALLELE <- paste0(REF,"/",ALT)
    INFO <- vcfR::getINFO(vcf)

    hap <- vcf2hap_data(vcf, REF = REF, ALT= ALT, ALLELE = ALLELE, POS = POS)
    hap <- reduce_hap(hap)

    # Drop hyb or N
    if(hyb_remove) {
        hap[hap == "H"] <- NA
        hap <- na.omit(hap)
    }
    if(na.drop) {
        hap[hap == "N"] <- NA
        hap <- na.omit(hap)
    }

    hap <- assign_hapID(hap, hap_prefix)

    # add infos
    meta <- rbind(
        c("CHR",CHR,""),
        c("POS",POS,""),
        c("INFO", INFO,""),
        c("ALLELE",ALLELE,""))
    colnames(meta) <- colnames(hap)
    hap <- rbind(meta, hap)
    rownames(hap) <- seq_len(nrow(hap))

    hap <- remove_redundancy_col(hap)
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
vcf2hap_data <- function(vcf, REF = REF, ALT= ALT, ALLELE = ALLELE, POS = POS){
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

    # convert Indel into +/-
    nr = nrow(hap)
    for(l in seq_len(nr)){
        if(stringr::str_length(ALLELE[l]) > 3){
            if(stringr::str_length(REF[l]) > stringr::str_length(ALT[l])){
                gety <- c("++", "+-", "-+", "--", "NN")
            } else {
                gety <- c("--", "-+", "+-", "++", "NN")
            }
            REFl <- REF[l]
            ALTl <- ALT[l]
            names(gety) <- paste(
                c(REFl,REFl,ALTl,ALTl,"N"),
                c(REFl,ALTl,REFl,ALTl,"N"),
                sep = "/")
            probe <- hap[l,]
            hap[l,] <- gety[probe]
        } else {
            hap[l,] <- stringr::str_remove_all(hap[l,], "[/]")
        }
    }

    hap <- t(hap)

    return(hap)
}


reduce_hap <- function(hap){
    # reform the genotypes
    # +/-,-/+ ->H
    # A/G,A/T,A/C -> H
    # C/A,C/G,C/T -> H
    # G/A,G/T,G/A -> H
    # T/A,T/G,T/C -> H
    # ++ -> +; -- -> -
    # A/A -> A; T/T ->T; C/C -> C; G/G ->G
    # N/N -> N
    probe_hyb <- c("AA","CC","GG","TT","++","--","NN")
    hap[!(hap %in% probe_hyb)] <- "H"

    for(i in probe_hyb) {
        hap[hap == i] <- stringr::str_sub(i,1,1)
    }
    return(hap)
}


assign_hapID <- function(hap, hap_prefix){
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
    hapnms <- paste0(hap_prefix, hapnms)
    names(hapnms) <- names(haps)
    hap$Hap <- hapnms[hap$Hap]
    hap <- hap[order(hap$Hap),]
    return(hap)
}


remove_redundancy_col <- function(hap){
    # removed Redundancy cols
    removecols <- c()
    nc_hap <- ncol(hap)
    for(c in seq_len(nc_hap)){
        namec <- colnames(hap)[c]
        if(!(namec %in% c("Hap", "Accession", "freq"))){
            gtc <- unique(hap[-c(1, 2, 3, 4),c])
            if(length(gtc) == 1) {
                removecols <- c(removecols, c)
            }
        }
    }
    if(!is.null(removecols)) hap <- hap[, -removecols]
    return(hap)
}

#' @name hap_result
#' @title generate hap results
#' @description summarize hap result and output a txt file
#' @usage hap_result(hap, hap_prefix = "H", out = FALSE, file = "hapResult.txt")
#' @examples
#'
#' data("quickHap_test")
#' # You can import your vcf data by `import_vcf()` here
#' hap <- get_hap(vcf, hyb_remove = TRUE, na.drop =TRUE)
#' hapResult <- hap_result(hap)
#' plotHapTable(hapResult)
#' @import tidyr
#' @importFrom utils write.table
#' @param hap hap
#' @param out write hap results to a txt file
#' @param file file path
#' @export
#' @return data.frame, first four rows are fixed to meta information: CHR, POS, INFO, ALLELE
#' Hap names were placed in col1, Accessions and freqs were placed at the last two cols.
hap_result <- function(hap, hap_prefix = "H", out = FALSE, file = "hapResult.txt"){
    requireNamespace('tidyr')
    hapResults <- hap %>% data.frame(check.names = FALSE)
    hapfre <- table(hapResults[,1])
    hapfre <- hapfre[stringr::str_starts(names(hapfre), hap_prefix)]
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
    return(hapResults)
}
