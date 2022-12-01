#' @name seqs2hap
#' @title Generate Hap Results from Seqs
#' @description generate hapResults from aligned and trimed sequences
#' @usage
#' seqs2hap(seqs,
#'          Ref = names(seqs)[1],
#'          hetero_remove = TRUE, na_drop = TRUE,
#'          maxGapsPerSeq = 0.25,
#'          hapPrefix = "H", ...)
#' @examples
#' \donttest{
#' data("geneHapR_test")
#' seqs
#' seqs <- trimSeqs(seqs,
#'          minFlankFraction = 0.1)
#' seqs
#' hapResult <- seqs2hap(seqs,
#'                       Ref = names(seqs)[1],
#'                       hetero_remove = TRUE, na_drop = TRUE,
#'                       maxGapsPerSeq = 0.25,
#'                       hapPrefix = "H")
#' }
#' @importFrom Biostrings letterFrequency width
#' @importFrom methods as
#' @param seqs object of DNAStringSet or DNAMultipleAlignment class
#' @param Ref the name of reference sequences.
#' Default as the name of the first sequence
#' @param hapPrefix prefix of hap names. Default as "H"
#' @param maxGapsPerSeq value in `[0, 1]` that indicates the maximum
#' fraction of gaps allowed in each seq after alignment (default as 0.25).
#' Seqs with gap percent exceed that will be dropped
#' @param hetero_remove whether remove accessions contains hybrid site or not.
#' Default as `TRUE`
#' @param na_drop whether drop sequeces contain "N"
#' Default as `TRUE`.
#' @param ... Parameters not used.
#' @inherit hap_summary examples
#' @return
#' object of hapResult class
#' @export
seqs2hap <- function(seqs,
                     Ref = names(seqs)[1],
                     hetero_remove = TRUE, na_drop = TRUE,
                     maxGapsPerSeq = 0.25,
                     hapPrefix = "H",
                     ...) {
    seqs <- as(seqs, "DNAStringSet")
    options <- c(hapPrefix = hapPrefix)
    RefSeq <- seqs[Ref]
    accAll <- names(seqs)
    if (names(seqs)[1] != Ref)
        seqs <- seqs[c(Ref, setdiff(names(seqs), Ref))]
    allS_new <- allS


    # remove seqs contains too much gaps
    rowGaps <- Biostrings::letterFrequency(seqs, "-")[, 1]
    nc <- Biostrings::width(seqs)[1]
    probe <- rowGaps / nc < maxGapsPerSeq
    if (FALSE %in% probe) {
        message(table(probe)["FALSE"],
                " seq(s) will be removed due to too much gaps")
        message("  ", paste(names(seqs)[!probe], collapse = ", "))
        seqs <- seqs[probe]
    }

    # seq2hapData
    hapData <-
        seq2hap_data(seqs, allS_new = allS_new, RefSeq = RefSeq)
    allS_new <- hapData$allS_new
    hap <- hapData$hap

    # Drop hyb or N
    if (hetero_remove) {
        hap[!hap %in% allS_new$homo] <- NA
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

    # set meta info
    POS <- colnames(hap)
    ALLELE <-
        apply(hap, 2, function(x)
            paste(unique(x), collapse = ","))
    ALLELE <- stringr::str_replace(ALLELE, ",", "/")
    INFO <- CHR <- rep(NA, ncol(hap))

    # assign hapID
    hap <- assign_hapID(hap, hapPrefix)
    hap <- rbind(
        CHR = c("CHR", CHR, ""),
        POS = c("POS", POS, ""),
        INFO = c("INFO", INFO, ""),
        ALLELE = c("ALLELE", ALLELE, ""),
        hap
    )

    hap <- remove_redundancy_col(hap)

    # set attributes
    class(hap) <- unique(c('hapResult', "data.frame"))
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


#' #' @name seqs2hap
#' #' @description allign imported sequences
#' #' @usage allignSeqs(seqs, ...)
#' #' @param ... parameters will pass to `muscle::muscle()`
#' #' @seealso
#' #' \code{\link[muscle:muscle]{muscle::muscle()}}
#' #' @export
#' allignSeqs <- function(seqs, ...) {
#'     muscle::muscle(seqs, ...)
#' }



#' @name seqs2hap
#' @usage
#' trimSeqs(seqs,
#'          minFlankFraction = 0.1)
#' @importFrom Biostrings as.matrix
#' @param minFlankFraction A value in `[0, 1]` that indicates the minimum
#' fraction needed to call a gap in the consensus string (default as 0.1).
#' @export
trimSeqs <- function(seqs,
                     minFlankFraction = 0.1) {
    mseqs <- Biostrings::as.matrix(seqs)
    nr <- nrow(mseqs)
    gaps <- apply(mseqs, 2, function(x)
        table(x)["-"])
    gaps <- gaps / nr

    trimLeft <- 0
    for (i in seq_len(length(gaps))) {
        if (gaps[i] <= minFlankFraction) {
            trimLeft <- i
            break
        }
    }

    trimRight <- 0
    for (j in rev(seq_len(length(gaps)))) {
        if (gaps[j] <= minFlankFraction) {
            trimRight <- j
            break
        }
    }

    if (trimLeft >= trimRight)
        stop(
            "Too many gaps found in flank of aligned seqs.
Please consider elevate 'minFlankFraction' OR align and trim seqs manually"
        )
    seqs <- as(seqs, "DNAStringSet")
    seqs <- Biostrings::subseq(seqs, trimLeft, trimRight)
    return(seqs)
}


# return data.frame,
seq2hap_data <- function(seqs,
                         allS_new = allS_new,
                         RefSeq = RefSeq) {
    # get refSeq trimed length
    Ref <- names(RefSeq)
    aRefSeq <- as(seqs, "DNAStringSet")[Ref]
    aRefSeq <- stringr::str_remove_all(aRefSeq, "-")
    range <- stringr::str_locate(as.character(RefSeq), aRefSeq)
    trimedLen <- range[1, "start"] - 1

    # convert DNAset into matrix
    RefSeq <- as.matrix(RefSeq)
    mseqs <- as.matrix(seqs)
    mRef <- mseqs[Ref, ]

    # deal with indels
    newmseqs <-
        matrix(
            nrow = nrow(mseqs),
            ncol = 0,
            dimnames = list(rownames(mseqs))
        )
    POS <- c()
    n <- 1
    for (i in seq_len(ncol(mseqs))) {
        if (!"-" %in% mseqs[, i]) {
            newmseqs <- cbind(newmseqs, mseqs[, i])
            POS <- c(POS, n)
        } else {
            newmseqs[, ncol(newmseqs)] = paste0(newmseqs[, ncol(newmseqs)], mseqs[, i])
        }
        if (mRef[i] != "-")
            n <- n + 1
    }
    colnames(newmseqs) <- POS

    # remove redudant cols
    probe <- apply(newmseqs, 2, function(x)
        length(unique(x)))
    newmseqs <- newmseqs[, probe != 1]

    # remove "-"
    newmseqs[, ] <- stringr::str_remove_all(newmseqs, "-")

    # update allSite information
    ALLELE <-
        apply(newmseqs, 2, function(x)
            paste(unique(x), collapse = ","))
    ALLELE <- stringr::str_replace(ALLELE, ",", "/")
    for (i in unique(ALLELE)) {
        ALi <- unlist(strsplit(i, "/"))
        allS_new <- update_allS(allS_new = allS_new,
                                REF = ALi[1],
                                ALT = ALi[2])
    }

    newmseqs <- t(newmseqs) %>%
        toupper() %>%
        t()
    return(list(hap = newmseqs, allS_new = allS_new))
}
