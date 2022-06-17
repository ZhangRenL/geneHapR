#' @name seqs2hap
#' @title generate hap results from seqs
#' @description generate hap results from seqs
#' @usage
#' seqs2hap(seqs,
#'          Ref = names(seqs)[1],
#'          hyb_remove = TRUE, na.drop = TRUE,
#'          needAlign = TRUE, needTrim = TRUE,
#'          minFlankFraction = 0.05, maxGapsPerSeq = 0.2,
#'          hapPrefix = "H", ...)
#' @importFrom muscle muscle
#' @importFrom Biostrings letterFrequency width
#' @importFrom methods as
#' @param seqs DNAStringSet or DNAMultipleAlignment
#' @param needAlign T
#' @param needTrim T
#' @param hapPrefix T
#' @param Ref T
#' @param hyb_remove T
#' @param na.drop T
#' @param minFlankFraction A value in `[0, 1]` that indicates the minimum fraction needed to call a gap in the consensus string (default is 0.1).
#' @param maxGapsPerSeq A value in `[0, 1]` that indicates the maximum gaps allowed in a seq after alignment
#' @param ... parameters will pass to musle::musle() for DNA multi alignment
#' @inherit hap_result examples
#' @export
seqs2hap <- function(seqs,
                     Ref = names(seqs)[1],
                     hyb_remove = TRUE, na.drop = TRUE,
                     needAlign = TRUE, needTrim = TRUE,
                     minFlankFraction = 0.05, maxGapsPerSeq = 0.2,
                     hapPrefix = "H", ...){
    seqs <- as(seqs, "DNAStringSet")
    options <- c(hapPrefix = hapPrefix)
    RefSeq <- seqs[Ref]
    accAll <- names(seqs)
    if(names(seqs)[1] != Ref) seqs <- seqs[c(Ref, setdiff(names(seqs), Ref))]
    allS_new <- allS

    # perform Multialignment
    if(needAlign) seqs <- allignSeqs(seqs) else {
        cat("Multialignment will not perform due to 'needAlign' is FALSE\n")
        seqs <- Biostrings::DNAMultipleAlignment(seqs)
    }

    # trim flanks
    if(needTrim){
        seqs <- trimSeqs(seqs, minFlankFraction = minFlankFraction)
    } else {
        cat("Seqs will not be trimed due to 'needTrim' is FALSE\n")
    }

    # remove seqs contains too much gaps
    rowGaps <- Biostrings::letterFrequency(seqs,"-")[,1]
    nc <- Biostrings::width(seqs)[1]
    probe <- rowGaps/nc < 0.1
    if(FALSE %in% probe){
        message(table(probe)["FALSE"], " seq(s) will be removed due to too much gaps")
        message("  ", paste(names(seqs)[!probe], collapse = ", "))
        seqs <- seqs[probe]
    }

    # seq2hapData
    hapData <- seq2hap_data(seqs, allS_new = allS_new, RefSeq = RefSeq)
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

    # set meta info
    POS <- colnames(hap)
    ALLELE <- apply(hap, 2, function(x) paste(unique(x), collapse = ","))
    ALLELE <- stringr::str_replace(ALLELE,",","/")
    INFO <- CHR <- rep(NA, ncol(hap))

    # assign hapID
    hap <- assign_hapID(hap, hapPrefix)
    hap <- rbind(CHR = c("CHR", CHR, ""),
                 POS = c("POS", POS, ""),
                 INFO = c("INFO", INFO, ""),
                 ALLELE = c("ALLELE", ALLELE, ""),
                 hap)

    hap <- remove_redundancy_col(hap)

    # set attributes
    class(hap) <- unique(c("haptypes", "data.frame"))
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

allignSeqs <- function(seqs, ...){
    muscle::muscle(seqs, ...)
}

#' @importFrom Biostrings as.matrix
trimSeqs <- function(seqs, minFlankFraction = minFlankFraction){
    mseqs <- Biostrings::as.matrix(seqs)
    nr <- nrow(mseqs)
    gaps <- apply(mseqs, 2, function(x) table(x)["-"])
    gaps <- gaps/nr

    trimLeft <- 0
    for(i in seq_len(length(gaps))){
        if(gaps[i] <= minFlankFraction) {
            trimLeft <- i
            break
        }
    }

    trimRight <- 0
    for(j in rev(seq_len(length(gaps)))){
        if(gaps[j] <= minFlankFraction) {
            trimRight <- j
            break
        }
    }

    if(trimLeft == trimRight | trimLeft == 0 | trimRight == 0)
        stop("Too many gaps found in flank of aligned seqs.
Please consider elevate 'minFlankFraction'")
    seqs <- as(seqs, "DNAStringSet")
    seqs <- Biostrings::subseq(seqs, trimLeft, trimRight)
    return(seqs)
}


# return data.frame,
seq2hap_data <- function(seqs, allS_new = allS_new, RefSeq = RefSeq){
    # get refSeq trimed length
    Ref <- names(RefSeq)
    aRefSeq <- as(seqs, "DNAStringSet")[Ref]
    aRefSeq <- stringr::str_remove_all(aRefSeq,"-")
    range <- stringr::str_locate(as.character(RefSeq), aRefSeq)
    trimedLen <- range[1, "start"] - 1

    # convert DNAset into matrix
    RefSeq <- as.matrix(RefSeq)
    mseqs <- as.matrix(seqs)
    mRef <- mseqs[Ref,]

    # deal with indels
    newmseqs <- matrix(nrow = nrow(mseqs), ncol = 0, dimnames = list(rownames(mseqs)))
    POS <- c()
    n <- 1
    for(i in seq_len(ncol(mseqs))){
        if(!"-" %in% mseqs[,i]) {
            newmseqs <- cbind(newmseqs, mseqs[,i])
            POS <- c(POS, n)
        } else {
            newmseqs[,ncol(newmseqs)] = paste0(newmseqs[,ncol(newmseqs)],mseqs[,i])
        }
        if(mRef[i] != "-") n <- n + 1
    }
    colnames(newmseqs) <- POS

    # remove redudant cols
    probe <- apply(newmseqs, 2, function(x) length(unique(x)))
    newmseqs <- newmseqs[, probe != 1]

    # remove "-"
    newmseqs[,] <- stringr::str_remove_all(newmseqs,"-")

    # update allSite information
    ALLELE <- apply(newmseqs, 2, function(x) paste(unique(x), collapse = ","))
    ALLELE <- stringr::str_replace(ALLELE, ",", "/")
    for(i in unique(ALLELE)){
        ALi <- unlist(strsplit(i, "/"))
        allS_new <- update_allS(allS_new = allS_new,
                                REF = ALi[1], ALT = ALi[2])
    }


    return(list(hap = newmseqs, allS_new = allS_new))
}

