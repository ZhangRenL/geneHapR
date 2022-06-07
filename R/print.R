#' @exportS3Method print haptypes
print.haptypes <- function(x){
    freq <- attr(x, "freq")
    accAll <- attr(x, "AccAll")
    accRemain <- attr(x, "AccRemain")
    accRemoved <- attr(x,"AccRemoved")

    cat("Accessions:")
    cat("\n All:\t ", length(accAll))
    cat("\n Removed:", length(accRemoved)," ")
    if(length(accRemoved) <= 10){
        cat(accRemoved, sep = ", ")
    } else {
        cat(head(accRemoved), "...", sep = ", ")
    }
    cat("\n Remain: ", length(accRemain))
    cat("\n\nhap frequences:")
    print(freq)

    options <- attr(x, "options")
    cat("\nOptions:")
    for(i in seq_len(length(options))){
        nn <- paste0(names(options)[i],":")
        cat("\n ",stringi::stri_pad_right(nn, 11)," ", sep = "")
        cat(options[i])
    }
}


#' @exportS3Method print hapSummary
print.hapSummary <- function(x){
    freq <- table(names(attr(x, "hap2acc")))
    ALLELE <- x[x$Hap %in% "ALLELE",
                        !names(x) %in% c("Hap","Accession", "freq")]
    nIndel <- table(is.indel.allele(ALLELE))

    cat("\nAccssions:", sum(freq))
    cat("\nSites: ", sum(nIndel),
        "\n  Indels: ", nIndel["TRUE"],
        "\n  SNPs:   ", nIndel["FALSE"], sep = "")
    cat("\n\nHaplotypes:", length(freq), "\n  ")
    s <- cbind(x$Hap,
               x$freq,
               sapply(x$Accession,
                      function(x) {
                          l = unlist(strsplit(x,";"))
                          if(length(l) > 6){
                              paste(c(head(l, 6), "..."),
                                    collapse = ", ", sep = "")
                          } else {
                              paste(l, collapse = ", ")
                          }
                      })
    )
    for(i in seq_len(nrow(s))) if(!NA %in% s[i,]) cat(s[i,], "\n  ", sep = "\t")
    options <- attr(x, "options")
    cat("\nOptions:")
    for(i in seq_len(length(options))){
        cat("\n ",names(options)[i],":\t", sep = "")
        cat(options[i])
    }
}


#' @exportMethod write
write.hap <- function(x, file = file, ...){
    nc <- ncol(x)
    nm <- names(x)
    if("Accession" %in% nm) x[1,nm == "Accession"] <- "Accession"
    if("freq" %in% nm) x[1,nm == "freq"] <- "freq"
    cat("",file = file, sep = "", append = FALSE)
    for(i in seq_len(nrow(x))){
        cat(as.matrix(x[i,]), file = "a.txt", sep = c(rep.int("\t", nc-1), "\n"),append = TRUE)
    }
}


#' @exportS3Method write.hap haptypes
write.hap.haptypes <- function(x, file = file, ...){
    write.hap(x, file = file, ...)
}


#' @exportS3Method write.hap hapSummary
write.hap.hapSummary <- function(x, file = file, ...){
    write.hap(x, file = file, ...)
}
