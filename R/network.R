#' @title get_hapNet
#' @description calculate hap net using hap result
#' @usage get_hapNet(hapResult)
#' @importFrom pegas haploNet
#' @param hapResult hapResult
#' @return haplonet class
#' @export
get_hapNet <- function(hapResult){
    if(!inherits(hapResult, "haplotype"))
        hapResult <- as.haplotype(hapResult)
    hapnet=pegas::haploNet(hapResult)
}



#' @title as.haplotype
#' @description calculat Hap net from a hapSummary object
#' @importFrom ape as.DNAbin
#' @param hapResult hapResult
#' @return haplotype class
#' @export
as.haplotype <- function(hapResult){
    if(!inherits(hapResult, "hapSummary"))
        stop("'hapResult' must be  of class 'hapSummary'")
    # get freq
    freq <- hapResult$freq
    names(freq) <- hapResult$Hap
    freq <- freq[!is.na(freq)]

    # get hap mattrix
    hapDNAset <- hap2DNAstring(hapResult = hapResult)
    hapBin <- ape::as.DNAbin(hapDNAset)
    hapmatt <- unlist(as.character(as.matrix(hapBin)))

    # set as haplotype
    class(hapmatt) <- c("haplotype", "character")
    rownames(hapmatt) <- names(freq)
    N <- nrow(hapmatt)
    attr(hapmatt, "index") <- lapply(1:N, function(i) seq_len(freq[i]))
    return(hapmatt)
}


#' @title hap2DNAstring
#' @description convert hapsummary to BioStringsset class
#' @importFrom stringr str_detect str_pad str_length
#' @importFrom Biostrings DNAStringSet
#' @param hapResult hapResult
#' @export
hap2DNAstring <- function(hapResult){
    colNms <- colnames(hapResult)
    if("Accession" %in% colNms)
        hapResult <- hapResult[,colnames(hapResult) != 'Accession']
    if("freq" %in% colNms)
        hapResult <- hapResult[,colnames(hapResult) != 'freq']
    if(nrow(hapResult) <= 5) stop("There is only one Hap ?")
    meta <- hapResult[1:4,]
    hap <- hapResult[5:nrow(hapResult),]

    # padding multiallelic indel sites
    ALLELE <- meta[meta[,1] == "ALLELE",]
    multi_probe <- is.indel.allele(ALLELE)

    for(c in seq_len(ncol(hap))){
        if(multi_probe[c]){
            maxLen <- max(stringr::str_length(hap[,c]))
            hap[,c] <- stringr::str_pad(hap[,c],
                                        width = maxLen,
                                        side = "right",
                                        pad = "-")
        }
    }

    # conneting strings
    DNAseqs <- c()
    for(r in seq_len(nrow(hap))){
        DNAseqr = paste(hap[r,-1],sep = "", collapse = "")
        names(DNAseqr) <- hap[r,1]
        DNAseqs <- c(DNAseqs, DNAseqr)
    }
    hapDNASet <- Biostrings::DNAStringSet(DNAseqs, use.names = T, start = 1)
    return(hapDNASet)
}



is.indel.allele <- function(allele.vector){
    probe1 <- stringr::str_detect(allele.vector,"/")
    probe <- lapply(stringr::str_split(allele.vector,"[,/]"),
                    stringr::str_length)
    probe1 & unlist(lapply(probe, function(i) max(i)>1))
}

is.biallelic.allele <- function(allele.vector){
    probe1 <- stringr::str_detect(allele.vector,"/")
    !stringr::str_detect(allele.vector,",") & probe1
}

is.multiallelic.allele <- function(allele.vector){
    probe1 <- stringr::str_detect(allele.vector,"/")
    stringr::str_detect(allele.vector,",") & probe1
}
