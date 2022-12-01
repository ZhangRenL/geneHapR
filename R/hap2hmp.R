#' @name hap2hmp
#' @title Convert hapResult object to hapmap (hmp) format,
#' for interact with other packages
#' @param hap object of "hapResult" class
#' @param hmp object of "data.frame" class in hapmap format
#' @param ... Parameters not used.
#' @inheritDotParams table2hap
#' @return a data.frame in hapmap format.
#' @usage
#' hap2hmp(hap)
#' @examples
#' \donttest{
#' data("geneHapR_test")
#' hmp <- hap2hmp(hapResult)
#' hap <- hmp2hap(hmp)
#' }
#' @export
hap2hmp <- function(hap){
    hap <- t(hap)
    hap <- hap[-1,]
    data.frame(hap[1:10,1:10])[1,]
    alleles <- hap[, 4]
    chrom <- hap[, 1]
    pos <- hap[, 2]
    rs <- paste(chrom, pos, sep = "_")
    colnames(hap) <- hap[nrow(hap),]
    hap <- hap[, -c(1:4)]
    hmp <- data.frame("rs#" = rs,
                      alleles = alleles,
                      chrom = chrom,
                      pos = pos,
                      strand = NA,
                      "assembly#" = NA,
                      center = NA,
                      protLSID = NA,
                      assayLSID = NA,
                      panelLSID = NA,
                      QCcode = NA, check.names = FALSE)
    hmp <- cbind(hmp, hap)
    hmp <- hmp[-nrow(hmp),]
    return(hmp)
}


#' @name hap2hmp
#' @inheritParams table2hap
#' @export
hmp2hap <- function(hmp,
                    hapPrefix = "H",
                    hetero_remove = TRUE,
                    na_drop = TRUE, ...){
    if(ncol(hmp) < 10) stop("hapmap file column numer should more than 10")
    data <- hmp[, c(3, 4, 2, 2, 2, 11 : ncol(hmp))]
    data[, 2] <- as.numeric(data[, 2])
    table2hap(data,
              hapPrefix = hapPrefix,
              hetero_remove = hetero_remove,
              na_drop = na_drop)
}
