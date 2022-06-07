#' @name seqs2hap
#' @title generate hap results from seqs
#' @description
#' @usage
#' @examples
#' @importFrom Biostrings width
#' @importFrom muscle muscle
#' @param seqs DNAStringSet or DNAMultipleAlignment
#' @param needAlign
#' @param needTrim
#' @param ... parameters will pass to musle::musle() for DNA multi alignment
seqs2hap <- function(seqs, needAlign = TRUE, needTrim = TRUE, ...){
    if(needAlign) seqs <- allignSeqs(seqs) else
        cat("Multialignment will not perform due to 'needAlign' is FALSE\n")

    if(needTrim){
        seqs <- trimSeqs(seq)
    } else {
        cat("Seqs will not be trimed due to 'needTrim' is FALSE\n")
        widths <- unique(Biostrings::width(seq))
        if(length(widths) != 1)
            stop("please make sure all seqs with same length")
    }

}

allignSeqs <- function(seqs, ...){
    muscle::muscle(seq)
    dnaseq <- Biostrings::maskedwidth(aligned)
}

trimSeqs <- function(seqs){
    as.
}
a = maskGaps(aligned,min.fraction=0.02, min.block.width=1)
detail(a)
Biostrings::writePairwiseAlignments(a)
class(a)
