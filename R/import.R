#' @name import_vcf
#' @title import vcf from file
#' @author Zhangrenl
#' @usage import_vcf(vcf_file = vcf_file, ...)
#' @description  import vcf file
#' @examples
#' \dontrun{
#'
#' vcf <- import_vcf(file = "test.vcf.gz")
#' }
#' @importFrom vcfR read.vcfR
#' @param vcf_file file path of vcf
#' @param ... pass to vcfR::read.vcfR
#' @usage import_vcf(vcf_file = vcf_file, ...)
#' @export
#' @return vcfR object
import_vcf <- function(vcf_file = vcf_file, ...) {
    vcf <- vcfR::read.vcfR(vcf_file, ...)
    return(vcf)
}


#' @name import_pheno
#' @title imports phenos from file
#' @usage import_pheno(phenoFile, comment.char = "#", ...)
#' @description first col should be Accessions
#' phenos should loacted begin second col，
#' phenoName should at the first row,
#' If "." located in pheno name, the former part of phenoName will be set as y axis name
#' and the latter part will be set as foot of the fig.
#' @examples
#' \dontrun{
#'
#' pheno <- import_pheno("test.pheno.txt")
#' }
#' @importFrom utils read.delim
#' @param phenoFile pheno file path, should be a table separated by tab
#' @param comment.char comment.char, start with comment.char will be ignored
#' @param ... parameters will pass to read.delim
#' @export
#' @return data.frame, Accession names were set as rownames and cols were
#' named by pheno names
import_pheno <- function(phenoFile, comment.char = "#", ...){
    phenos <- utils::read.delim(phenoFile,
                                check.names = FALSE,
                                row.names = 1,
                                comment.char = comment.char, ...)
    return(phenos)
}


#' @name import_gff
#' @title  import_gff
#' @usage import_gff(gffFile, format = "GFF")
#' @examples
#' \dontrun{
#'
#'     gff <- import_gff("your.gff", format = "GFF")
#' }
#' @importFrom rtracklayer import
#' @param gffFile gff file path
#' @param format defalt as GFF
#' @export
#' @return GRange object
import_gff <- function(gffFile, format = "GFF"){
    rtracklayer::import(gffFile, format = "GFF")
}

#' @name import_seqs
#' @title  import_seqs
#' @usage import_seqs(file, format = "fasta", ...)
#' @examples
#' \dontrun{
#'    geneSeqs <- import_seqs(file = "fastaFilePath", format = "fasta")
#' }
#' @importFrom rtracklayer import
#' @param file gff file path
#' @param format Either "fasta" (the default) or "fastq"
#' @param ... Others parameters supported by Biostrings::readDNAStringSet
#' @export
#' @return DNAstringSet
import_seqs <- function(file, format = "fasta", ...){
    Biostrings::readDNAStringSet(filepath = file, format = format, ...)
}

#' @name import_hapResult
#' @title  import_hapResult
#' @usage import_hapResult(file, ...)
#' @description
#' This function could be used for import hap result or hap summary result.
#' The type of returned object is decided by hap result format, see details.
#' @details
#' The hap result and hap summary result have common features.
#'   The common features of these two types are:
#'     First four rows contains extra information: CHR, POS, INFO and ALLELE
#'     Hap names were in the first column.
#'   The differences are:
#'     Hap summary result have a freq column while hap result not.
#'     Rows represent haplotypes in hap summary result, while rows represent accessions in hap result.
#'     In addtion, the accessions of each haplotype in hap summary result were separated by ';'.
#' @examples
#' \dontrun{
#'
#' hapResult <- import_hapResult("your_hapResult_file.txt")
#' }
#' @param file hapResult file path
#' @param ... extras will pass to `read.delim()`
#' @export
#' @return hapSummary or haptypes
import_hapResult <- function(file, ...){
    hapResult <- read.delim(file, header = F, ...)

    # check rows format
    if(nrow(hapResult) == 5)
        warning("There is only one haplotype?") else
            if(nrow(hapResult) < 5)
                stop("Please check your input file.")

    # get POS
    POS <- suppressWarnings(as.numeric(hapResult[hapResult[,1] == "POS",]))
    POS <- na.omit(POS)
    if(length(POS) == 1)
        warning("There is only one loci?") else
            if(length(POS) < 1)
                stop("Please check your input file")

    # check columns
    cn <- c("Hap", POS)
    if(ncol(hapResult) - length(POS) == 2){
        cn <- c(cn, "Accession")
        class(hapResult) <- c("haptypes", "data.frame")
    } else if(ncol(hapResult) - length(POS) == 3){
        if(is.numeric(hapResult[,ncol(hapResult)])){
            cn <- c(cn, "Accession", "freq")
        } else if(is.numeric(hapResult[,ncol(hapResult) - 1])){
            cn <- c(cn, "freq", "Accession")
        } else {
            stop("Can't find Please check your input file.")
        }
        class(hapResult) <- c("hapSummary", "data.frame")
    } else {
        stop("Please check your input file.")
    }

    # set colnames
    colnames(hapResult) <- cn

    # set attr options
    attr(hapResult, "options") <- c(Source = "Read from file")
    return(hapResult)
}


#' @title save hap results on disk
#' @name write.hap
#' @usage write.hap(x, file = file, sep = "\t")
#' @description
#' This function will write hap result into a txt file.
#' @inherit import_hapResult details
#' @examples
#' \dontrun{
#'
#' write.hap(hap, file = "hap.txt")
#' }
#' @param x objec of haplotypes or hapSummary
#' @param file file path, where to save the hap result
#' @param sep the field separator string. Values within each row of x are separated by this string.
#' @export
write.hap <- function(x, file = file, sep = "\t"){
    nc <- ncol(x)
    nm <- names(x)
    if("Accession" %in% nm) x[1,nm == "Accession"] <- "Accession"
    if("freq" %in% nm) x[1,nm == "freq"] <- "freq"
    cat("",file = file, sep = "", append = FALSE)
    for(i in seq_len(nrow(x))){
        cat(as.matrix(x[i,]), file = file, sep = c(rep.int(sep, nc-1), "\n"),append = TRUE)
    }
}


# import pips
`%>%` <- magrittr::`%>%`
`%over%` <- IRanges::`%over%`

