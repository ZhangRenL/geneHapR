#' @name import_vcf
#' @title import vcf from file
#' @author Zhangrenl
#' @usage import_vcf(vcf_file = vcf_file, ...)
#' @description  import vcf file
#' @examples
#' \dontrun{
#' data("quickHap_test")
#' vcfR::write.vcf(vcf, file = "test.vcf.gz")
#' vcf <- import_vcf(file = "test.vcf.gz")
#' unlink("test.vcf.gz")
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
#' data("quickHap_test")
#' write.table(pheno, file = "test.pheno.txt", quote = FALSE, sep = "\t")
#' pheno <- import_pheno("test.pheno.txt")
#' unlink("test.pheno.txt")
#'}
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



#### Import the pipe operator from magrittr ####
#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @importFrom IRanges %over%
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

