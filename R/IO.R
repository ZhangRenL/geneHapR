# File Checked

# Checked
#' @name import_vcf
#' @title import vcf from file
#' @author Zhangrenl
#' @usage import_vcf(vcf_file = vcf_file, ...)
#' @description Read and files in the `*.vcf` structured text format,
#' as well as the compressed `*.vcf.gz` format.
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
#' @seealso
#' \code{\link[vcfR:read.vcfR]{vcfR::read.vcfR()}}
#' @return vcfR object
import_vcf <- function(vcf_file = vcf_file, ...) {
    vcf <- vcfR::read.vcfR(vcf_file, ...)
    return(vcf)
}


# Checked
#' @name import_AccINFO
#' @title imports accession information from file
#' @usage
#' import_AccINFO(file, comment.char = "#",
#'                check.names = FALSE, row.names = 1, ...)
#' @description import accession information including phenotype data,
#' accession group, location from a tab delimed table file
#' @details
#' First column should be Accessions;
#' phenos should begin from second col,
#' phenoName/group/locations should located at the first row,
#' If a dot '.' is located in pheno name, then
#' the part before the dot will be set as y axis name
#' while the followed will be set as foot of the fig.
#' @examples
#' \dontrun{
#'
#' pheno <- import_AccINFO("test.pheno.txt")
#' }
#' @importFrom utils read.delim
#' @param file file path, this file should be a tab delimed table
#' @param comment.char comment.char, start with comment.char will be ignored
#' @param ... parameters will pass to read.delim
#' @inheritParams utils::read.delim
#' @export
#' @return data.frame, Accession names were set as rownames and cols were
#' named by pheno names
import_AccINFO <- function(file, comment.char = "#",
                           check.names = FALSE, row.names = 1, ...) {
    phenos <- utils::read.delim(
        file,
        check.names = check.names,
        row.names = row.names,
        comment.char = comment.char,
        ...
    )
    return(phenos)
}


# Checked
#' @name import_gff
#' @title  import_gff
#' @description import genome annotations in `gff/gff3` format
#' @usage import_gff(gffFile, format = "GFF")
#' @examples
#' \dontrun{
#'
#'     gff <- import_gff("your.gff", format = "GFF")
#' }
#' @importFrom rtracklayer import
#' @param gffFile the gff file path
#' @param format should be one of "gff", "gff1", "gff2", "gff3", "gvf",
#' or "gtf". Default as GFF
#' @export
#' @return GRange object
import_gff <- function(gffFile, format = "GFF") {
    rtracklayer::import(gffFile, format = format)
}


# Checked
#' @name import_seqs
#' @title  import_seqs
#' @usage import_seqs(filepath, format = "fasta")
#' @examples
#' \dontrun{
#'    geneSeqs <- import_seqs(filepath = "fastaFilePath", format = "fasta")
#' }
#' @param filepath A character vector containing the path(s) to the file(s)
#' to read or write.
#' Reading files in gzip format (which usually have the '.gz' extension) is
#' supported.
#' *Note* that only DNA supported here.
#' @param format Either \code{"fasta"} (the default) or \code{"fastq"}
#' @export
import_seqs <- function(filepath, format = "fasta") {
    Biostrings::readDNAStringSet(filepath = filepath, format = format)
}


# Checked
#' @name import_MultipleAlignment
#' @title import MultipleAlignment
#' @usage import_MultipleAlignment(filepath, format = "fasta", type = "DNA")
#' @examples
#' \dontrun{
#'    geneSeqs <- import_MultipleAlignment(filepath = "fastaFilePath",
#'                                         format = "fasta",
#'                                          type = "DNA")
#'    geneSeqs <- import_MultipleAlignment(filepath = "fastaFilePath",
#'                                         format = "fasta",
#'                                          type = "Protein")
#' }
#' @param type one of 'DNA' and 'Protein'
#' @inheritParams Biostrings::readDNAMultipleAlignment
#' @export
import_MultipleAlignment <- function(filepath,
                                     format = "fasta",
                                     type = "DNA") {
    type <- toupper(type)
    if (type == "DNA") {
        Biostrings::readDNAMultipleAlignment(filepath = filepath,
                                             format = format)
    } else if (type == "PROTEIN") {
        Biostrings::readAAMultipleAlignment(filepath = filepath,
                                            format = format)
    } else if (type == "RNA") {
        Biostrings::readRNAMultipleAlignment(filepath = filepath,
                                             format = format)
    } else
        stop("type must be one of 'DNA' and 'Protein'")
}


# Checked
#' @name import_hap
#' @title  import_hap
#' @usage import_hap(file, ...)
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
#' hap <- import_hap("hapSummary_OR_hapResult_file.txt")
#' }
#' @param file hapSummary or hapResult file path
#' @param ... extras will pass to `read.delim()`
#' @export
#' @return hapSummary or hapResult
import_hap <- function(file, ...) {
    hap <- read.delim(file, header = F, ...)

    # check rows format
    if (nrow(hap) == 5)
        warning("There is only one haplotype?")
    else
        if (nrow(hap) < 5)
            stop("Please check your input file.")

    # get POS
    POS <-
        suppressWarnings(as.numeric(hap[hap[, 1] == "POS", ]))
    POS <- na.omit(POS)
    if (length(POS) == 1) {
        warning("There is only one loci?")
    } else {
        if (length(POS) < 1)
            stop("Please check your input file")
    }

    # check columns
    colnms <- c("Hap", POS)
    if (ncol(hap) - length(POS) == 2) {
        colnms <- c(colnms, "Accession")
        class(hap) <- c("hapResult", "data.frame")
    } else if (ncol(hap) - length(POS) == 3) {
        if (is.numeric(hap[, ncol(hap)])) {
            colnms <- c(colnms, "Accession", "freq")
        } else if (is.numeric(hap[, ncol(hap) - 1])) {
            colnms <- c(colnms, "freq", "Accession")
        } else {
            stop("Your file seems like 'hapSummary',
the 'freq' column (ususlly the last column) contains nonnumeric data.")
        }
        class(hap) <- c("hapSummary", "data.frame")
    } else {
        stop("Please check your input file.")
    }

    # set colnames
    colnames(hap) <- colnms

    # set attr options
    attr(hap, "options") <- c(Source = "Read from file")
    return(hap)
}


# Checked
#' @title save hap results on disk
#' @name write.hap
#' @usage write.hap(x, file = file, sep = "\t")
#' @description
#' This function will write hap result into a txt file.
#' @inherit import_hap details
#' @examples
#' \dontrun{
#'
#' write.hap(hap, file = "hap.txt")
#' }
#' @param x objec of `haplotypes` or `hapSummary` class
#' @param file file path, where to save the hap result
#' @param sep the field separator string. Values within each row of x are separated by this string.
#' @export
write.hap <- function(x, file = file, sep = "\t") {
    nc <- ncol(x)
    nm <- names(x)
    # if ("Accession" %in% nm)
    #     x[1, nm == "Accession"] <- "Accession"
    # if ("freq" %in% nm)
    #     x[1, nm == "freq"] <- "freq"
    cat("",
        file = file,
        sep = "",
        append = FALSE)
    for (i in seq_len(nrow(x))) {
        cat(
            as.matrix(x[i, ]),
            file = file,
            sep = c(rep.int(sep, nc - 1), "\n"),
            append = TRUE
        )
    }
}


# import pips
#' @importFrom magrittr `%>%`
#'
`%>%` <- magrittr::`%>%`
`%over%` <- IRanges::`%over%`
