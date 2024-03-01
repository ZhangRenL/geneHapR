#' @name import_vcf
#' @title Import VCF from File
#' @author Zhangrenl
#' @usage import_vcf(file = file, ...)
#' @description import *.vcf structured text format,
#' as well as the compressed `*.vcf.gz` format.
#' @examples
#'
#' vcfPath <- system.file("extdata", "var.vcf.gz", package = "geneHapR")
#' vcf <- import_vcf(file = vcfPath)
#' vcf
#'
#' @importFrom vcfR read.vcfR
#' @param file file path of VCF file
#' @param ... pass to `vcfR::read.vcfR()`
#' @usage import_vcf(file = file, ...)
#' @export
#' @seealso
#' \code{\link[vcfR:read.vcfR]{vcfR::read.vcfR()}}
#' @return vcfR object
import_vcf <- function(file = file, ...) {
    vcf <- vcfR::read.vcfR(file, ...)
    return(vcf)
}


#' @title import_plink.pedmap
#' @name import_plink.pedmap
#' @description used for import regular p.link file stored in map and ped format
#' @param root The file name without suffix. This function only support p.link file
#'   format stored in "map" and "ped" format, the file names after removed suffix
#'   should be same with each other.
#' @param sep_ped a character indicate the separation of ped file
#' @param sep_map a character indicate the separation of map file
#' @param pedfile,mapfile if `root` is missing then `pedfile` and `mapfile` are needed
#' @return list, contains map information stored in data.frame and ped
#'   information stored in data.frame
#' @importFrom utils read.table
#' @usage
#'   import_plink.pedmap(root = root,
#'                       sep_ped = "\t", sep_map = "\t",
#'                       pedfile = pedfile, mapfile = mapfile)
#' @inherit plink.pedmap2hap examples
#' @export
import_plink.pedmap <- function(root = root,
                                sep_ped = "\t", sep_map = "\t",
                                pedfile = pedfile, mapfile = mapfile){
    if(!missing(root)){
        pedfile <- paste0(root, ".ped")
        mapfile <- paste0(root, ".map")
    }
    ped <- read.table(file <- pedfile, sep = sep_ped, # nrows = 100,
                      header = FALSE)
    map <- read.table(file = mapfile, header = FALSE, sep = sep_map)
    if((ncol(ped) - 6) / nrow(map) != 2)
        warning("ped or map file may corrupted")
    return(list(map = map, ped = ped))
}



#' @name import_AccINFO
#' @title Import Accession Information from File
#' @usage
#' import_AccINFO(file, comment.char = "#",
#'                check.names = FALSE, row.names = 1, ...)
#' @description import accession information including phenotype data,
#' accession group, location from a tab delimed table file
#' @details
#' First column should be Accessions;
#' phenos/accession information should begin from second column,
#' phenoName/group/locations should located at the first row,
#' If a dot '.' is located in pheno name, then
#' the part before the dot will be set as y axis name
#' and the latter will be set as foot when plot figures.
#' @examples
#'
#' oldDir <- getwd()
#' temp_dir <- tempdir()
#' if(! dir.exists(temp_dir))
#'   dir.create(temp_dir)
#' setwd(temp_dir)
#' data("geneHapR_test")
#' write.table(pheno, file = "test.pheno.txt", sep = "\t")
#' pheno <- import_AccINFO("test.pheno.txt")
#' pheno
#' setwd(oldDir)
#'
#' @importFrom utils read.delim
#' @param file file path, this file should be a tab delimed table
#' @inheritParams utils::read.delim
#' @export
#' @return data.frame, Accession names were set as rownames and columns were
#' named by pheno/info names
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



#' @name import_gff
#' @title Import Annotations from GFF Format File
#' @description import genome annotations in GFF/GFF3 format
#' @usage import_gff(gffFile, format = "gff")
#' @examples
#'
#' gff.Path <- system.file("extdata", "annotation.gff", package = "geneHapR")
#' gff <- import_gff(gff.Path, format = "gff")
#' gff
#'
#' @importFrom rtracklayer import
#' @param gffFile the gff file path
#' @param format should be one of "gff", "gff1", "gff2", "gff3", "gvf",
#' or "gtf". Default as gff
#' @export
#' @return GRange object
import_gff <- function(gffFile, format = "gff") {
    gff <- rtracklayer::import(gffFile, format = format)
    gff$Parent <- lapply(gff$Parent,
                          function(x) if(identical(x, character(0))) "" else x)
    gff$Parent <- unlist(gff$Parent)
    return(gff)
}


#' @name import_bed
#' @title import annotation files in BED format
#' @description import bed files contains annotations into R as GRanges object
#' @details
#' If there is no genome annotation file in GFF format for your interest species,
#' a BED file is convenient to custom a simple annotation file for a gene.
#' Here we suggest two type of BED format: BED6 and BED4.
#'
#' As the definition of
#' [UCSC](http://genome.ucsc.edu/FAQ/FAQformat.html#format1).
#' The BED6 contains 6 columns, which are 1) chrom, 2) chromStart, 3) chromEnd,
#' 4) name, 5) score and 6) strand. The BED4 format contains the first 4 column
#' of BED6 format.
#'
#' However, in gene haplotype statistics, we only care about the type of each site.
#' Thus we use the fourth column to definition the
#' transcripts name and "CDS" or "URTs", separated by a space, eg.:
#'
#'   \code{Chr8   678   890   HD1.1 CDS   .   -}
#'
#'   \code{Chr8   891   989   HD1.1 five_prime_UTR   .   -}
#'
#'   \code{Chr8   668   759   HD1.2 CDS   .   -}
#'
#'   \code{Chr8   908   989   HD1.2 CDS   .   -}
#'
#' This example indicate a small gene named as HD1 have two transcripts,
#' named as HD1.1 and HD1.2, separately. HD1 has a CDS and a UTR region;
#' while HD1.2 has two CDS region.
#'
#' @inheritParams rtracklayer::import.bed
#' @param quite whether show message
#' @importFrom rtracklayer import.bed
#' @importFrom stringr str_count
#' @usage
#' import_bed(con, quite = FALSE)
#' @examples
#'
#' bed.Path <- system.file("extdata", "annotation.bed6", package = "geneHapR")
#' bed <- import_bed(bed.Path)
#' bed
#'
#' @return GRange object
#' @export
import_bed <- function(con, quite = FALSE){
    bed <- rtracklayer::import.bed(con)
    if(length(bed)  == 0) stop("the bed file seems has no annotations")
    probe <- stringr::str_count(bed$name, " ") %>% unique()

    m <- "the fourth column should contains transcript name and region type, seprated by a space"
    if(length(probe) != 1) {
        stop(m)
    } else {
        if(probe != 1) stop(m)
    }

    if(quite){
        cat("Please note that:")
        cat("  this function was only intend to import customed annotations\n")
        cat("  in BED4/BED6 format for haplotype statistics")
    }

    name <- sapply(bed$name, function(x) strsplit(x, " ")[[1]][1]) %>% unlist()
	type <- sapply(bed$name, function(x) strsplit(x, " ")[[1]][2]) %>% unlist()
    bed$Parent <- bed$Name <- name
    bed$type <- type
    gff <- bed
    gff$Parent <- lapply(gff$Parent,
                         function(x) if(identical(x, character(0))) "" else x)
    gff$Parent <- unlist(gff$Parent)
    return(gff)
}


#' @name import_seqs
#' @title Import Sequences
#' @description import DNA sequences in FASTA format
#' @usage import_seqs(filepath, format = "fasta")
#' @examples
#'
#' seqPath <- system.file("extdata", "seqs.fa", package = "geneHapR")
#' geneSeqs <- import_seqs(filepath = seqPath, format = "fasta")
#'
#' @param filepath A character vector containing the path to the DNA sequences file.
#' Reading files in gzip format (which usually have the '.gz' extension) is
#' supported.
#' *Note* that only DNA supported here.
#' @param format Either "fasta" (the default) or "fastq"
#' @return object of DNAStringSet class
#' @export
import_seqs <- function(filepath, format = "fasta") {
    Biostrings::readDNAStringSet(filepath = filepath, format = format)
}



#' @name import_MultipleAlignment
#' @title Import MultipleAlignment Result
#' @description import sequences algned results
#' @usage import_MultipleAlignment(filepath, format = "fasta", type = "DNA")
#' @examples
#' aliSeqPath <- system.file("extdata", "seqs.fa", package = "geneHapR")
#'
#' geneSeqs <- import_MultipleAlignment(filepath = aliSeqPath,
#'                                      format = "fasta",
#'                                      type = "DNA")
#' geneSeqs <- import_MultipleAlignment(filepath = aliSeqPath,
#'                                      format = "fasta",
#'                                      type = "Protein")
#'
#' @param type one of "DNA" and "Protein"
#' @inheritParams Biostrings::readDNAMultipleAlignment
#' @return object of DNAMultipleAlignment
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



#' @name import_hap
#' @title Import hapResult/hapSummary
#' @description
#' This function could be used for import hap result or hap summary result.
#' The type of returned object is decided by input file, see details.
#' @details
#' The hap result and hap summary result have common features.
#'   The common features of these two types are:
#'     First four rows contains extra information: CHR, POS, INFO and ALLELE
#'     Hap names were in the first column.
#'   The differences are:
#'     Hap summary result have a freq column while hap result not.
#'     Rows represent haplotypes in hap summary result, while rows represent accessions in hap result.
#'     In addtion, the accessions of each haplotype in hap summary result were separated by ";".
#' @examples
#'
#' oldDir <- getwd()
#' temp_dir <- tempdir()
#' if(! dir.exists(temp_dir))
#'   dir.create(temp_dir)
#' setwd(temp_dir)
#' data("geneHapR_test")
#' write.hap(hapResult, file = "test.pheno.txt", sep = "\t")
#' hap <- import_hap("test.pheno.txt")
#' hap
#' setwd(oldDir)
#'
#' @param file hapSummary or hapResult file path.
#' @param type file type, if not auto, should be one of hapSummary or hapResult
#' @param type the content type of imported file, should be one of c("hapResult", "hapSummary")
#' @param ... extras will pass to `read.delim()`
#' @export
#' @return hapSummary or hapResult
import_hap <- function(file, type = "auto", ...) {
    hap <- read.delim(file, header = F, ...)

    # check rows format
    if (nrow(hap) == 5)
        warning("There is only one haplotype?")
    else
        if (nrow(hap) < 5)
            stop("Please check your input file.")

    # set data class
    if (tolower(type) == "auto"){
        cfreq <- where_('freq', hap, 'c')
        cAccession <- where_('Accession', hap, 'c')
        if(cfreq == -1){
            if(cAccession == -1){
                stop("Accession and freq column wasn't found")
            } else {
                class(hap) <- c("hapResult", "data.frame")
                hap <- hap[,seq_len(cAccession)]
            }
        } else {
            if(cAccession == -1){
                stop("Accession column wasn't found")
            } else {
                class(hap) <- c("hapSummary", "data.frame")
            }
        }
    } else if (tolower(type == "hapsummary")){
        class(hap) <- c("hapSummary", "data.frame")
    } else if (tolower(type == "hapresult")){
        class(hap) <- c("hapResult", "data.frame")
    } else stop("type should be one of 'auto', 'hapSummary', 'hapResult'")

     # get POS
    POS <-
        suppressWarnings(as.numeric(hap[hap[, 1] == "POS", ]))
    if (length(na.omit(POS)) == 1) {
        warning("There is only one loci?")
    } else {
        if (length(na.omit(POS)) < 1)
            stop("Please check your input file")
    }

    # if(is.na(POS[length(POS)]) & is.na(POS[length(POS) - 1])){
    #     # the last two vectors were not numeric
    #     # hapSummary exported by old version
    #     class(hap) <- c("hapSummary", "data.frame")
    # } else if(is.na(POS[length(POS)]) & !is.na(POS[length(POS) - 1])){
    #     # the last one of POS is NA and the next to last is numeric
    #     # hapResult
    #     class(hap) <- c("hapResult", "data.frame")
    # } else if(!is.na(POS[length(POS)]) & is.na(POS[length(POS) - 1])){
    #     # the last one of POS is numeric and the next to last is na
    #     # hapSummary exported by new version
    #     class(hap) <- c("hapSummary", "data.frame")
    # } else {
    #     stop("")
    # }

    # set column names
    POS[1] <- "Hap"
    if(inherits(hap, "hapSummary")) {
        POS[(length(POS) - 1) : length(POS)] <- c("Accession", "freq")
        hap[1:4, ncol(hap)] <- NA
        hap[1:4, ncol(hap) - 1] <- ""
    }
    if(inherits(hap, "hapResult")) {
        POS[length(POS)] <- c("Accession")
        hap[1:4, ncol(hap)] <- ""
    }
    # set colnames
    colnames(hap) <- POS


    # set attr options
    if(inherits(hap, "hapResult")){
        accAll <- hap$Accession
        attr(hap, "AccAll") <- accAll
        hap2acc <- hap$Accession[-c(1:4)]
        names(hap2acc) <- hap$Hap[-c(1:4)]
        attr(hap, "hap2acc") <- hap2acc
        attr(hap, "freq") <- table(hap$Hap[-c(1:4)])

    } else {
        accAll <- c()
        for(i in seq_len(nrow(hap))){
            if(i <= 4) next
            accAlli <- strsplit(hap$Accession[i], ";") %>% unlist()
            names(accAlli) <- rep(hap$Hap[i], length(accAlli))
            accAll <- c(accAll, accAlli)
        }
        attr(hap, "hap2acc") <- accAll
        names(accAll) <- NULL
        attr(hap, "AccRemain") <- attr(hap, "AccAll") <- accAll
        attr(hap, "freq") <- hap$freq[-c(1:4)]
    }
    attr(hap, "options") <- c(Source = "Read from file")
    return(hap)
}



#' @title Save Haplotype Results to Disk
#' @name write.hap
#' @description
#' This function will write hap result into a txt file.
#' @inherit import_hap details
#' @examples
#'
#'
#' oriDir <- getwd()
#'  temp_dir <- tempdir()
#'  if(! dir.exists(temp_dir))
#'    dir.create(temp_dir)
#'  setwd(temp_dir)
#' data("geneHapR_test")
#' write.hap(hapResult, file = "hapResult.txt")
#' setwd(oriDir)
#'
#' @param x objec of hapResult or hapSummary class
#' @param file file path, where to save the hap result/summary
#' @param sep the field separator string. Values within each row of x are separated by this string.
#' Default as "`\t`"
#' @param pheno,phenoName the phenotype data.frame, only used for export hapResult object.
#' @return No return value
#' @export write.hap
write.hap <- function(x, file = file, sep = "\t", pheno = pheno, phenoName = phenoName) {
    # add Summaries
    hap2acc <- attr(x, "hap2acc")
    haps <- names(hap2acc) %>% unique() %>% length()
    if(inherits(x, "hapResult")){
        x$Accession[1] <- paste0("Haplotypes: ", haps)
        x$Accession[2] <- paste0("Individuals: ", length(hap2acc))
        x$Accession[3] <- paste0("Variants: ", ncol(x) - 2)
        x$Accession[4] <- "Accession"

        # add phenotype
        if(! missing(pheno)){
            if("INFO" %in% row.names(pheno))
                stop("There is a sample named as 'INFO', which is a reserved word")
            if("CHR" %in% row.names(pheno))
                stop("There is a sample named as 'CHR', which is a reserved word")
            if("POS" %in% row.names(pheno))
                stop("There is a sample named as 'POS', which is a reserved word")
            if("ALLELE" %in% row.names(pheno))
                stop("There is a sample named as 'ALLELE', which is a reserved word")
            if(! missing(phenoName)){
                if(! phenoName %in% names(pheno))
                    stop("'", phenoName, "' was not found in `pheno`")
                x <- cbind(x, pheno[row.names(x), phenoName])
            }
            x <- cbind(x, pheno[row.names(x), ])
        }
    }

    if(inherits(x, "hapSummary")){
        x$Accession[1] <- "Haplotypes: "
        x$freq[1] <- haps
        x$Accession[2] <- "Individuals: "
        x$freq[2] <- length(hap2acc)
        x$Accession[3] <- "Variants: "
        x$freq[3] <- ncol(x) - 3
        x$Accession[4] <- "Accession"
        x$freq[4] <- "freq"
    }
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

# find the position of `x` from m,
# t: r, return the row number
# t: c, return the column number
# if not found, return -1
where_ <- function(x, m, t = 'r'){
    m <- as.matrix(m)
    p <- x == m
    p[is.na(p)] <- FALSE
    if(any(p)){
        switch (t,
                'r' = which(p) %% nrow(m),
                'c' = which(p) %/% nrow(m) + 1
        )
    } else {
        -1
    }
}
#' @importFrom tidyr `%>%` unite pivot_wider matches chop
#' @importFrom IRanges `%over%`
`%>%` <- tidyr::`%>%`
`%over%` <- IRanges::`%over%`
