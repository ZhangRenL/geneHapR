% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf2hap.R
\name{vcf2hap}
\alias{vcf2hap}
\title{generat haps from vcf}
\usage{
vcf2hap(vcf,
        hapPrefix = "H",
        filter_Chr = FALSE,
        Chr = Chr,
        filter_POS = FALSE,
        startPOS = startPOS,
        endPOS = endPOS,
        hyb_remove = TRUE,
        na.drop = TRUE, file = file)
}
\arguments{
\item{vcf}{vcfR object imported by \code{import_vcf}}

\item{hapPrefix}{Prefix of hap names, default as "H"}

\item{filter_Chr}{logical, whether filter vcf by chromosome or not, default
as \code{FALSE}. If set as \code{TRUE}, \code{Chr} is needed}

\item{Chr}{Chromosome name, needed when \code{filter_Chr} was set as \code{TRUE}}

\item{filter_POS}{logical, whether filter vcf by position or not. Default
as \code{FALSE}. If set as \code{TRUE}, \code{startPOS} and \code{endPOS} are needed}

\item{startPOS, endPOS}{Start and end position, needed when \code{filter_POS} was
set as \code{TRUE}. In addition, \code{startPOS} must less than \code{endPOS}}

\item{hyb_remove}{whether remove accessions contains hybrid site or not.
Default as \code{TRUE}}

\item{na.drop}{whether drop accessions contains unknown allele site or not
Default as \code{TRUE}.}

\item{file}{file path where to save the hap result. If missing, nothing will
be saved to disk}
}
\value{
\code{haptypes} object, first four rows are meta information:
\code{CHR}, \code{POS}, \code{INFO}, \code{ALLELE.}
\code{Hap} names were located in the first column, \code{Accessions} were placed at
the last column.
}
\description{
Generate hap result from \code{vcfR} object.
A simple filter by position was provided in this function,
however it's prefer to filter vcf (\code{vcfR} object) through
\code{\link[geneHapR:filter_vcf]{filter_vcf()}}.
}
\details{
The class of returned object is \code{data.frame} class with some
additional \code{attributes} generated during processing,
such as \code{options}: options used during generated haplotypes,
\code{AccRemoved}: names of removed Accessions due to \code{NA} or hybrid sites,
\code{AccRemain}: names of remained Accessions in haplotype results,
\code{AccAll}: names of all Accessions in souce vcf file,
\code{freq}: frequency of each haplotype.
Though, only summary view of those \code{attributes} will be displayed, the full
information could be displayed with \code{attr(x, "options")} or \code{attributes(x)}.
}
\note{
inhet test
description  generate hap format from vcf.
}
\examples{
\dontrun{

#### 1. generate hap result from vcf file ####
# imoprt vcf data
vcf <- import_vcf("myVcfFile.vcf.gz")

# calculate hap result with vcfR
hap <- vcf2hap(vcf, hyb_remove = TRUE, na.drop =TRUE)


#### 2. calculate hap from DNA seqs ####
# import DNA seqs from fasta file
seqs <- import_seqs("myseq.fa", format = "fasta")

# aline multi-seqs
aliedSeqs <- allignSeqs(seqs)

# trim alined seqs
trimedSeqs <- trimSeqs(seqs,
                       minFlankFraction = 0.1)

# calculate hap result with trimed seqs
hap <- seqs2hap(trimedSeqs,
                Ref = names(seqs)[1],
                hyb_remove = TRUE,
                maxGapsPerSeq = 0.25,
                hapPrefix = "H")



#### 3. summary and visualizatiion ####
# summary hap result for plot hap table
hapResult <- hap_result(hap)

# plot hap table
plotHapTable(hapResult)

# user could change/add themes, title, ...
p <- plotHapTable(hapResult)
p + ggplot2::theme_classic()
p + ggplot2::labs(title = "new title")
}
}
\seealso{
extract genotype from vcf:
\code{\link[vcfR:extract_gt_tidy]{vcfR::extract_gt_tidy()}},
import vcf files:
\code{\link[geneHapR:import_vcf]{import_vcf()}} (preferred) and
\code{\link[vcfR:read.vcfR]{vcfR::read.vcfR()}},
filter vcf according \strong{position} and \strong{annotations}:
\code{\link[geneHapR:filter_vcf]{filter_vcf()}}
}
\author{
Zhangrenl
}