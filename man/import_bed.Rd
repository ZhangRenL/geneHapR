% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IO.R
\name{import_bed}
\alias{import_bed}
\title{import annotation files in BED format}
\usage{
import_bed(con, quite = FALSE)
}
\arguments{
\item{con}{A path, URL, connection or \code{BEDFile} object. For the
    functions ending in \code{.bed}, \code{.bedGraph} and \code{.bed15},
    the file format is indicated by the function name. For the base
    \code{export} and \code{import} functions, the format must be
    indicated another way. If \code{con} is a path, URL or connection,
    either the file extension or the \code{format} argument needs to be
    one of \dQuote{bed}, \dQuote{bed15}, \dQuote{bedGraph},
    \dQuote{bedpe}, \dQuote{narrowPeak}, or
    \dQuote{broadPeak}. Compressed
    files (\dQuote{gz}, \dQuote{bz2} and
    \dQuote{xz}) are handled transparently.
  }

\item{quite}{whether show message}
}
\value{
GRange object
}
\description{
import bed files contains annotations into R as GRanges object
}
\details{
If there is no genome annotation file in GFF format for your interest species,
a BED file is convenient to custom a simple annotation file for a gene.
Here we suggest two type of BED format: BED6 and BED4.

As the definition of
\href{http://genome.ucsc.edu/FAQ/FAQformat.html#format1}{UCSC}.
The BED6 contains 6 columns, which are 1) chrom, 2) chromStart, 3) chromEnd,
4) name, 5) score and 6) strand. The BED4 format contains the first 4 column
of BED6 format.

However, in gene haplotype statistics, we only care about the type of each site.
Thus we use the fourth column to definition the
transcripts name and "CDS" or "URTs", separated by a space, eg.:

\code{Chr8   678   890   HD1.1 CDS   .   -}

\code{Chr8   891   989   HD1.1 five_prime_UTR   .   -}

\code{Chr8   668   759   HD1.2 CDS   .   -}

\code{Chr8   908   989   HD1.2 CDS   .   -}

This example indicate a small gene named as HD1 have two transcripts,
named as HD1.1 and HD1.2, separately. HD1 has a CDS and a UTR region;
while HD1.2 has two CDS region.
}
\examples{

bed.Path <- system.file("extdata", "annotation.bed6", package = "geneHapR")
bed <- import_bed(bed.Path)
bed

}
