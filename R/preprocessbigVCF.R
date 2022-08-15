# File Checked
#' @name filterLargeVCF
#' @title Pre-process of Large VCF File(s)
#' @description
#' Filter/extract one or multiple gene(s)/range(s) from a large
#' `*.vcf/*.vcf.gz` file.
#' @usage
#'  filterLargeVCF(VCFin = VCFin, VCFout = VCFout,
#'                 Chr = Chr,
#'                 POS = POS,
#'                 start = start,
#'                 end = end,
#'                 override = TRUE)
#' @param VCFin Path of input `*.vcf/*.vcf.gz` file.
#' @param VCFout Path(s) of output `*.vcf/*.vcf.gz` file.
#' @param Chr a single CHROM name or CHROM names vector.
#' @param POS,start,end provide the range should be extract from orignal vcf.
#'  `POS`: a vector consist with start and end position or a list
#'  with length equal to `Chr`, eg.: `list(c(1,200), c(300,500), c(300,400))`
#'  indicates 3 ranges (1~200, 300~500 and 300~400).
#'  if `POS` is `NULL`, `start` and `end` are needed, eg.:
#'  `start = c(1, 30)` and `end = c(200, 150)` indicates 2 ranges
#'  (1~200 and 30~150)
#' @param override whether override existed file or not, default as `TRUE`.
#' @details
#' This package import VCF files with 'vcfR' which is more efficient to
#' import/manipulate VCF files in 'R'. However, import a large VCF file is time and
#' memory consuming. It's suggested that filter/extract variants in target
#' range with `filterLargeVCF()`.
#'
#' When filter/extract multi genes/ranges, the parameter of `Chr` and `POS`
#' must have equal length. Results will save to a single file if the user
#' provide a single file path or save to multiple VCF file(s) when a equal length
#' vector consist with file paths is provided.
#'
#' However, if you have hundreds gene/ranges need to extract from very
#' large VCF file(s), it's prefer to process with other linux tools in a script
#' on server, such as: 'vcftools' and 'bcftools'.
#' @examples
#' \donttest{
#'  # The filteration of small vcf should be done with `filter_vcf()`.
#'  # however, here, we use a mini vcf instead just for example
#'
#'  vcfPath <- system.file("extdata", "var.vcf.gz", package = "geneHapR")
#'
#'  oldDir <- getwd()
#'  setwd(tempdir())
#'  # extract a single gene/range from large vcf
#'  filterLargeVCF(VCFin = vcfPath, VCFout = "filtered.vcf.gz",
#'                 Chr = "scaffold_1", POS = c(4300,500), override = TRUE)
#'
#'  # extract multi genes/ranges from large vcf
#'  filterLargeVCF(VCFin = vcfPath,
#'                 VCFout = c("filtered1.vcf.gz",
#'                            "filtered2.vcf.gz",
#'                            "filtered3.vcf.gz"),
#'                 Chr = rep("scaffold_1", 3),
#'                 POS = list(c(4300, 5000),
#'                            c(5000, 6000),
#'                            c(5000, 7000)),
#'                 override = TRUE)
#'
#' setwd(oldDir)
#' }
#' @return No return value
#' @export
filterLargeVCF <- function(VCFin = VCFin,
                           VCFout = VCFout,
                           Chr = Chr,
                           POS = NULL,
                           start = start,
                           end = end,
                           override = TRUE)
{
    if (missing(Chr))
        stop("Chr is missing")
    if (is.null(POS) && (missing(start) | missing(end)))
        stop("POS is null, please provide start(s) and end(s)")
    if (!file.exists(VCFin))
        stop("Can't find ", VCFin, ", please check input")
    if (!override)
        if (TRUE %in% file.exists(VCFout))
            stop(paste(VCFout[file.exists(VCFout)], collapse = " "),
                 " existed, please check 'VCFout'")

    if(is.null(POS)){
        if(length(start) != length(end) | length(start) != length(POS))
            stop("length of 'POS', 'start', 'end' should be equal")
        if(!is.numeric(start) | !is.numeric(end))
            stop("'start' and 'end' should be numeric")
        if(length(Chr) == 1) POS <- c(start, end) else{
            POS <- list()
            for(i in seq_len(length(start))) POS <- c(POS, c(start[i], end[i]))
        }
    } else {
        if (length(Chr) == 1)
            if (length(POS) != 2 | !is.numeric(POS)) {
                stop("'POS' should be interger vector consist with start and end position")
                if (!is.numeric(POS))
                    stop("'POS' should be a numeric vector")
            }
        if (length(Chr) > 1)
            if (length(POS) != length(Chr) | ! inherits(POS,"list"))
                stop("'POS' should be a list have length equal with 'Chr'")
    }
    # number of conditions and file paths
    nCond <- length(Chr)

    if (nCond == 1) {
        filterLargeVCF_One(
            VCFin = VCFin,
            VCFout = VCFout,
            Chr = Chr,
            POS = POS,
            override = override
        )
    } else if (nCond > 1) {
        if (length(VCFout) != nCond)
            stop("length of 'VCFout' should be equal with 'Chr' and 'POS'")
        filterLargeVCF_Multi(
            VCFin = VCFin,
            VCFout = VCFout,
            Chr = Chr,
            POS = POS,
            override = override
        )
    }

}



filterLargeVCF_One <- function(VCFin = VCFin,
                               VCFout = VCFout,
                               Chr = Chr,
                               POS = POS,
                               override = override) {
    start <- min(POS)
    end <- max(POS)
    t0 <- proc.time()
    # set gz
    if (endsWith(VCFin, "gz"))
        iz <- gzcon(file(VCFin, "rb"))
    else
        if (endsWith(VCFin, "vcf"))
            iz <- file(VCFin, "rb")
    else
        stop("Input should with surfix 'vcf' or '.vcf.gz'")
    on.exit(close(oz))
    if (endsWith(VCFout, "gz"))
        oz <- gzcon(file(VCFout, "wb"))
    else
        if (endsWith(VCFout, "vcf"))
            oz <- file(VCFout, "wb")
    else
        stop("Output should with surfix 'vcf' or '.vcf.gz'")
    on.exit(close(iz), add = TRUE)

    nl <- 0

    # process meta information
    message("Processing meta information")
    while (TRUE) {
        l <- readLines(iz, n = 1)
        if (substr(l, 1, 1) != "#")
            break
        writeLines(l, con = oz)
        nl <- nl + 1
    }

    # process variant information
    subn <- nchar(Chr) + nchar(max(POS)) + 1
    while (TRUE) {
        if (length(l) == 0)
            break
        if (nl %% 1e4 == 0)
            message("Processed ", nl, " lines,", " used ", proc.time() - t0, )
        lc <- unlist(strsplit(substr(l, 1, subn), split = "\t"))
        Chrl <- lc[1]
        POSl <- as.numeric(lc[2])
        if (Chrl == Chr & POSl >= start & POSl <= end)
            writeLines(l, con = oz)

        l <- readLines(iz, n = 1)
        nl <- nl + 1

    }
    message("Processed ", nl, " lines. \nExit")
}



filterLargeVCF_Multi <- function(VCFin = VCFin,
                                VCFout = VCFout,
                                Chr = Chr,
                                POS = POS,
                                override = override) {
    nC <- length(Chr)
    nF <- length(VCFout)
    if (nC != nF)
        stop("Length of 'Chr'is not equal to 'VCFout'")

    # set conditions
    for (i in seq_len(nC)) {
        assign(paste0("start", i), min(POS[[i]]))
        assign(paste0("end", i), max(POS[[i]]))
        assign(paste0("Chr", i), min(Chr[i]))
    }
    t0 <- proc.time()

    # set gzcon
    if (endsWith(VCFin, "gz"))
        iz <- gzcon(file(VCFin, "rb"))
    else
        if (endsWith(VCFin, "vcf"))
            iz <- file(VCFin, "rb")
    else
        stop("Input should with surfix 'vcf' or '.vcf.gz'")
    on.exit(close(oz))

    for (i in seq_len(nF)) {
        if (endsWith(VCFout, "gz")) {
            assign(paste0("oz", i), gzcon(file(VCFout, "wb")))
        } else if (endsWith(VCFout, "vcf")) {
            assign(paste0("oz", i), file(VCFout, "wb"))
        } else {
            stop("Output should with surfix 'vcf' or '.vcf.gz'")
        }
    }
    for (i in seq_len(nF))
        on.exit(expr =  close(get(paste0("oz", i))),
                add = TRUE)

    # init number of rows
    nl <- 0

    # process meta information
    message("Processing meta information")
    while (TRUE) {
        l <- readLines(iz, n = 1)
        if (substr(l, 1, 1) != "#")
            break

        # save lines
        for (i in seq_len(nF))
            writeLines(l, con = get(paste0("oz", i)))

        nl <- nl + 1
    }

    # process variant information
    subn <- nchar(Chr) + nchar(max(POS)) + 10
    while (TRUE) {
        if (length(l) == 0)
            break
        if (nl %% 1e4 == 0)
            message("Processed ", nl, " lines,", " used ", proc.time() - t0)
        lc <- unlist(strsplit(substr(l, 1, subn), split = "\t"))
        Chrl <- lc[1]
        POSl <- as.numeric(lc[2])

        # save lines
        for (i in seq_len(nF)) {
            if (Chrl == get(paste0("Chr", i)) &
                POSl >= get(paste0("start", i)) &
                POSl <= get(paste0("end", i)))
                writeLines(l, con = get(paste0("oz", i)))
        }

        l <- readLines(iz, n = 1)
        nl <- nl + 1
    }

    # close con
    message("Processed ", nl, " lines. \nExit")
}
