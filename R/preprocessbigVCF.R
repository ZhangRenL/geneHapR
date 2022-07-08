# File Checked
#' @name filterLargeVCF
#' @title Pre-process large vcf
#' @description
#' Filter/extract one or multiple gene(s)/range(s) from a large
#' `*.vcf/*.vcf.gz` file.
#' @usage
#'  filterLargeVCF(VCFin = "", VCFout = VCFout,
#'                 Chr = Chr, POS = POS, override = TRUE)
#' @param VCFin Path of input `*.vcf/*.vcf.gz` file.
#' @param VCFout Path or path vector of output `*.vcf/*.vcf.gz` file.
#' @param Chr a single `CHROM` name or `CHROM` names vector.
#' @param POS a vector consist with `start` and `end` position or a list
#'  with length equal to `Chr`, each element consist with `start` and `end`
#'  position.
#' @param override whether override existed file or not, default as `TRUE`.
#' @details
#' This package import vcfs with `vcfR` which is more efficient to
#' import/manipulate vcf in R. However, import a large vcf file is time and
#' memory consuming. It's suggested that filter/extract variants in target
#' range with `filterLargeVCF`.
#'
#' When filter/extract multi genes/ranges, the parameter of `Chr` and `POS`
#' must have equal length. Results will save to a single file if the user
#' provide a single file path or save to multiple vcf files when a equal length
#' vector consist with file paths is provided.
#'
#' However, if you have hundreds gene/ranges need to extract from a very
#' large VCF, it's prefer to process with other linux tools in a script
#' on server such as: `vcftools` and `bcftools`.
#' @examples
#' \dontrun{
#'  # extract a single gene/range from large vcf
#'  filterLargeVCF(VCFin = "Ori.vcf.gz", VCFout = "filtered.vcf.gz",
#'                 Chr = "scaffold_8", POS = c(19802,24501), override = TRUE)
#'
#'  # extract multi genes/ranges from large vcf
#'  filterLargeVCF(VCFin = "Ori.vcf.gz",
#'                 VCFout = c("filtered1.vcf.gz",
#'                            "filtered2.vcf.gz",
#'                            "filtered3.vcf.gz"),
#'                 Chr = c("scaffold_8",
#'                         "scaffold_8",
#'                         "scaffold_7"),
#'                 POS = list(c(19802,24501),
#'                            c(27341,28949),
#'                            c(38469,40344)),
#'                 override = TRUE)
#'
#' }
#' @export
filterLargeVCF <- function(VCFin = "",
                           VCFout = VCFout,
                           Chr = Chr,
                           POS = POS,
                           override = TRUE)
{
    if (missing(Chr))
        stop("Chr is missing")
    if (missing(POS))
        stop("POS is missing")
    if (!file.exists(VCFin))
        stop("Can't find ", VCFin, ", please check input")
    if (!override)
        if (TRUE %in% file.exists(VCFout))
            stop(paste(VCFout[file.exists(VCFout)], collapse = " "),
                 " existed, please check 'VCFout'")
    if (length(Chr) == 1)
        if (length(POS) != 2 | !is.numeric(POS)) {
            stop("'POS' should be interger vector consist with start and end position")
            if (!is.numeric(POS))
                stop("'POS' should be a numeric vector")
        }
    if (length(Chr) > 1)
        if (length(POS) != length(Chr) | class(POS) != "list")
            stop("'POS' should be a list have length equal with 'Chr'")

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
    if (endsWith(VCFout, "gz"))
        oz <- gzcon(file(VCFout, "wb"))
    else
        if (endsWith(VCFout, "vcf"))
            oz <- file(VCFout, "wb")
    else
        stop("Output should with surfix 'vcf' or '.vcf.gz'")

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
    close(iz)
    close(oz)
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

    for (i in seq_len(nF)) {
        if (endsWith(VCFout, "gz")) {
            assign(paste0("oz", i), gzcon(file(VCFout, "wb")))
        } else if (endsWith(VCFout, "vcf")) {
            assign(paste0("oz", i), file(VCFout, "wb"))
        } else {
            stop("Output should with surfix 'vcf' or '.vcf.gz'")
        }
    }


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
    close(iz)
    for (i in seq_len(nF))
        close(get(paste0("oz", i)))
}
