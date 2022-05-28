#' @title filterVCFin
#' @param VCFin Path of input vcf file
#' @param VCFout Path of output vcf file
#' @param Chr CHROM name
#' @param POS a vector consist with start and end postion
#' @param overWrite Overwrite output file or not, defalt is TRUE
#' @export
#' @usage
#'  bigVcfFilter(VCFin = "", VCFout = "",
#'               Chr = Chr, POS = POS, overWrite = TRUE)
#' @examples
#' \dontrun{
#'  bigVcfFilter(VCFin = "Ori.vcf.gz", VCFout = "filtered.vcf.gz",
#'               Chr = "scaffold_8", POS = c(19802,24501), overWrite = TRUE)
#'               }
bigVcfFilter <- function(VCFin = "", VCFout = "",
                         Chr = Chr, POS = POS, overWrite = TRUE)
{
    if(missing(Chr)) stop("Chr is missing")
    if(missing(POS)) stop("POS is missing")
    if(missing(subn)) subn <- nchar(Chr) + nchar(max(POS)) + 1
    if(!file.exists(VCFin)) stop("Can't find ", VCFin, ", please check input")
    if(!overWrite)
        if(file.exists(VCFout))
            stop(VCFout, " existed, please check inputs")
    start <- min(POS)
    end <- max(POS)
    t0 <- proc.time()
    # set gz
    if(endsWith(VCFin,"gz"))
        iz <- gzcon(file(VCFin, "rb")) else
            if(endsWith(VCFin,"vcf"))
                iz <- file(VCFin, "rb") else
                    stop("Input should with surfix 'vcf' or '.vcf.gz'")
    if(endsWith(VCFout,"gz"))
        oz <- gzcon(file(VCFout, "wb")) else
            if(endsWith(VCFout,"vcf"))
                oz <- file(VCFout, "wb") else
                    stop("Output should with surfix 'vcf' or '.vcf.gz'")

    nl <- 0

    message("Processing meta information")
    while(TRUE){
        l <- readLines(iz, n = 1)
        if(substr(l,1,1) != "#") break
        writeLines(l, con = oz)
        nl <- nl + 1
    }
    while(TRUE){
        if(length(l) == 0) break
        if(nl %% 1e4 == 0) message("Processed ",nl, " lines,", " used ", proc.time() - t0,)
        lc <- unlist(strsplit(substr(l,1,subn),split = "\t"))
        Chrl <- lc[1]
        POSl <- as.numeric(lc[2])
        if(Chrl == Chr & POSl >= start & POSl <= end)
                writeLines(l, con = oz)

        l <- readLines(iz, n = 1)
        nl <- nl + 1

    }
    message("Processed ",nl, " lines. \nExit")
    close(iz)
    close(oz)
}
