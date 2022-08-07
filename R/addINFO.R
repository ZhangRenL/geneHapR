#' @name addINFO
#' @title Summary hap results
#' @description add annotations to INFO fields used for `plotHapTable()`
#' @usage
#' addINFO(hap,
#'         tag = "", values = values,
#'         replace = FALSE, sep = ";")
#' @examples
#'
#' data("geneHapR_test")
#'
#' # length of values must be equal with number of sites in hap result
#' values <- paste0("newInfo",c(1:9))
#' hapResult <- addINFO(hapResult, tag = "new", values = values, replace = TRUE)
#'
#' @seealso
#' \code{\link[geneHapR:plotHapTable]{plotHapTable()}}
#' @param hap object of hapResult or hapSummary class
#' @param tag tag names, usually is a single word used before "="
#' @param values annotation for each site.
#' Length of values must be equal with sites in hapResult
#' @param replace whether replace origin INFOs in hapResult or not.
#' Default as FALSE
#' @inheritParams base::paste
#' @return object of hapSummary or hapResult class with added/replaced INFOs
#' @export
addINFO <- function(hap,
                    tag = "",
                    values = values,
                    replace = FALSE,
                    sep = ";") {
    # get POS in hap
    POS <- t(hap[hap$Hap == "POS", ])[, 1]
    POS <- suppressWarnings(as.numeric(names(POS)))
    probe <- !is.na(POS)
    POS <- POS[probe]
    if (length(values) != length(POS))
        stop("Length of 'values' should be equal with sites")

    # set tag_value pairs
    probe_na <- is.na(values)
    values <- paste0(tag, "=", values)
    values[probe_na] <- ""

    if (replace) {
        INFO <- values
    } else {
        # get INFO in hap
        INFO <- t(hap[hap$Hap == "INFO", ])[, 1]
        INFO <- INFO[probe]
        INFO[is.na(INFO)] <- ""

        # set INFOs
        for (i in seq_len(length(INFO))) {
            if(nchar(values[i]) == 0) next
            if (nchar(INFO[i]) == 0) {
                INFO[i] <- values[i]
            } else {
                INFO[i] <- paste(INFO[i], values[i], sep = sep)
            }
        }
    }

    n <- ncol(hap) - length(INFO) - 1
    if(n == 1) {
        hap[hap$Hap == "INFO", ] <- c("INFO", INFO, "")
    } else if(n == 2){
        hap[hap$Hap == "INFO", ] <- c("INFO", INFO, "", NA)
    }
    return(hap)
}

#' @name addINFO
#' @usage
#' sites(hap)
#' @examples
#' \donttest{
#' # check how many sites were concluded in hapResult/hapSummary
#' sites(hapResult)
#' }
#' @seealso
#' \code{\link[geneHapR:plotHapTable]{plotHapTable()}}
#' @param hap object of `hapResult` or `hapSummary` class
#' @export
sites <- function(hap){
    # get POS in hap
    POS <- t(hap[hap$Hap == "POS", ])[, 1]
    POS <- suppressWarnings(as.numeric(names(POS)))
    probe <- !is.na(POS)
    POS <- POS[probe]
    length(POS)
}
