# sites definition
homo <- c("A", "C", "G", "T")
names(homo) <- paste(homo, homo, sep = "/")
unknown <- c(".", "N")
names(unknown) <- paste(unknown, unknown, sep = "/")
hetero <- c("A|C",
            "A|G",
            "A|T",
            "C|A",
            "C|T",
            "C|G",
            "G|A",
            "G|C",
            "G|T",
            "T|A",
            "T|C",
            "T|G")
names(hetero) <- stringr::str_replace(hetero, "[|]", "/")


allS <- list(
    homo = homo,
    unknown = unknown,
    hetero = hetero,
    all = c(homo, hetero)
)


# B: C G T
# D: A G T
# H: A C T
# V：A C G
# N: A C G T
update_allS <- function(allS_new, REF = REF, ALT = ALT) {
    ALT <- stringr::str_split(ALT, pattern = ",")
    ALT <- unlist(ALT)
    # update homo
    for (l in c(REF, ALT)) {
        if (!l %in% allS_new$homo) {
            homonew <- l
            allS_new$homo <- c(allS_new$homo, homonew)
        }
    }
    allS_new$homo <- unique(allS_new$homo)
    names(allS_new$homo) <-
        paste(allS_new$homo, allS_new$homo, sep = "/")

    # update hetero
    for (i in c(REF, ALT)) {
        for (j in c(REF, ALT)) {
            if (i != j) {
                hetero_ij <- paste(i, j, sep = "|")
                allS_new$hetero <- c(allS_new$hetero, hetero_ij)
            }
        }
    }
    allS_new$hetero <- unique(allS_new$hetero)
    names(allS_new$hetero) <-
        stringr::str_replace(allS_new$hetero, "[|]", "/")

    # update allS_new$all
    allS_new$all <-
        c(allS_new$homo, allS_new$hetero, allS_new$unknown)

    return(allS_new)
}


is.indel.allele <- function(allele.vector) {
    probe1 <- stringr::str_detect(allele.vector, "/")
    probe <- lapply(stringr::str_split(allele.vector, "[,/]"),
                    stringr::str_length)
    probe <- probe1 & unlist(lapply(probe, function(i)
        max(i) > 1))
    probe[is.na(probe)] <- FALSE
    return(probe)
}


is.biallelic.allele <- function(allele.vector) {
    probe1 <- stringr::str_detect(allele.vector, "/")
    probe <- ! stringr::str_detect(allele.vector, ",") & probe1
    probe[is.na(probe)] <- FALSE
    return(probe)
}


is.multiallelic.allele <- function(allele.vector) {
    probe1 <- stringr::str_detect(allele.vector, "/")
    probe <- stringr::str_detect(allele.vector, ",") & probe1
    probe[is.na(probe)] <- FALSE
    return(probe)
}

#' @importFrom stats IQR quantile
removeOutlier <- function(x){
    outlier_limup <-
        quantile(x, 3 / 4, na.rm = TRUE, names = FALSE) +
        3 * IQR(x, na.rm = TRUE) # Q3+k(Q3-Q1)

    outlier_limdown <-
        quantile(x, 1 / 4, na.rm = TRUE, names = FALSE) -
        3 * IQR(x , na.rm = TRUE) # Q1-k(Q3-Q1)

    x[x >= outlier_limup | x <= outlier_limdown] = NA
    return(x)
}
