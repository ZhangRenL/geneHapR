#' @title get_hapNet
#' @description calculate hap net using hap result
#' @usage get_hapNet(hapResult, accGroup = accGroup, groupName = groupName)
#' @importFrom pegas haploNet
#' @importFrom stringdist stringdist
#' @param hapResult hapResult
#' @param accGroup data.frame specified groups of each accessions,
#' used for pie plot. If missing, pie will not draw in plotHapNet.
#' Or you can supplied a hap.group mattrix with plot(hapNet, pie = hap.group).
#' @param groupName one colname of accGroup
#' @return haplonet class
#' @export
get_hapNet <- function(hapResult, accGroup = accGroup, groupName = groupName){
    if(!inherits(hapResult,"hapSummary"))
        stop("'hapResult' should be of 'hapSummary' class")
    d <- getStringDist(hapResult)
    hapNet <- pegas::haploNet(as.haplotype(hapResult), d)
    if(!missing(accGroup)){
        if(missing(groupName))
            groupName <- colnames(accGroup)[1]
        hapGroup <- getHapGroup(hapResult,
                               accGroup = accGroup,
                               groupName = groupName)
        attr(hapNet, "hapGroup") <- hapGroup
    }
    return(hapNet)
}


#' @title plotHapNet
#' @param lenendPOS one of "bottomright", "bottom", "bottomleft", "left",
#' "topleft", "top", "topright", "right" and "center".
#' @param size circle size
#' @param cex character expansion factor relative to current par("cex").
#' Used for text, and provides the default for pt.cex.
#' @param pieCol colors
#' @param ... pass to plot function
#' @export
plotHapNet <- function(hapNet,
                       size = "freq", scale.ratio = 1, cex = 0.8,
                       col.link = 1, link.width = link.width,
                       show.mutation = 1, lwd = 1,
                       pieCol = pieCol, pie = hapGroup,
                       legendPOS = "left",  ...){
    if(!inherits(hapNet, "haploNet"))
        stop("'hapNet' must be of 'haploNet' class")
    hapGroup <- attr(hapNet, "hapGroup")
    if(size == "freq") size <- attr(hapNet, "freq") else
        if(!is.numeric(size))
            stop("'size' should be 'freq' or a given vector")


    if(!is.null(hapGroup)){
        if(missing(pieCol))
            pieCol <- rainbow(ncol(hapGroup))
        plot(hapNet, col.link = col.link,
             size = size, scale.ratio = scale.ratio, cex = cex,
             show.mutation = show.mutation, lwd = link.width,
             bg = pieCol, pie = hapGroup, legend = TRUE, ...)

        #legend(x = legendPOS, legend = colnames(hapGroup), col = pieCol)

    } else {
        if(missing(pieCol))
            pieCol <- "grey90"
        plot(hapNet, col.link = col.link,
             size = size, scale.ratio = scale.ratio, cex = cex,
             show.mutation = show.mutation, lwd = link.width,legend = TRUE,
             ...)
    }
}

# plotHapNet(hapNet,legendPOS = "left")
#
# if (legend[1]) {
#     if (is.logical(legend)) {
#         cat("Click where you want to draw the legend")
#         xy <- unlist(locator(1))
#         cat("\nThe coordinates x = ", xy[1], ", y = ", xy[2], " are used\n", sep = "")
#     } else {
#         if (!is.numeric(legend) || length(legend) < 2)
#             stop("wrong coordinates of legend")
#         xy <- legend
#     }
#     if (length(SZ <- unique(size)) > 1) {
#         SZ <- unique(c(min(SZ), floor(median(SZ)), max(SZ)))
#         SHIFT <- max(SZ) / 2
#         vspace <- strheight(" ")
#         if (any(shape == "circles")) {
#             for (sz in SZ) {
#                 seqx <- seq(-sz / 2, sz / 2, length.out = 100)
#                 seqy <- sqrt((sz /2)^2 - seqx^2)
#                 seqx <- seqx + xy[1] + SHIFT
#                 seqy <- xy[2] + seqy - SHIFT
#                 lines(seqx, seqy)
#                 text(seqx[100], seqy[100], sz, adj = c(0.5, 1.1))
#             }
#             xy[2] <- xy[2] - SHIFT - 2 * vspace # update 'y'
#         }
#         if (any(shape == "squares") || any(shape == "diamonds")) {
#             sqrtPIon4 <- sqrt(pi) / 4
#             ## center of the squares:
#             orig.x <- xy[1] + max(SZ) * sqrtPIon4
#             orig.y <- xy[2] - max(SZ) * sqrtPIon4
#             for (sz in SZ) {
#                 TMP <- sz * sqrtPIon4
#                 lines(orig.x + c(-TMP, -TMP, TMP, TMP),
#                       orig.y + c(0, TMP, TMP, 0))
#                 text(orig.x + TMP, orig.y, sz, adj = c(0.5, 1.1))
#             }
#             xy[2] <- xy[2] - SHIFT - 2 * vspace # update
#         }
#     }
#     if (!is.null(pie)) {
#         nc <- ncol(pie)
#         ##co <- if (length(bg) == 1 && bg == "white") rainbow(nc) else rep(bg, length.out = nc)
#         co <- if (is.function(bg)) bg(nc) else rep(bg, length.out = nc)
#         w <- diff(par("usr")[3:4]) / 40
#         TOP <- seq(xy[2], by = -w, length.out = nc)
#         BOTTOM <- TOP + diff(TOP[1:2]) * 0.9
#         LEFT <- rep(xy[1], nc)
#         RIGHT <- LEFT + w
#         rect(LEFT, BOTTOM, RIGHT, TOP, col = co)
#         text(RIGHT, (TOP + BOTTOM) /2, colnames(pie), adj = -0.5)
#     }
# }


#' @title as.haplotype
#' @description calculat Hap net from a hapSummary object
#' @importFrom ape as.DNAbin
#' @param hapResult hapResult
#' @return haplotype class
#' @export
as.haplotype <- function(hapResult){
    if(!inherits(hapResult, "hapSummary"))
        stop("'hapResult' must be  of class 'hapSummary'")
    # get freq
    freq <- hapResult$freq
    names(freq) <- hapResult$Hap
    freq <- freq[!is.na(freq)]

    # get hap mattrix
    hapDNAset <- hap2string(hapResult = hapResult, type = "DNA")
    hapBin <- ape::as.DNAbin(hapDNAset)
    hapmatt <- unlist(as.character(as.matrix(hapBin)))

    # set as haplotype
    class(hapmatt) <- c("haplotype", "character")
    rownames(hapmatt) <- names(freq)
    N <- nrow(hapmatt)
    attr(hapmatt, "index") <- lapply(1:N, function(i) seq_len(freq[i]))
    return(hapmatt)
}


#' @title hap2DNAstring
#' @description convert hapsummary to BioStringsset class
#' @importFrom stringr str_detect str_pad str_length
#' @importFrom Biostrings DNAStringSet
#' @param hapResult hapResult
#' @export
hap2string <- function(hapResult, type = "DNA"){
    colNms <- colnames(hapResult)
    if("Accession" %in% colNms)
        hapResult <- hapResult[,colnames(hapResult) != 'Accession']
    if("freq" %in% colNms)
        hapResult <- hapResult[,colnames(hapResult) != 'freq']
    if(nrow(hapResult) <= 5) stop("There is only one Hap ?")
    meta <- hapResult[1:4,]
    hap <- hapResult[5:nrow(hapResult),]
    ALLELE <- meta[meta[,1] == "ALLELE",]

    # padding multiallelic indel sites
    if(type == "DNA"){
        multi_probe <- is.indel.allele(ALLELE)
        for(c in seq_len(ncol(hap))){
            if(multi_probe[c]){
                if(stringr::str_sub(ALLELE[c],2,2) == "/") side = "right" else
                    side = "left"
                maxLen <- max(stringr::str_length(hap[,c]))
                hap[,c] <- stringr::str_pad(hap[,c],
                                            width = maxLen,
                                            side = side,
                                            pad = "-")
            }
        }

        # conneting strings
        DNASeqs <- sapply(seq_len(nrow(hap)),
                          function(i) paste0(hap[i,-1],collapse = ""))
        names(DNASeqs) <- hap$Hap
        hapString <- Biostrings::DNAStringSet(DNASeqs, use.names = T, start = 1)
    } else {
        for(c in 2:ncol(hap)){
            ALc <- ALLELE[c]
            ALs <- unlist(stringr::str_split(ALc, "[,/]"))
            ALn <- LETTERS[seq_len(length(ALs))]
            names(ALn) <- ALs
            hap[,c] <- ALn[hap[,c]]
        }

        hapString <- sapply(seq_len(nrow(hap)),
                            function(i) paste0(hap[i,-1],collapse = ""))
        names(hapString) <- hap$Hap
    }
    return(hapString)
}


#' @importFrom stringdist stringdist
getStringDist <- function(hapResult){
    hapStrings <- hap2string(hapResult, type = "LETTER")
    n <- names(hapStrings)
    l <- length(hapStrings)
    d <- matrix(nrow = l, ncol = l,dimnames = list(n,n))
    for(i in seq_len(length(hapStrings))){
        for(j in seq_len(length(hapStrings))){
            d[i,j] <- stringdist::stringdist(hapStrings[i],
                                             hapStrings[j],
                                             method = "lv")
        }
    }
    d <- d[lower.tri(d)]
    attr(d, "Size") <- l
    attr(d, "Labels") <- n
    attr(d, "method") <- "lv"
    attr(d, "Diag") <- attr(d, "Upper") <- FALSE
    class(d) <- "dist"
    return(d)
}


getHapGroup <- function(hapResult,
                        accGroup = accGroup,
                        groupName = groupName){
    # get indvidual group number of each hap
    hap2acc <- attr(hapResult, "hap2acc")
    acc2hap <- names(hap2acc)
    names(acc2hap) <- hap2acc
    accGroup$Hap <- acc2hap[rownames(accGroup)]
    with(
        accGroup[,c("Hap", groupName)],
        table(hap = accGroup$Hap, group=accGroup[,groupName]))
}


is.indel.allele <- function(allele.vector){
    probe1 <- stringr::str_detect(allele.vector,"/")
    probe <- lapply(stringr::str_split(allele.vector,"[,/]"),
                    stringr::str_length)
    probe1 & unlist(lapply(probe, function(i) max(i)>1))
}

is.biallelic.allele <- function(allele.vector){
    probe1 <- stringr::str_detect(allele.vector,"/")
    !stringr::str_detect(allele.vector,",") & probe1
}

is.multiallelic.allele <- function(allele.vector){
    probe1 <- stringr::str_detect(allele.vector,"/")
    stringr::str_detect(allele.vector,",") & probe1
}

