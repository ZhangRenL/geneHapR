#' @title generate haplotype net relationshop with haplotype result
#' @description get_hapNet computes a haplotype network with haplotype summary result.
#' @seealso plot, plothapNet, hap_result
#' @usage get_hapNet(hapResult, accGroup = accGroup, groupName = groupName)
#' @importFrom pegas haploNet
#' @importFrom stringdist stringdist
#' @references
#' Article{,
#' title = {The stringdist package for approximate string matching},
#' author = {M.P.J. {van der Loo}},
#' year = {2014},
#' journal = {The {R} {J}ournal},
#' volume = {6},
#' issue = {1},
#' url = {https://CRAN.R-project.org/package=stringdist},
#' pages = {111-122},
#' }
#' Article{,
#'     title = {pegas: an {R} package for population genetics with an integrated--modular approach},
#'   author = {E. Paradis},
#'     journal = {Bioinformatics},
#'     year = {2010},
#'     volume = {26},
#'     pages = {419-420},
#' }
#' @param hapResult hapResult
#' @param accGroup data.frame, specified groups of each accession.
#' Used for pie plot. If missing, pie will not draw in plotHapNet.
#' Or you can supplied a hap.group mattrix with plot(hapNet, pie = hap.group).
#' @param groupName the group name used for pie plot, should be one of accGroup colnames, defalt as the first column name
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
#' @importFrom graphics legend
#' @usage
#' plotHapNet(hapNet,
#'            size = "freq", scale = TRUE, scale.ratio = 1, cex = 0.8,
#'            col.link = 1, link.width = 1,
#'            show.mutation = 1, lwd = 1,
#'            pieCol = pieCol, pieData = hapGroup,
#'            addLegend = TRUE,
#'            legendPosition = "left", ...)
#' @param hapNet an object of class "haploNet".
#' @param size a numeric vector giving the diameter of the circles representing the haplotypes: this is in the same unit than the links and eventually recycled.
#' @param scale.ratio the ratio of the scale of the links representing the number of steps on the scale of the circles representing the haplotypes. It may be needed to give a value greater than one to avoid overlapping circles.
#' @param col.link a character vector specifying the colours of the links; eventually recycled.
#' @param link.width a numeric vector giving the width of the links; eventually recycled.
#' @param show.mutation an integer value: if 0, nothing is drawn on the links; if 1, the mutations are shown with small segments on the links; if 2, they are shown with small dots; if 3, the number of mutations are printed on the links.
#' @param lwd a numeric vector giving the width of the links; eventually recycled.
#' @param pieCol color vector
#' @param pieData a matrix used to draw pie charts for each haplotype; its number of rows must be equal to the number of haplotypes.
#' @param addLegend a logical specifying whether to draw the legend, or a vector of length two giving the coordinates where to draw the legend; FALSE by default. If TRUE, the user is asked to click where to draw the legend.
#' @param cex character expansion factor relative to current par("cex").
#' Used for text, and provides the default for pt.cex.
#' @param pieCol colors
#' @param scale TRUE
#' @param legendPosition "left", "right"
#' @param ... pass to plot function
#' @export
plotHapNet <- function(hapNet,
                       size = "freq", scale = TRUE, scale.ratio = 1, cex = 0.8,
                       col.link = 1, link.width = 1,
                       show.mutation = 1, lwd = 1,
                       pieCol = pieCol, pieData = hapGroup,
                       addLegend = TRUE, legendPosition = "left", ...){
    if(!inherits(hapNet, "haploNet"))
        stop("'hapNet' must be of 'haploNet' class")
    if(missing(pieData)){
        hapGroup <- attr(hapNet, "hapGroup")
    } else {
        hapGroup <- pieData
    }
    if(size == "freq") size <- attr(hapNet, "freq") else
        if(!is.numeric(size))
            stop("'size' should be 'freq' or a given vector")

    if(scale) size <- (log10(size + 1) * 10) %/% 1


    if(!is.null(hapGroup)){
        if(missing(pieCol))
            pieCol <- rainbow(ncol(hapGroup))
        plot(hapNet, col.link = col.link,threshold = 10,
             size = size, scale.ratio = scale.ratio, cex = cex,
             show.mutation = show.mutation, lwd = link.width,
             bg = pieCol, pie = hapGroup)
        if(addLegend)
            legend(x = legendPosition, legend = colnames(hapGroup), fill = pieCol,cex = 0.6,...)


    } else {
        if(missing(pieCol))
            pieCol <- "grey90"
        plot(hapNet, col.link = col.link,
             size = size, scale.ratio = scale.ratio, cex = cex,
             show.mutation = show.mutation, lwd = link.width,legend = addLegend,
             ...)
    }
}


#' @title as.haplotype
#' @description convert hapSummary into haplotype (pegas)
#' @note It's not advised for hapSummary with indels, due to indels will replaced by equal length of SNPs.
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
#' @param type return a charactot vector or a DNAbin object
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

