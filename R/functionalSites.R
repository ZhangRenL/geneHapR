#' @name siteFunEff
#' @title Find site functional effective
#' @param hap object of "hapResult" class
#' @param pheno phenotype data, with column names as pheno name
#' and row name as accessions.
#' @param phenoNames pheno names used for analysis, if missing,
#' will use all pheno names in `pheno`
#' @return data.frame, with column name as pheno name and
#'  row name as site position
#' @importFrom stats t.test
#' @usage
#' siteFunEff(hap, pheno, phenoNames)
#' @examples
#' data("geneHapR_test")
#'
#' # calculate site functional effect
#' funResults <- siteFunEff(hapResult, pheno, names(pheno))
#' plotSiteFunEff(funResults, facet = FALSE)
#' @export
siteFunEff <- function(hap, pheno, phenoNames){
    if(missing(phenoNames)) phenoNames <- names(pheno)
    if(!inherits(hap, "hapResult"))
        stop("hap should be object of 'hapResult' class")

    # get positions
    POS <- hap[hap$Hap == "POS",]
    POS <- suppressWarnings(as.numeric(POS))
    POS <- POS[! is.na(POS)]

    # extract genotype data
    hapData <- hap[! hap$Hap %in% c("POS","CHR","ALLELE","INFO"),]

    # get accession list
    accessions <- hapData[,names(hapData) == "Accession"]

    # preset of results
    results <- data.frame()

    # processing
    for(phynoname in phenoNames){
        res <- c()
        for(pos in POS){

            # get alleles
            alleles <- hapData[,as.character(pos)]
            Als <- unique(alleles)
            Aln <- length(unique(Als))

            # get accessions of each genotype
            for(i in seq_len(Aln)){
                probe <- c(alleles == Als[i])
                assign(paste0("acc",i), accessions[probe])
            }
            p <- c()
            if(Aln > 1)
                for(i in seq_len(Aln)){
                    for(j in rev(seq_len(Aln))){
                        if(i >= j) next
                        acci <- get(paste0("acc",i))
                        accj <- get(paste0("acc",j))
                        phenoi <- pheno[acci, phynoname]
                        phenoj <- pheno[accj, phynoname]
                        pij <- try(t.test(phenoi, phenoj), silent = TRUE)
                        if(inherits(pij, "htest")){
                            pij <- pij$p.value
                        } else {
                            pij <- 1
                        }
                        p <- c(p, pij)
                    }
                }
            p <- min(p)
            res <- c(res, p)
        }
        results <- rbind(results, res)
    }
    names(results) <- POS
    results <- -log10(results)
    row.names(results) <- phenoNames
    # results <- cbind(pheno = phenoNames, results)
    return(t(results))
}



#' @name siteFunEff
#' @param results results of site functional effect analysis
#' @param facet use facet or not, default as `FALSE`
#' @inheritParams ggplot2::labs
#' @inheritParams ggplot2::facet_wrap
#' @import ggplot2
#' @importFrom rlang .data
#' @usage
#' plotSiteFunEff(results,
#'                title = title,
#'                caption = caption,
#'                facet = FALSE,
#'                ...)
#' @export
plotSiteFunEff <- function(results, title = title, caption = caption, facet = FALSE, ...){
    data <- suppressMessages(reshape2::melt(results))
    colnames(data) <- c("Position", "pheno", "value")
    data$value <- round(data$value, digits = 2)
    p <- ggplot2::ggplot(data = data,
                         mapping = aes_(x = ~ Position,
                                        y = ~ value,
                                        label = ~ value,
                                        color = ~ pheno))
    p <- p + ggplot2::geom_line() +
        ggplot2::xlab("Position") +
        ggplot2::ylab(expression(-log[10](italic(P)~value))) +
        theme_bw()
    if(!missing(title)) p <- p + ggplot2::labs(title = title, caption = caption)
    if(facet) p <- p + facet_wrap(facets = "pheno", ncol = 1)
    p
}
