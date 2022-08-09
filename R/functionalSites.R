#' @name siteEff
#' @title Find site functional effective
#' @param hap object of "hapResult" class
#' @param pheno phenotype data, with column names as pheno name
#' and row name as accessions.
#' @param phenoNames pheno names used for analysis, if missing,
#' will use all pheno names in `pheno`
#' @return matrix, with column name are pheno names and
#'  row name are site position
#' @importFrom stats t.test
#' @usage
#' siteEff(hap, pheno, phenoNames)
#' @examples
#' \donttest{
#' data("geneHapR_test")
#'
#' # calculate site functional effect
#' siteEFF <- siteEff(hapResult, pheno, names(pheno))
#' plotEff(siteEFF, gff = gff, Chr = "scaffold_1")
#' }
#' @export
siteEff <- function(hap, pheno, phenoNames){
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
        cat("\n\t", phynoname)
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



#' @title plotEff
#' @name plotEff
#' @importFrom graphics par strwidth rect points
#' @usage
#' plotEff(siteEff, gff = gff,
#'         Chr = Chr, start = start, end = end,
#'         showType = c("five_prime_UTR", "CDS", "three_prime_UTR"),
#'         CDS.height = 1, cex = 0.1, col = col, pch = 20,
#'         title = title, legend.cex = 0.8, legend.ncol = legend.ncol,
#'         markMutants = TRUE, mutants.col = 1, mutants.type = 1)
#' @inherit siteEff examples
#' @param siteEff matrix, column name are pheno names and row name are site position
#' @param gff gff
#' @param Chr the chromosome name
#' @param start start postion
#' @param end end position
#' @param showType character vector, eg: "CDS", "five_prime_UTR",
#' "three_prime_UTR"
#' @param CDS.height numeric indicate the height of CDS in gene model,
#' range: `[0,1]`
#' @param cex a numeric control the size of point
#' @param col vector controls points color, see
#' \code{\link[graphics:points]{points()}}
#' @param pch vector controls points type, see
#' \code{\link[graphics:par]{par()}}
#' @param title main title
#' @param legend.cex a numeric control the legend size
#' @param legend.ncol the number of columns in which to set the legend items
#' @param markMutants whether mark mutants on gene model, default as `TRUE`
#' @param mutants.col color of lines which mark mutants
#' @param mutants.type a vector of line types
#' @export
plotEff <- function(siteEff, gff = gff,
                    Chr = Chr, start = start, end = end,
                    showType = c("five_prime_UTR", "CDS", "three_prime_UTR"),
                    CDS.height = 1, cex = 0.1, col = col, pch = 20,
                    title = title, legend.cex = 0.8, legend.ncol = legend.ncol,
                    markMutants = TRUE, mutants.col = 1, mutants.type = 1){
    EFF <- as.matrix(siteEff)
    POS <- suppressWarnings(as.numeric(row.names(EFF)))

    if(missing(Chr))
        stop("Chr is missing")
    if(missing(start))
        start <- min(POS) - 0.05 * diff(range(POS))
    if(missing(end))
        end <- max(POS) + 0.05 * diff(range(POS))

    # get GFF ranges for display
    gr <- GenomicRanges::GRanges(seqnames = Chr,
                                 ranges = IRanges::IRanges(start = start,
                                                           end = end))
    gff <- gff[IRanges::`%over%`(gff, gr)]
    gff <- gff[gff$type %in% showType]


    # reset of par
    oldPar.fig <- par("fig")
    oldPar.mar <- par("mar")
    on.exit(par(fig = oldPar.fig, mar = oldPar.mar))

    # plot genemodel
    # set of fig.h
    Parents <- unique(unlist(gff$Parent))
    nsplicement <- length(Parents)
    fig.h <- ifelse(nsplicement >= 5, 0.5, 0.1 * (1 + nsplicement))

    # set of par
    par.mar <- oldPar.mar
    par.mar[3] <- 0
    par(fig = c(0, 1, 0, fig.h), mar = par.mar)


    # just plot
    plot(x = c(start, end), y = c(0, nsplicement * 1.1),
         yaxt = "n", type = "n", xlab="", ylab ="",
         frame.plot = FALSE)

    # add legend
    xy <- par("usr")
    if(missing(legend.ncol))
        legend.ncol <- ifelse(nsplicement <= 3, nsplicement, 3)

    legend(x = 0.5 * (xy[1] + xy[2]), y = xy[3] - 4 * strheight(""), legend = Parents,
           fill = rainbow(length(Parents)), xjust = 0.5,cex = legend.cex,
           ncol = legend.ncol, xpd = TRUE)

    ln <- -0.6
    n <- 1
    for(s in Parents){
        gffs <- gff[unlist(gff$Parent) == s]
        anno <- ifelse(gffs@strand[1] == "-", "3'<-5'", "5'->3'")


        ln <- ln + 1.1
        lines(c(start,end),c(ln,ln))
        text(start - strwidth(anno), ln, anno, xpd = T )
        s.col <- rainbow(nsplicement)[n]
        n <- n + 1
        for(i in seq_len(length(gffs))){
            gffi <- gffs[i]
            h <- ifelse(gffi$type == "CDS", CDS.height, CDS.height * 0.5) * 0.5
            xl <- gffi@ranges@start
            xr <- xl + gffi@ranges@width - 1
            rect(xleft = xl, xright = xr, ybottom = ln - h, ytop = ln + h, col = s.col)
        }
    }

    # plot EFFs
    # set of mar
    par.mar <- oldPar.mar
    par.mar[1] <- 0

    # set of par and plot frame
    par(fig = c(0,1,fig.h + 0.01, 1), mar = par.mar, new = TRUE)
    plot(x = POS[1], y = EFF[1,1], type = "n",
         xlim = c(start, end), ylim = c(0, max(EFF)),
         col = 3, cex = 0.5,
         xaxt = "n", xlab = "", ylab = expression("-log"[10]~italic(p)~"Value"))
    cols <- rainbow(length(POS))

    if(missing(col)) col <- seq_len(nrow(EFF)) else
        col <- if(length(col) == 1) rep(col, nrow(EFF)) else col

    if(missing(pch)) pch <- 20
    pch <- if(length(pch) != nrow(EFF)) rep(pch, nrow(EFF)) else pch


    # plot points indicate EFFs
    for(i in seq_len(nrow(EFF))){
        points(x = rep(POS[i], ncol(EFF)),
               y = EFF[i,],
               cex = cex, col = col[i], pch = pch[i])
    }

    # markMutants
    if(markMutants){
        par.mar <- oldPar.mar
        par.mar[3] <- 0
        par(fig = c(0, 1, 0.01, fig.h + 0.01), mar = par.mar, new = TRUE)
        plot(start, xlim = c(start, end), ylim = c(0, nsplicement * 1.1),
             type = "n", xaxt = "n", yaxt = "n",
             xlab = "", ylab = "", frame.plot = FALSE)
        for(pos in POS){
            lines(c(pos,pos), c(0.4, ln+2.1),
                  col = mutants.col, lty = mutants.type)
        }
    }

    # add title
    if(!missing(title))
        title(main = title)
}
