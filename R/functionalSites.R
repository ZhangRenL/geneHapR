#' @name siteEFF
#' @title Calculation of Sites Effective
#' @param hap object of "hapResult" class
#' @param pheno phenotype data, with column names as pheno name
#' and row name as accessions.
#' @param phenoNames pheno names used for analysis, if missing,
#' will use all pheno names in `pheno`
#' @param quality bool type, indicate whther the type of phenos are quality or
#' quantitative. Length of `quality` could be 1 or equal with length of
#' `phenoNames`. Default as `FALSE`
#' @param method character or character vector with length equal with
#' `phenoNames` indicate which method should be performed towards each
#' phenotype. Should be one of "t.test", "chi.test", "anova" and "auto".
#' Default as "auto", see details.
#' @param p.adj character, indicate correction method.
#' Could be "BH", "BY", "none"
#' @details
#' The site **EFF** was determinate by the phenotype difference between each
#' site geno-type.
#'
#' The *p* was calculated with statistical analysis method as designated by the
#' parameter `method`. If `method` set as "auto", then
#' chi.test will be
#' selected for quantity phenotype, eg.: color;
#' for quantity phynotype, eg.: height, with at least 30 observations per
#' geno-type and fit Gaussian distribution t.test will be performed or
#' anova will be performed.
#'
#'
#' @return a list containing two matrix names as "p" and "EFF",
#' with column name are pheno names and row name are site position.
#' The matrix names as "p" contains all *p*-value.
#' The matrix named as "EFF" contains scaled difference between each geno-types
#' per site.
#' @importFrom stats t.test chisq.test p.adjust shapiro.test
#' @usage
#' siteEFF(hap, pheno, phenoNames, quality = FALSE, method = "auto",
#'         p.adj = "none")
#' @examples
#' \donttest{
#' data("geneHapR_test")
#'
#' # calculate site functional effect
#' siteEFF <- siteEFF(hapResult, pheno, names(pheno))
#' plotEFF(siteEFF, gff = gff, Chr = "scaffold_1")
#' }
#' @export
siteEFF <- function(hap, pheno, phenoNames, quality = FALSE, method = "auto",
                    p.adj = "none"){
    if(missing(phenoNames)) phenoNames <- names(pheno)
    m <- "'quality' length should be equal with 'phenoNames'"
    if(length(quality) == 1)
        quality <- rep(quality, length(phenoNames)) else
            stopifnot("'quality' length should be equal with 'phenoNames'" =
                          length(quality[1:10]) == length(phenoNames))
    names(quality) <- phenoNames

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
    results.p <- data.frame()
    results.d <- data.frame()

    echo <- FALSE
    t <- Sys.time()
    # processing
    for(phynoname in phenoNames){
        # whether echo pheno name
        # if(echo) cat("\n\t", phynoname) else {
        #     if(t){
        #         if((Sys.time() - t) > 5)
        #             echo <- TRUE
        #         t <- FALSE
        #     }
        # }


        # is.quality
        is.quality <- quality[phynoname]

        # scale phenos if not quality
        pheno.n <- pheno[, phynoname]
        if(!is.quality) pheno.n <- pscale(pheno.n)
        names(pheno.n) <- rownames(pheno)

        # EFF and pValue calculate
        res.p <- c()
        res.d <- c()


        for(pos in POS){

            # get alleles
            alleles <- hapData[,as.character(pos)]
            Als <- unique(alleles)
            Aln <- length(unique(Als))

            # get accessions of each genotype
            phenos <- list()
            for(i in seq_len(Aln)){
                probe <- c(alleles == Als[i])
                accs <- accessions[probe]
                phenos[[i]] <- pheno.n[accs]
            }
            # test start


            if(method == "auto"){
                if(is.quality) {
                    # quility pheno
                    res.ps <- chisq.test.ps(phenos)
                } else { # quantity pheno
                    # shaporo.test
                    sha.p <- sapply(phenos,
                                    function(x) {
                                        x <- na.omit(x)
                                        if(length(x) < 3) return(0)
                                        if(length(x) > 5000) x <- sample(x, 5000)
                                        shapiro.test(x)$p.value
                                    }
                                    )
                    if(min(sha.p, na.rm = TRUE) >= 0.05){
                        # all sub data set fit normal distribution
                        res.ps <- t.test.ps(phenos)
                    } else {
                        # not all sub data set fit normal distribution
                        res.ps <- t.test.ps(phenos)
                    }
                }
            } else {
                switch (method,
                    "chisq.test" = chisq.test.ps(phenos),
                    "t.test" = t.test.ps(phenos)
                )
            }


            # test end
            p <- min(res.ps$p)
            d <- max(res.ps$d)
            res.p <- c(res.p, p)
            res.d <- c(res.d, d)

        }

        results.p <- rbind(results.p, res.p)
        results.d <- rbind(results.d, res.d)
    }

    if(p.adj != "none"){
        results.p <- matrix(p.adjust(as.matrix(results.p), method = p.adj),
                            nrow = nrow(results.p))
    }
    names(results.p) <- POS
    row.names(results.p) <- phenoNames
    names(results.d) <- POS
    results.d <- (results.d)
    row.names(results.d) <- phenoNames
    # results <- cbind(pheno = phenoNames, results)
    return(list(p = t(results.p), EFF = t(results.d)))
}


t.test.ps <- function(phenos){
    p <- c()
    d <- c()
    l = length(phenos)
    for(i in seq_len(l)){
        for(j in rev(seq_len(l))){
            if(i >= j) next
            phenoi <- phenos[[i]]
            phenoj <- phenos[[j]]

            # t.test or chisqure test or anova analysis
            pij.res <- try(t.test(phenoi, phenoj,
                              alternative = "two.sided"),
                       silent = TRUE)
            if(inherits(pij.res, "htest")){
                pij <- pij.res$p.value
                dij <- abs(diff(pij.res$estimate))
            } else {
                pij <- NA
                dij <- NA
            }

            p <- c(p, pij)
            d <- c(d, dij)
        }
    }
    list(p = p, d = d)
}


chisq.test.ps <- function(phenos){
    nms <- phenos %>%
        unlist() %>%
        na.omit() %>%
        unique() %>%
        as.character()
    nms <- nms[order(nms)]
    l <- length(phenos)
    ptable <- matrix(ncol = length(nms),
                     nrow = l,
                     dimnames = list(seq_len(l),
                                     nms))
    for(i in seq_len(l)) {
        freqi <- table(phenos[[i]])
        ptable[i,] <- freqi[nms]
    }
    ptable[is.na(ptable)] <- 0
    p <- chisq.test(t(ptable))
    p <- p$p.value
    ptable.f <- matrix(nrow = nrow(ptable), ncol = ncol(ptable))
    for(i in seq_len(nrow(ptable)))
        ptable.f[i,] <- ptable[i,]/sum(ptable[i,])
    d <- 0
    for(i in seq_len(ncol(ptable)))
        d <- (max(ptable.f[,i]) - min(ptable.f[,i])) / 2

    list(p = p, d = d)
}

# TODO
# add delta EFF plot function
#' @title plotEFF
#' @name plotEFF
#' @importFrom graphics par strwidth rect points
#' @importFrom grDevices heat.colors
#' @usage
#' plotEFF(siteEFF, gff = gff,
#'         Chr = Chr, start = start, end = end,
#'         showType = c("five_prime_UTR", "CDS", "three_prime_UTR"),
#'         CDS.height = 1, cex = 0.1, col = col, pch = 20,
#'         main = main, legend.cex = 0.8, legend.ncol = legend.ncol,
#'         markMutants = TRUE, mutants.col = 1, mutants.type = 1,
#'         ylab = "effect")
#' @inherit siteEFF examples
#' @param siteEFF matrix, column name are pheno names and row name are site position
#' @param gff gff
#' @param Chr the chromosome name
#' @param start start postion
#' @param end end position
#' @param showType character vector, eg.: "CDS", "five_prime_UTR",
#' "three_prime_UTR"
#' @param CDS.height numeric indicate the height of CDS in gene model,
#' range: `[0,1]`
#' @param cex a numeric control the size of point
#' @param col vector controls points color, see
#' \code{\link[graphics:points]{points()}}
#' @param pch vector controls points type, see
#' \code{\link[graphics:par]{par()}}
#' @param main main title
#' @param legend.cex a numeric control the legend size
#' @param legend.ncol the number of columns in which to set the legend items
#' @param markMutants whether mark mutants on gene model, default as `TRUE`
#' @param mutants.col color of lines which mark mutants
#' @param mutants.type a vector of line types
#' @param ylab character, yaxis label
#' @return No return value, called for side effects
#' @export
plotEFF <- function(siteEFF, gff = gff,
                    Chr = Chr, start = start, end = end,
                    showType = c("five_prime_UTR", "CDS", "three_prime_UTR"),
                    CDS.height = 1, cex = 0.1, col = col, pch = 20,
                    main = main, legend.cex = 0.8, legend.ncol = legend.ncol,
                    markMutants = TRUE, mutants.col = 1, mutants.type = 1,
                    ylab = "effect"){
                    # ylab = expression("-log"[10]~italic(p)~"Value")){
        EFF <- as.matrix(siteEFF$EFF)
    p <- siteEFF$p
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
    oldPar.mar.m <- oldPar.mar <- par("mar")
    on.exit(par(fig = oldPar.fig, mar = oldPar.mar))


    # plot genemodel
    # set of fig.h
    Parents <- unique(unlist(gff$Parent))
    nsplicement <- length(Parents)
    fig.h <- ifelse(nsplicement >= 5, 0.5, 0.1 * (1.2 + nsplicement))
    ln <- -0.6


    # set of par
    oldPar.mar.m[4] <- 0
    par.mar <- oldPar.mar.m
    par.mar[3] <- 0
    par(fig = c(0, 0.82, 0, fig.h), mar = par.mar)


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

    # markMutants
    if(markMutants){
        par.mar <- oldPar.mar.m
        par.mar[3] <- 0
        par(fig = c(0, 0.82, 0.01, fig.h + 0.01), mar = par.mar, new = TRUE)
        plot(start, xlim = c(start, end), ylim = c(0, nsplicement * 1.1),
             type = "n", xaxt = "n", yaxt = "n",
             xlab = "", ylab = "", frame.plot = FALSE)
        for(pos in POS){
            y.up <- ln + 1.1 * length(Parents) + 2.1
            lines(c(pos, pos), c(0.4, y.up),
                  col = mutants.col, lty = mutants.type)
        }
    }

    n <- 1
    for(s in Parents){
        gffs <- gff[unlist(gff$Parent) == s]
        anno <- ifelse(gffs@strand[1] == "-", "3'<-5'", "5'->3'")


        ln <- ln + 1.1
        lines(c(start,end),c(ln,ln))
        text(start - strwidth(anno), ln, anno, xpd = TRUE)
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
    par.mar <- oldPar.mar.m
    par.mar[1] <- 0

    # set of par and plot frame
    par(fig = c(0, 0.82, fig.h + 0.01, 1), mar = par.mar, new = TRUE)
    plot(x = POS[1], y = EFF[1,1], type = "n",
         xlim = c(start, end), ylim = c(0, max(EFF, na.rm = TRUE)),
         col = 3, cex = 0.5,
         xaxt = "n", xlab = "", ylab = ylab)
    cols <- rainbow(length(POS))

    if(missing(col)) col <- seq_len(nrow(EFF)) else
        col <- if(length(col) == 1) rep(col, nrow(EFF)) else col

    if(missing(pch)) pch <- 20
    pch <- if(length(pch) != nrow(EFF)) rep(pch, nrow(EFF)) else pch


    # plot points indicate EFFs
    # TODO
    # 1. add color for pValue
    # 2. height for EFF
    heatcols <- rev(heat.colors(1100))
    p <- -log10(p)
    p.max <- max(p)
    p.min <- min(p)
    cols <- round((p - p.min + 1) / (p.max - p.min + 1) * 1000)
    cols[,] <- heatcols[cols]
    for(i in seq_len(nrow(EFF))){
        points(x = rep(POS[i], ncol(EFF)),
               y = EFF[i,],
               cex = 1, col = cols[i,], pch = pch[i])
    }

    # add title
    if(!missing(main))
        title(main = main)


    # add color legend
    # set of mar
    par.mar <- oldPar.mar
    par.mar[1] <- 0
    par.mar[2] <- 1.5
    par(mar = par.mar, fig = c(0.82, 0.98, fig.h + 0.01, 1), new = TRUE)
    plot(y = 9,
         x = 1,
         xlim = c(0, 1),
         ylim = c(0, 1000),
         xaxt = 'n',
         yaxt = 'n',
         type = 'n',
         frame.plot  = FALSE)
    rect(xleft = rep(0, 1000),
         ybottom = seq_len(1000),
         xright = rep(1, 1000),
         ytop = seq_len(1000) + 1,
         col = heatcols,
         border = NA)
    t1 <- p.max - (p.max - p.min) / 4 * 1
    t2 <- p.max - (p.max - p.min) / 4 * 2
    t3 <- p.max - (p.max - p.min) / 4 * 3
    xy <- par("usr")
    text(xy[2] + strwidth(" "), 1000, round(p.max), xpd = TRUE, adj = 0)
    text(xy[2] + strwidth(" "), 750, round(t1), xpd = TRUE, adj = 0)
    text(xy[2] + strwidth(" "), 500, round(t2), xpd = TRUE, adj = 0)
    text(xy[2] + strwidth(" "), 250, round(t3), xpd = TRUE, adj = 0)
    text(xy[2] + strwidth(" "), 0, round(p.min), xpd = TRUE, adj = 0)
    text(xy[1] - strwidth("  "), 1000, expression("-log"[10]~italic(p)~"Value"),
         xpd = TRUE, cex = 1, adj = 1, srt = 90)

}


# phenos scale function here
pscale <- function(x){
    # remove outlier
    x <- removeOutlier(x)
    # scale
    x.max <- max(x, na.rm = TRUE)
    x.min <- min(x, na.rm = TRUE)
    x <- (x - x.min)/(x.max - x.min)
    return(100 * x)
}

#' @importFrom stats IQR quantile
removeOutlier <- function(x){
    outlier_limup <-
        3 * IQR(x, na.rm = TRUE) +
        quantile(x, 3 / 4, na.rm = TRUE, names = FALSE)# Q3+k(Q3-Q1)

    outlier_limdown <-
        quantile(x, 1 / 4, na.rm = TRUE, names = FALSE) -
        3 * IQR(x , na.rm = TRUE) # Q1-k(Q3-Q1)

    x[x >= outlier_limup | x <= outlier_limdown] = NA
    return(x)
}

