#' @name hapVsPheno
#' @title hapVsPheno
#' @usage
#' hapVsPheno(hap,
#'            pheno,
#'            phenoName, hapPrefix = "H",
#'            title = "",
#'            comparisons = comparisons,
#'            method.args = list(),
#'            symnum.args = list(),
#'            mergeFigs = FALSE,
#'            angle = angle,
#'            minAcc = 5, outlier.rm = TRUE,
#'            method = "t.test", ...)
#' @examples
#'
#' \donttest{
#' data("geneHapR_test")
#' # plot the figs directly
#' hapVsPheno(hap = hapResult,
#'            pheno = pheno,
#'            phenoName = "GrainWeight.2021",
#'            minAcc = 3)
#'
#' # do not merge the files
#' results <- hapVsPheno(hap = hapResult,
#'                       pheno = pheno,
#'                       phenoName = "GrainWeight.2021",
#'                       minAcc = 3,
#'                       mergeFigs = FALSE)
#' plot(results$fig_pvalue)
#' plot(results$fig_Violin)
#' }
#' @param hap object of hapResult class, generate with`vcf2hap()` or
#' `seqs2hap()`
#' @param pheno object of data.frame class, imported by `import_pheno()`
#' @param phenoName pheno name for plot, should be one column name of pheno
#' @param hapPrefix prefix of hapotypes, default as "H"
#' @param title a charater which will used for figure title
#' @param mergeFigs bool type, indicate whether merge the heat map and box
#' plot or not. Default as `FALSE`
#' @param minAcc If observations number of a Hap less than this number will
#' not be compared with others or be ploted. Should not less than 3 due to the
#' t-test will meaninglessly. default as 5
#' @param outlier.rm whether remove ouliers, default as TRUE
#' @param angle the angle of x labels
#' @param comparisons a list contains comparison pairs
#' eg. `list(c("H001", "H002"), c("H001", "H004"))`,
#' or a character vector contains haplotype names for comparison,
#' or "none" indicates do not add comparisons.
# @param method a character string indicating which method to be used for comparing means.
# @param ... options will pass to `ggpubr`
#' @inheritDotParams ggpubr::ggviolin
#' @inheritParams ggpubr::stat_compare_means
#' @importFrom stats na.omit t.test
#' @importFrom rlang .data
#' @export
#' @return list. A list contains a character vector with Haps were applied
#' student test, a mattrix contains p-value of each compare of Haps and a
#' ggplot2 object named as figs if mergeFigs set as `TRUE`, or two ggplot2
#' objects names as fig_pvalue and fig_Violin
hapVsPheno <- function(hap,
                       pheno,
                       phenoName,
                       hapPrefix = "H",
                       title = "",
                       comparisons = comparisons,
                       method.args = list(),
                       symnum.args = list(),
                       mergeFigs = FALSE,
                       angle = angle,
                       minAcc = 5,
                       outlier.rm = TRUE,
                       method = "t.test",
                       ...)
{
    if(! inherits(hap, "hapResult"))
        stop("hap should be object of 'hapResult' class")
    if (missing(phenoName)) {
        warning("phenoName is null, will use the first pheno")
        phenoName <- colnames(pheno)[1]
    }
    if (!(phenoName %in% colnames(pheno))) {
        stop("Could not find ", phenoName, " in colnames of pheno")
    }
    result <- list()
    hap <- hap[stringr::str_starts(hap[, 1], hapPrefix),]
    Accessions <- hap[, colnames(hap) == "Accession"]
    haps <- hap[, 1]
    names(haps) <- Accessions

    pheno$Hap <- haps[row.names(pheno)]
    phenop <- pheno[, c("Hap", phenoName)]

    # remove outliers
    if(outlier.rm)
        phenop[, phenoName] <- removeOutlier(phenop[, phenoName])
    phenop <- na.omit(phenop)
    if (nrow(phenop) == 0)
        stop(
            "After removed NAs, accession with certain Hap have no
    observations of ",
            phenoName,
            ". Please check your pheno file."
        )

    hps <- table(phenop$Hap)

    # filter Haps for plot
    if (max(hps) < minAcc)
        stop("there is no haps to plot (no Haps with observations more than ",
             minAcc,
             ")")

    hps <- hps[hps >= minAcc]

    hpsnm <- names(hps)
    hps <- paste0(names(hps), "(", hps, ")")
    names(hps) <- hpsnm

    # T test
    plotHap <- c()
    my_comparisons <- list()
    T.Result <- matrix(nrow = length(hpsnm), ncol = length(hpsnm))
    colnames(T.Result) <- hpsnm
    row.names(T.Result) <- hpsnm
    nr = nrow(T.Result)
    for (m in seq_len(nr)) {
        for (n in nr:m) {
            i <- hpsnm[m]
            j <- hpsnm[n]
            hapi <- phenop[phenop$Hap == i, phenoName]
            hapj <- phenop[phenop$Hap == j, phenoName]
            if (length(hapi) >= minAcc & length(hapj) >= minAcc) {
                pvalue <- try(t.test(hapi, hapj)$p.value, silent = TRUE)

                T.Result[j, i] <- ifelse(is.numeric(pvalue) &
                                             !is.na(pvalue),
                                         pvalue, 1)
                T.Result[i, j] <- T.Result[j, i]
                plotHap <- c(plotHap, i, j)
                if (T.Result[i, j] < 0.05) {
                    my_comparisons <- c(my_comparisons,
                                        list(hps[c(i, j)]))
                }
            }
        }
    }

    result$plotHap <- unique(plotHap)
    result$T.Result <- T.Result
    plotHap <- unique(plotHap)
    if (is.null(plotHap))
        stop("there is no haps to plot (no Haps with observations more than ",
             minAcc)

    if (length(plotHap) > 1) {
        T.Result <- T.Result[!is.na(T.Result[, 1]), !is.na(T.Result[1, ])]

        # ggplot

        if (nrow(T.Result) > 1)  {
            # get upper or lower tri
            T.Result[lower.tri(T.Result)] = NA
        }
        melResult <- reshape2::melt(T.Result, na.rm = TRUE)

        melResult$label <- ifelse(melResult$value > 1,
                                  1,
                                  ifelse(
                                      melResult$value < 0.001,
                                      0.001,
                                      round(melResult$value, 3)
                                  ))

        fig1 <- ggplot2::ggplot(data = melResult,
                                mapping = ggplot2::aes_(
                                    x =  ~ Var1,
                                    y =  ~ Var2,
                                    fill =  ~ value
                                )) +
            ggplot2::geom_tile(color = "white") +
            ggplot2::ggtitle(label = title, subtitle = phenoName) +
            ggplot2::scale_fill_gradientn(
                colours = c("red", "grey", "grey90"),
                limit = c(0, 1.000001),
                name = parse(text = "italic(p)~value")
            ) +
            ggplot2::geom_text(
                ggplot2::aes_(
                    x =  ~ Var1,
                    y =  ~ Var2,
                    label =  ~ label
                ),
                color = "black",
                size = 4
            ) +
            ggplot2::theme(
                axis.title.x =  ggplot2::element_blank(),
                axis.title.y =  ggplot2::element_blank(),
                panel.grid.major =  ggplot2::element_blank(),
                panel.border =  ggplot2::element_blank(),
                panel.background =  ggplot2::element_blank(),
                axis.ticks =  ggplot2::element_blank(),
                plot.subtitle = ggplot2::element_text(hjust = 0.5),
                plot.title = ggplot2::element_text(hjust = 0.5)
            ) +
            ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top",
                                                           title.hjust = 0.5))
    } else
        fig1 <- ggplot2::ggplot() + ggplot2::theme_minimal()

    # boxplot
    data <- phenop[phenop$Hap %in% plotHap,]
    data <- data[order(data$Hap, decreasing = FALSE),]

    data$Hap <- hps[data$Hap]
    capt <- stringr::str_split(phenoName, "[.]")[[1]][2]
    if (is.na(capt))
        fig2 <- ggpubr::ggviolin(
            data,
            x = "Hap",
            y = phenoName,
            color = "Hap",
            legend = "right",
            legend.title = "",
            add = "boxplot",
            ...
        )
    else
        fig2 <- ggpubr::ggviolin(
            data,
            x = "Hap",
            y = phenoName,
            color = "Hap",
            caption = capt,
            legend = "right",
            legend.title = "",
            add = "boxplot",
            ...
        )

    if(missing(angle))
        angle <- ifelse(length(hps) >= 6, 45, 0)
    fig2 <- fig2 +  # do not modify here
        #    stat_compare_means(label.y = max(data[,2]))+
        #    no comparision by remove this line (Kruskal-Wallis test)
        ggplot2::ggtitle(label = title) +
        ggplot2::theme(
            plot.subtitle = ggplot2::element_text(hjust = 0.5),
            axis.text.x = ggplot2::element_text(
                angle = angle,
                hjust = ifelse(length(hps) >= 6, 1, 0.5)
            ),
            plot.title = ggplot2::element_text(hjust = 0.5)
        ) +
        ggplot2::ylab(stringr::str_split(phenoName, "[.]")[[1]][1])

    if(! missing(comparisons)){
        if(inherits(comparisons, "list")) {
            for(i in seq_len(length(comparisons)))
                comparisons[[i]] <- hps[comparisons[[i]]]
            my_comparisons <- comparisons
        } else if(inherits(comparisons, "character")){
            if(comparisons[1] == "none") {
                my_comparisons <- list()
            } else {
                comparisons <- hps[comparisons]
                probe <- c()
                for(i in my_comparisons)
                    probe <- c(probe, TRUE %in% (i %in% comparisons))
                my_comparisons <- my_comparisons[probe]
            }
        }
    }

    if (length(my_comparisons) > 0) {
        fig2 <- fig2 + # 添加箱线图
            ggpubr::stat_compare_means(
                comparisons = unique(my_comparisons),
                paired = FALSE,
                method.args = method.args,
                symnum.args = symnum.args,
                method = method
            )
    }

    if (mergeFigs)  {
        fig3 <- ggpubr::ggarrange(fig1,
                                  fig2,
                                  nrow = 1,
                                  labels = c("A", "B"))
        result$figs <- fig3
    } else {
        result$fig_pvalue <- fig1
        result$fig_Violin <- fig2
    }

    return(result)
}



#' @name hapVsPhenos
#' @title hapVsPhenos
#' @usage
#' hapVsPhenos(hap, pheno,
#'          outPutSingleFile = TRUE,
#'          hapPrefix = "H",
#'          title = "Seita.0G000000",
#'          width = 12,
#'          height = 8,
#'          res = res,
#'          compression = "lzw",
#'          filename.prefix = filename.prefix,
#'          filename.surfix = "pdf",
#'          filename.sep = "_",
#'          outlier.rm = TRUE,
#'          ...)
#' @param outPutSingleFile `TRUE` or `FALSE` indicate whether put all figs
#' into to each pages of single file or generate multi-files.
#' Only worked while file type is pdf
#' @param width manual option for determining the output file width in inches.
#' (default: 12)
#' @param height manual option for determining the output file height in inches.
#' (default: 8)
#' @param res The nominal resolution in ppi which will be recorded in the
#' bitmap file, if a positive integer. Also used for units other than the
#' default, and to convert points to pixels
#' @inheritParams grDevices::tiff
#' @param filename.prefix,filename.surfix,filename.sep
#' if multi files generated, file names will be formed by
#' prefix `filename.prefix`, a seperate charcter `filename.sep`,
#' pheno name, a dot and surfix `filename.surfix`,
#' and file type was decide by `filename.surfix`;
#' if single file was generated, file name will be formed by
#' prefix `filename.prefix`, a dot and surfix `filename.surfix`
#' @inheritParams hapVsPheno
#' @inheritDotParams hapVsPheno
#' @examples
#' data("geneHapR_test")
#'
#' oriDir <- getwd()
#' setwd(tempdir())
#' # analysis all pheno in the data.frame of pheno
#' hapVsPhenos(hapResult,
#'             pheno,
#'             outPutSingleFile = TRUE,
#'             hapPrefix = "H",
#'             title = "Seita.0G000000",
#'             filename.prefix = "test",
#'             width = 12,
#'             height = 8,
#'             res = 300)
#' setwd(oriDir)
#' @importFrom stats na.omit t.test
#' @import grDevices
#' @export
#' @return No return value
hapVsPhenos <- function(hap,
                        pheno,
                        outPutSingleFile = TRUE,
                        hapPrefix = "H",
                        title = "Seita.0G000000",
                        width = 12,
                        height = 8,
                        res = res,
                        compression = "lzw",
                        filename.prefix = filename.prefix,
                        filename.surfix = "pdf",
                        filename.sep = "_",
                        outlier.rm = TRUE,
                        ...) {

    # pheno association
    if (missing(hap))
        stop("hap is missing!")

    if (missing(pheno))
        stop("pheno is missing!")

    if(missing(filename.prefix))
        stop("filename.prefix is missing!")

    if (!filename.surfix %in% c("pdf", "png", "tif", "tiff", "jpg", "jpeg", "bmp"))
        stop("The file type should be one of pdf, png, tiff, jpg and bmp")

    if (filename.surfix != "pdf")
        outPutSingleFile <- FALSE

    probe <- ifelse(filename.surfix == "pdf",
                    ifelse(outPutSingleFile,
                           TRUE,
                           FALSE),
                    FALSE)
    if (probe) {
        filename <- paste0(filename.prefix,".",filename.surfix)
        message("
        File type is pdf and 'outPutSingleFile' set as TRUE,
        all figs will plot in ",
                filename)
        pdf(filename, width = width, height = height)
        on.exit(dev.off())
    }
    if (!is.data.frame(pheno))
        stop("pheno should be a data.frame object")
    if (ncol(pheno) == 1)
        warning("There is only one col detected in pheno, 'hapVsPheno' is prefered")
    phenoNames <- colnames(pheno)
    steps <- 0
    for (phenoName in phenoNames) {
        if (!probe) {
            filename <- paste0(filename.prefix,
                               filename.sep,
                               phenoName,
                               ".",
                               filename.surfix)
            switch(
                filename.surfix,
                "pdf" = pdf(filename, width = width, height = height),
                "png" = png(
                    filename = filename,
                    width = width,
                    height = height,
                    units = "in",
                    res = res
                ),
                "bmp" = bmp(
                    filename = filename,
                    width = width,
                    height = height,
                    units = "in",
                    res = res
                ),
                "jpg" = png(
                    filename = filename,
                    width = width,
                    height = height,
                    units = "in",
                    res = res
                ),
                "jpeg" = png(
                    filename = filename,
                    width = width,
                    height = height,
                    units = "in",
                    res = res
                ),
                "tif" = tiff(
                    filename = filename,
                    width = width,
                    height = height,
                    units = "in",
                    res = res,
                    compression = compression
                ),
                "tiff" = tiff(
                    filename = filename,
                    width = width,
                    height = height,
                    units = "in",
                    res = res,
                    compression = compression
                )
            )
        }
        steps <- steps + 1
        message("Total: ",
                ncol(pheno),
                "; current: ",
                steps,
                ";\tphynotype: ",
                phenoName, appendLF = FALSE)
        cat("\tfile: ", filename, "\n", sep = "")

        resulti <- try(hapVsPheno(hap = hap,
                                  pheno = pheno,
                                  phenoName = phenoName,
                                  hapPrefix = hapPrefix,
                                  title = title,
                                  mergeFigs = TRUE,
                                  ...))
        if(!inherits(resulti, "try-error")) plot(resulti$figs) else resulti
        if (!probe)
            dev.off()
        resulti <- NULL
    }
}


# getSurFix <- function(Name) {
#     parts <- strsplit(Name, ".", fixed = TRUE,)
#     lparth <- length(parts[[1]])
#     surFix <- parts[[1]][lparth]
#     surFix <- stringr::str_to_lower(surFix)
#     return(surFix)
# }
