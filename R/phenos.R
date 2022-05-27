#' @name hapVsPheno
#' @title hapVsPheno
#' @usage hapVsPheno(hap,
#' pheno,
#' phenoName,hapPrefix = "H",
#' geneID = "test1G0387",
#' mergeFigs = TRUE,
#' minAcc = 5,
#' ...)
#' @examples
#'
#' data("quickHap_test")
#' hap <- get_hap(vcf,hyb_remove = TRUE, na.drop = TRUE)
#' # plot the figs direactively
#' hapVsPheno(hap = hap,pheno = phenos,phenoName = "GrainWeight.2021",minAcc = 3)
#'
#' #do not merge the files
#' results <- hapVsPheno(hap = hap,
#'                       pheno = phenos,
#'                       phenoName = "GrainWeight.2021",
#'                       minAcc = 3,
#'                       mergeFigs = FALSE)
#' plot(results$fig_pvalue)
#' plot(results$fig_Violin)
#'
#' \dontrun{
#' # plot multi figs in a 'for' work folw
#' pheno$GrainWeight.2022 = pheno$GrainWeight.2021 + c(1:nrow(pheno))
#' for(i in colnames(pheno)){
#'     results <- hapVsPheno(hap = hap,pheno = phenos,phenoName = i,minAcc = 3)
#'     png(paste0(i,".png"))
#'     plot(results$figs)
#' dev.off()
#' }
#' }
#' @param hap hap
#' @param pheno pheno
#' @param phenoName pheno name
#' @param hapPrefix hap  prefix
#' @param geneID gene ID for plot title
#' @param mergeFigs merge heatmap and box plot or not
#' @param minAcc If observations numberof a Hap less than this number will
#' not be compared with others or ploted
#' @param ... options will pass to ggpubr
#' @importFrom stats na.omit t.test
#' @importFrom rlang .data
#' @export
#' @return list. A list contains a character vector with Haps were applied
#' student test, a mattrix contains p-value of each compare of Haps and a
#' ggplot2 object named as figs if mergeFigs set as TRUE, or two ggplot2
#' objects names as fig_pvalue and fig_Violin
hapVsPheno <- function(
    hap,
    pheno,
    phenoName,
    hapPrefix = "H",
    geneID = "test1G0387",
    mergeFigs = TRUE,
    minAcc = 5,
    ...)
{
    if(missing(phenoName)) {
        warning("phenoName is null, will use the first pheno")
        phenoName <- colnames(pheno)[1]
    }
    if(!(phenoName %in% colnames(pheno))) {
        stop("Could not find ", phenoName, " in colnames of pheno")
    }
    result <- list()
    hap <- hap[stringr::str_starts(hap[,1],hapPrefix),]
    Accessions <- hap[,colnames(hap) == "Accession"]
    haps <- hap[,1]
    names(haps) <- Accessions

    pheno$Hap <- haps[row.names(pheno)]
    phenop <- pheno[,c("Hap",phenoName)]
    phenop <- na.omit(phenop)
    if(nrow(phenop) == 0) stop(
    "After removed NAs, accession with certain Hap have no
    observations of ", phenoName,". Please check your pheno file.")

    hps <- table(phenop$Hap)

    # filter Haps for plot
    if(max(hps) < minAcc) stop(
        "there is no haps to plot (no Haps with observations more than ",
        minAcc, ")")

    hps <- hps[hps >= minAcc]

    hpsnm <- names(hps)
    hps <- paste0(names(hps),"(",hps,")")
    names(hps) <- hpsnm

    # T 检验
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
                pvalue <- try(t.test(hapi,hapj)$p.value,silent = TRUE)

                T.Result[j, i] <- ifelse(
                    is.numeric(pvalue) & !is.na(pvalue),
                    pvalue, 1)
                T.Result[i, j] <- T.Result[j, i]
                plotHap <- c(plotHap, i, j)
                if(T.Result[i, j] < 0.05) {
                    my_comparisons <- c(
                        my_comparisons,
                        list(hps[c(i,j)]))
                }
            }
        }
    }

    result$plotHap <- unique(plotHap)
    result$T.Result <- T.Result
    plotHap <- unique(plotHap)
    if(is.null(plotHap)) stop(
        "there is no haps to plot (no Haps with observations more than ",
        minAcc)

    if(length(plotHap) > 1){
        T.Result <- T.Result[!is.na(T.Result[, 1]), !is.na(T.Result[1, ])]

        # ggplot

        if(nrow(T.Result) > 1)  { # 获得矩阵的上三角或下三角
            T.Result[lower.tri(T.Result)] = NA
        }
        melResult <- reshape2::melt(T.Result, na.rm = TRUE)

        melResult$label <- ifelse(
            melResult$value>1,
            1,
            ifelse(
                melResult$value<0.001,
                0.001,
                round(melResult$value,3)))

        fig1 <- ggplot2::ggplot(
            data = melResult,
            mapping = ggplot2::aes_(x =~Var1, y =~Var2, fill =~value)) +
            ggplot2::geom_tile(color = "white") +
            ggplot2::ggtitle(label = geneID, subtitle = phenoName)+
            ggplot2::scale_fill_gradientn(
                colours = c("red","grey", "grey90"),
                limit = c(0, 1.000001),
                name = parse(text = "italic(p)~value")) +
            ggplot2::geom_text(
                ggplot2::aes_(
                    x=~Var1, y=~Var2,
                    label =~ label),
                color = "black",
                size = 4) +
            ggplot2::theme(
                axis.title.x =  ggplot2::element_blank(),
                axis.title.y =  ggplot2::element_blank(),
                panel.grid.major =  ggplot2::element_blank(),
                panel.border =  ggplot2::element_blank(),
                panel.background =  ggplot2::element_blank(),
                axis.ticks =  ggplot2::element_blank(),
                plot.subtitle = ggplot2::element_text(hjust = 0.5),
                plot.title = ggplot2::element_text(hjust = 0.5)) +
            ggplot2::guides(
                fill = ggplot2::guide_colorbar(
                    title.position = "top",
                    title.hjust = 0.5))
    } else fig1 <- ggplot2::ggplot() + ggplot2::theme_minimal()

    # 作图，箱线图
    data <- phenop[phenop$Hap %in% plotHap,]
    data <- data[order(data$Hap,decreasing = FALSE),]

    data$Hap <- hps[data$Hap]

    fig2 <- ggpubr::ggviolin(
        data,
        x = "Hap",
        y = phenoName,
        color = "Hap",
        caption = stringr::str_split(phenoName,"[.]")[[1]][2],
        legend = "right", legend.title = "",
        add = "boxplot", ...)  +  # 不要动
        #    stat_compare_means(label.y = max(data[,2]))+
        # 去掉这行就没有比较了(Kruskal-Wallis test)
        ggplot2::ggtitle(label = geneID) +
        ggplot2::theme(
            plot.subtitle = ggplot2::element_text(
                hjust = 0.5),
            axis.text.x = ggplot2::element_text(
                angle=ifelse(length(hps) >= 6, 45, 0),
                hjust = ifelse(length(hps) >= 6, 1, 0.5)),
            plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::ylab(stringr::str_split(phenoName,"[.]")[[1]][1])
    if(length(my_comparisons) > 0) {
        fig2 <- fig2 + # 添加箱线图
            ggpubr::stat_compare_means(
                comparisons = unique(my_comparisons),
                paired = FALSE, method = "t.test")
    }

    if(mergeFigs)  {
        fig3 <- ggpubr::ggarrange(fig1, fig2, nrow = 1, labels = c("A","B"))
        result$figs <- fig3
    } else {
        result$fig_pvalue <- fig1
        result$fig_Violin <- fig2
    }

    return(result)
}



#' @name hapVsPhenos
#' @title hapVsPhenoS
#' @usage
#' hapVsPhenos(hap,
#'             phenos,
#'             outPutSingleFile = TRUE,
#'             hapPrefix = "H",
#'             geneID = "Seita.0G000000",
#'             mergeFigs = TRUE,
#'             file = file,
#'             width = 12,
#'             height = 8, ...)
#' @param hap hap
#' @param phenos phenos
#' @param outPutSingleFile only worked while file type is pdf.
#' Out all figs in a single pdf file, figs will plot in each pages.
#' @param hapPrefix hap  prefix
#' @param geneID gene ID forplot title
#' @param mergeFigs merge heatmapand box plot or not
#' @param file out put file name. File type could be pdf, png, tiff, jpg, bmp.
#' @param width manual option for determining the output file width in inches.
#' (defalt: 12)
#' @param height manual option for determining the output file height in inches.
#' (defalt: 8)
#' @param ... options will pass to ggpubr
#' @importFrom stats na.omit t.test
#' @import grDevices
#' @export
#' @return NULL
hapVsPhenos <- function(hap,
                        phenos,
                        outPutSingleFile = TRUE,
                        hapPrefix = "H",
                        geneID = "Seita.0G000000",
                        mergeFigs = TRUE,
                        file = file,
                        width = 12,
                        height = 8,
                        ...){
    # 表型关联
    if(missing(hap)) stop("hap is missing!")
    if(missing(phenos)) stop("phenos is missing!")
    if(missing(file)) stop("Error: Output file name/path is missing!")

    surFix <- getSurFix(file)
    if(!surFix %in% c("pdf", "png", "tif", "tiff", "jpg", "jpeg", "bmp"))
        stop("The file type should be one of pdf, png, tiff, jpg and bmp")
    if(surFix != "pdf") outPutSingleFile <- FALSE

    probe <- ifelse(surFix == "pdf",
                    ifelse(outPutSingleFile,
                           TRUE,
                           FALSE),
                    FALSE)
    if(probe){
        message("File type is pdf and 'outPutSingleFile' set as TRUE, all figs will plot in ", file)
        pdf(file, width = width, height = height)
    }
    if(!is.data.frame(phenos)) stop("phenos should be a data.frame object")
    if(ncol(phenos) == 1)
        warning("There is only one col detected in phenos, 'hapVsPheno' is prefered")
    phenoNames <- colnames(phenos)
    for (phenoName in phenoNames){
        if(!probe){
            switch(
                surFix,
                "pdf" = pdf(file, width = width, height = height),
                "png" = png(filename = file,
                            width = width, height = height, units = "in", res = 300),
                "bmp" = bmp(filename = file,
                            width = width, height = height, units = "in", res = 300),
                "jpg" = png(filename = file,
                            width = width, height = height, units = "in", res = 300),
                "jpeg" = png(filename = file,
                             width = width, height = height, units = "in", res = 300),
                "tif" = tiff(filename = file,
                             width = width, height = height, units = "in", res = 300),
                "tiff" = tiff(filename = file,
                              width = width, height = height, units = "in", res = 300))
        }
        resulti <- hapVsPheno(hap,
                              pheno = phenos,
                              phenoName = phenoName,
                              hapPrefix=hapPrefix,
                              geneID = geneID,
                              mergeFigs = mergeFigs,
                              ...)

        if(mergeFigs){
            plot(resulti$figs)
        } else {
            plot(resulti$fig_pvalue)
            plot(resulti$fig_Violin)
        }
        if(!probe) dev.off()
    }
    if(probe) dev.off()
}


getSurFix <- function(Name){
    parts <- strsplit(Name, ".", fixed = TRUE,)
    lparth <- length(parts[[1]])
    surFix <- parts[[1]][lparth]
    surFix <- stringr::str_to_lower(surFix)
    return(surFix)
}
