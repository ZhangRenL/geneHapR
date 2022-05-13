#' @name hapVsPheno
#' @title hapVsPheno
#' @usage hapVsPheno(hap,
#' pheno,
#' phenoName,hapPrefix = "H",
#' geneID = "test1G0387",
#' mergeFigs = TRUE,
#' minAcc = 5)
#' @examples
#'
#' data("quickHap_test")
#' hap <- get_hap(vcf,hyb_remove = TRUE, na.drop = TRUE)
#' # plot the figs direactively
#' hapVsPheno(hap = hap,pheno = pheno,phenoName = "GrainWeight.2021",minAcc = 3)
#'
#' #do not merge the files
#' results <- hapVsPheno(hap = hap,
#'                       pheno = pheno,
#'                       phenoName = "GrainWeight.2021",
#'                       minAcc = 3,
#'                       mergeFigs = FALSE)
#' plot(results$fig_pvalue)
#' plot(results$fig_Violin)
#' # plot multi figs in a `for` work folw
#' pheno$GrainWeight.2022 = pheno$GrainWeight.2021 + c(1:nrow(pheno))
#' for(i in colnames(pheno)){
#'     results <- hapVsPheno(hap = hap,pheno = pheno,phenoName = i,minAcc = 3)
#'     png(paste0(i,".png"))
#'     plot(results$figs)
#' dev.off()
#' }
#' @param hap hap
#' @param pheno pheno
#' @param phenoName pheno name
#' @param hapPrefix hap  prefix
#' @param geneID gene ID for plot title
#' @param mergeFigs merge heatmap and box plot or not
#' @param minAcc If observations numberof a Hap less than this number will
#' not be compared with others or ploted
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
    minAcc = 5)
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
            ggplot2::ggtitle(label = geneID,subtitle = phenoName)+
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
        add = "boxplot")  +  # 不要动
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
    fig3 <- ggpubr::ggarrange(fig1, fig2, nrow = 1, labels = c("A","B"))
    if(mergeFigs)  {
        result$figs <- fig3
    } else {
        result$fig_pvalue <- fig1
        result$fig_Violin <- fig2
    }

    return(result)
}


# @name hapVsPhenos
# @title hapVsPhenoS
# @usage hapVsPhenos(hap, phenos,hapPrefix = "H",geneID = "Seita.0G000000",
# mergeFigs = T)
# @param hap hap
# @param phenos phenos
# @param hapPrefix hap  prefix
# @param  geneID gene ID forplot title
# @param mergeFigs merge heatmapand box plot ornot
# @importFrom stats na.omit t.test
# @export
#hapVsPhenos <- function(hap, phenos,hapPrefix = "H",geneID = "Seita.0G000000",
# mergeFigs = T){
#  # 表型关联
#  phenoNames <- colnames(phenos)
#
#   Accs <- row.names(phenos)
#   results = list()
#   for (phenoName in phenoNames){
#     pheno <- phenos[,phenoName]
#     pheno <- as.data.frame(pheno,row.names = Accs)
#     colnames(pheno) <- phenoName
#     resulti <- hapVsPheno(hap, pheno, phenoName =  phenoName,
#     hapPrefix=hapPrefix, geneID = geneID, mergeFigs = mergeFigs)
#     results = c(results, resulti)
#   }
#
#   return(results)
# }
