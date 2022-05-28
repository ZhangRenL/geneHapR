#' @name plotHapTable
#' @title plotHapTable
#' @usage
#' plotHapTable(hapResult,
#'              hapPrefix = "H",
#'              geneID = "",
#'              displayIndelSize = 0, angle = c(0,45,90),
#'              replaceMultiAllele = TRUE,
#'              title.color = "grey90")
#' @description plot you hap result as a table
#' @examples
#'
#' data("quickHap_test")
#' hap <- get_hap(vcf,hyb_remove = TRUE, na.drop = TRUE)
#' hapResult <- hap_result(hap)
#' plotHapTable(hapResult)
#' @importFrom randomcoloR randomColor
#' @importFrom stringr str_starts
#' @importFrom stringr str_length
#' @import tidyr
#' @import ggplot2
#' @param hapResult hapResult
#' @param title.color title.color, defalt is grey90
#' @param hapPrefix hapPrefix
#' @param geneID geneID used as mainTitle
#' @param displayIndelSize display indels with max size of displayIndelSize,
#' If set as 0, all indels will convert into 'i*'.
#' @param angle angle of coordinate, should be one of 0, 45 and 90
#' @param replaceMultiAllele Whether to replace MultiAllele with T*,
#' defalt as TRUE
#' @export
#' @return ggplot2 object
plotHapTable <- function(
    hapResult,
    hapPrefix = "H",
    geneID = "",
    displayIndelSize = 0, angle = c(0,45,90), replaceMultiAllele = TRUE,
    title.color = "grey90")
{
    requireNamespace('tidyr')
    if("Accession" %in% colnames(hapResult)) {
        hapResult <- hapResult[,colnames(hapResult) != 'Accession']
    }

    hps <- hapResult[stringr::str_starts(hapResult[,1],hapPrefix),] %>%
        as.matrix()
    if(nrow(hps) <= 1)
        stop("please check 'hapResult' and 'hapPrefix'")

    # foot and labs
    ALLELE <- hapResult[hapResult[,1] == "ALLELE",]
    probe_indel <- is.indel.allele(ALLELE)
    probe_mula <- is.multiallelic.allele(ALLELE)

    footi <- ""
    if(TRUE %in% probe_indel){
        # set notes
        displayIndelSize <- displayIndelSize + 1
        allIndel <- ALLELE[probe_indel]
        allIndel <- unlist(stringr::str_split(allIndel,"[,/]"))
        allIndel <- allIndel[stringr::str_length(allIndel) > displayIndelSize] %>%
            unique()
        notes <- paste0("i",seq_len(length(allIndel)))
        names(notes) <- allIndel
        m <- paste(names(notes),"->", notes, collapse = "; ",sep = "")
        message("Indel replcements are:\n",m)
        # replace Indel by notes
        hps[hps %in% allIndel] <- notes[hps[hps %in% allIndel]]

        for(i in seq_len(length(ALLELE))){
            if(probe_indel[i]){
                ALi <- ALLELE[i]
                ALi <- unlist(stringr::str_split(ALi, "[,/]"))
                p <- ALi %in% names(notes)
                ALi[p] <- notes[ALi[p]]
                ALi <- paste(ALi, collapse = ",")
                ALi <- stringr::str_replace(ALi,",","/")
                ALLELE[i] <- ALi
            }
        }
        # set footi
        footi <- paste(notes, names(notes),sep = ":", collapse = "; ")
    }

    # replace multiallele title
    footT <- ""
    if(replaceMultiAllele & TRUE %in% probe_mula){
        rept <- ALLELE[probe_mula]
        noteT <- paste0("T",seq_len(length(rept)))
        names(noteT) <- ALLELE[probe_mula]
        ALLELE[probe_mula] <- noteT
        footT <- paste(noteT,names(noteT), sep = ":", collapse = "; ")
    }

    if(nchar(footi) > 0 & nchar(footT) > 0)
        foot <- paste(footT, footi, sep = "\n") else
            foot <- paste0(footT, footi)

    # set labs for plot
    hps <- rbind(ALLELE, hps)
    meltHapRes <- reshape2::melt(hps, 1)
    colnames(meltHapRes) <- c('Var1','Var2',"value")
    lab <- meltHapRes
    meltHapRes$value[meltHapRes$Var1 == "ALLELE"] <- NA
    meltHapRes$value[meltHapRes$Var2 == "freq"] = NA


    levels <- as.vector(unique(meltHapRes$Var1))
    levels <- levels[order(levels, decreasing = TRUE)]
    meltHapRes$Var1 <- factor(meltHapRes$Var1, levels = levels)
    if(length(angle) > 1) warning("using first 'angle': ",angle[1])
    angle <- angle[1]
    if(angle == 0) {vjust <- 0.5; hjust  <- 0.5} else
        if(angle == 45) {vjust <- 0.1; hjust  <- 0.1} else
            if(angle == 90) {vjust <- 1; hjust  <- 1} else
                stop("angle should be one of 0, 45 and 90")
    fig0 <- ggplot2::ggplot(
        data = meltHapRes,
        mapping = ggplot2::aes_(x=~Var2,
                                y=~Var1,
                                fill=~value)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::geom_text(
            ggplot2::aes_(
                x=~Var2,
                y=~Var1,
                label = lab$value)) +
        ggplot2::scale_fill_discrete(na.value = title.color) +
        ggplot2::labs(caption = foot) +
        ggplot2::ggtitle(label = geneID) +  ggplot2::scale_y_discrete() +
        ggplot2::scale_x_discrete(
            guide = ggplot2::guide_axis(position = "top")) +

        ggplot2::theme(
            legend.position = "none",
            axis.title.x =  ggplot2::element_blank(),
            axis.title.y =  ggplot2::element_blank(),
            panel.grid.major =  ggplot2::element_blank(),
            panel.border =  ggplot2::element_blank(),
            panel.background =  ggplot2::element_blank(),
            axis.ticks =  ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = angle,
                                                vjust = vjust,
                                                hjust = hjust),
            plot.subtitle = ggplot2::element_text(hjust = 0.5),
            plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::guides(
            fill = ggplot2::guide_colorbar(
                title.position = "top",
                title.hjust = 0.5))
    return(fig0)
}






#' @name plotGeneStructure
#' @title plotGeneStructure
#' @usage plotGeneStructure(gff, hapResult,
#' Chr,
#' startPOS, endPOS,
#' type = "pin", cex = 0.7,
#' CDS_h = 0.05, fiveUTR_h = 0.02, threeUTR_h = 0.01)
#' @examples
#'
#' data("quickHap_test")
#' hap <- get_hap(vcf,hyb_remove = TRUE, na.drop = TRUE)
#' hapResult <- hap_result(hap)
#' plotGeneStructure(gff, hapResult,
#'                   startPOS = 4100,
#'                   endPOS = 8210,
#'                   cex = 0.75)
#' @importFrom trackViewer lolliplot
#' @importFrom  GenomicRanges GRanges
#' @importFrom  GenomicRanges strand
#' @importFrom IRanges IRanges `%over%`
#' @import tidyr
#' @param gff gff
#' @param hapResult ouput of hap_result(), a data.frame consistant
#' with 4 metarows: Chr, POS, Allele and ANN, and each genotype of each hap
#' @param Chr Chr, it will use the first Chr in the meta rows by defalt
#' @param startPOS If missing, will use the min position
#'  in the second meta rows by defalt
#' @param endPOS If missing, will use the max position
#'  in the second meta rows by defalt
#' @param cex cex will control the size of circle.
#' @param type character. Could be circle, pie, pin, pie.stack or flag
#' @param CDS_h,fiveUTR_h,threeUTR_h The height of CDS 5'UTR and 3'UTR
#' @export
#' @return NULL
plotGeneStructure <- function(
    gff, hapResult,
    Chr, startPOS, endPOS,
    type = "pin", cex = 0.7,
    CDS_h = 0.05,
    fiveUTR_h = 0.02,
    threeUTR_h = 0.01)
{
    # lolliplot
    requireNamespace("trackViewer")
    requireNamespace("tidyr")
    if(missing(gff)) {
        stop("gff is missing!")}
    if(missing('hapResult')) {
        stop("hapResult is missing!")}

    meta <- hapResult[seq_len(4),-1]
    if("Accession" %in% colnames(meta)) {
        meta <- meta[,colnames(meta) != "Accession"]
    }
    if("freq" %in% colnames(meta)) {
        meta <- meta[,colnames(meta) != "freq"]
    }

    POS = meta[2,] %>% as.numeric()
    if(missing(Chr)) Chr <- meta[1,1]
    if(missing(startPOS)) startPOS <- min(POS)
    if(missing(endPOS)) endPOS <- max(POS)
    SNP <- meta[4,]

    #meta <- hapResult[1:4,]
    #meta <- meta[,-ncol(meta)]
    #meta[meta == ""] = NA
    #meta <- meta[,!is.na(meta[1,])]
    #POS <- as.numeric(meta[2,])


    geneElement <- c("CDS","three_prime_UTR","five_prime_UTR")

    SNP.gr <- GenomicRanges::GRanges(
        Chr, IRanges::IRanges(
            POS, width=1,
            names = paste0(POS,"(",SNP,")")),
        color = sample.int(6, length(SNP), replace=TRUE),
        score = sample.int(5, length(SNP), replace = TRUE),
        angle = 45)


    # set plot ranges
    gene <- GenomicRanges::GRanges(
        Chr,
        IRanges::IRanges(
            start = min(startPOS,endPOS),
            end = max(startPOS,endPOS)))

    over <- gff[gff %over% gene]
    over$height[over$type == "CDS"] <- CDS_h
    over$height[over$type == "three_prime_UTR"] <- threeUTR_h
    over$height[over$type == "five_prime_UTR"] <- fiveUTR_h

    features <- over[over$type %in% geneElement]
    strands <- as.character(GenomicRanges::strand(features))
    layerID <- unlist(features$Parent)
    layerID <- stringr::str_remove_all(layerID,".v2.2")
    layerID <- paste0(
        layerID, "(",
        ifelse(strands == "+", "5'->3'","3'<-5'"), ")")
    features$featureLayerID <- layerID
    names(features) <- features$featureLayerID
    l <- length(unique(names(features)))
    if (l < 6) {
        fillc <- c(seq_len(l) + 1)
    }else{
        fillc <- randomcoloR::randomColor(l)
    }
    names(fillc) <- unique(names(features))
    features$fill <- fillc[names(features)]


    # set ranges of features
    trackViewer::lolliplot(
        SNP.gr, features, gene,
        type = type, jitter = NULL,
        cex = cex,
        ylab = "", yaxis = FALSE)
}

