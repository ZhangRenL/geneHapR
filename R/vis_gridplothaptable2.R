#' @title plotHapTable2
#' @description
#' plot the hapResult in table like style using grid system.
#' @examples
#' #
#' data(geneHapR_test)
#' plotHapTable2(hapResult)
#' plotHapTable2(hapResult, gff = gff)
#' @param hapSummary the hapSummary or hapResult object
#' @param gff gff or bed annotation
#' @param title title of plot
#' @param CDS_height a numeric vector specified the height of CDS, and the height of utr is half of that,
#'  only useful when gff is provided,
#' @param equal_col_width a bool or numeric vector specified whether column with should be equal
#' @param annot_fill,head_fill,labels_fill the fill color of annotation, head and label row or columns
#' @param footbar the foot notes
#' @param gene_model_fill fill color of CDS ad UTR
#' @param cell_fill a color vector or function or named vector specified cell fill color
#' @param replaceMultiAllele replace multi-allele title by 'T1,T2,...' or not
#' @param col_annots,head_anno,headrows,row_annots,row_labels the column or row number of annotation or labels or heads
#' @param show_INFO show annotation field or not, default as `FALSE`
#' @param show_chr_name show chromosome name at left-top cell or not
#' @param show_indel_size the Indel length longer will be replaced by "i1,i2,i3,..."
#' @param Chr,start,end which range should be plotted in gene model
#' @param gene_model_height,table_height,space_height the plotting range height of gene model, table and spacer
#' @export
plotHapTable2 <- function(hapSummary, show_indel_size = 1,replaceMultiAllele = TRUE,
                          show_INFO = FALSE,title = "",
                          gff = NULL,show_chr_name = TRUE,
                          Chr = NULL, start = NULL, end = NULL,gene_model_fill = "white",
                          gene_model_height = 0.2, space_height = 0.1, table_height = NULL, CDS_height = 0.1,
                          headrows = 1, equal_col_width = FALSE,
                          head_anno = 1, col_annots = 0,
                          row_labels = 1, row_annots = 1,
                          labels_fill = "white", annot_fill = "grey90", head_fill = NULL,
                          cell_fill = NULL, style = gpar(fontfamily = "sans", fontface = 1, cex = 0.7),
                          footbar = ""){
    # function start
    if(is.null(cell_fill)) cell_fill <- OPTS$pie.colors.function
    if (inherits(hapSummary, "hapResult")) hapSummary <- hap_summary(hapSummary)

    {
        # 将hapSummary转为普通表格
        table <- data.frame(hapSummary)
        table <- table[,names(table)!= "Accession"]
        table[2,1] <- if(show_chr_name) table[1,2] else ""
        table[2,ncol(table)] = "frequency"
        table = table[-1,]
    }

    {
        # 处理Indel和多等位位点
        if(! show_indel_size ) show_indel_size <- 1
        ALLELE <- as.vector(as.matrix(table[table$Hap == "ALLELE", ]))
        probe_indel <- is.indel.allele(ALLELE)
        probe_mula <- is.multiallelic.allele(ALLELE)
        als <- unique(as.vector(as.matrix(hapSummary)[-c(1:4), seq_len(sites(hapSummary)) + 1]))
        als <- als[nchar(als) > show_indel_size]
        ft_i <- c()
        for (i in seq_len(length(als))){
            al_i <- als[i]
            table[table == al_i] <- paste0("i",i )
            ALLELE <- stringr::str_replace(ALLELE, al_i, paste0("i",i ))
            ft_i <- c(ft_i, paste0("i",i,":",al_i))
        }
        ft_i <- paste(ft_i, collapse = ";")
        if(replaceMultiAllele){
            ft_t <- paste("T",seq_len(sum(probe_mula)),":",ALLELE[probe_mula], sep = "", collapse = ";")
            ALLELE[probe_mula] <- paste0("T",seq_len(sum(probe_mula)))
        }
        table[table$Hap == "ALLELE",] = ALLELE

    }

    {
        # 生成填充表
        table[is.na(table)] <- ""
        if(! show_INFO) {
            table <- table[table[,1] != "INFO", ]
            head_anno <- 1
        } else head_anno <- 2
        nr <- nrow(table)
        nc <- ncol(table)

        # 生成填充表
        fill_table <- table
        fill_table[! is.na(fill_table)] <- NA
        head_fill <- if (is.null(head_fill)) labels_fill else head_fill
        if (col_annots >= 1) fill_table[(nr - seq_len(col_annots) + 1),] <- annot_fill
        if (row_annots >= 1) fill_table[,(nc - seq_len(row_annots) + 1)] <- annot_fill
        if (head_anno >= 1) fill_table[seq_len(headrows + head_anno),] <- annot_fill
        if (row_labels >= 1) fill_table[,seq_len(row_labels)] <- labels_fill
        if (headrows >= 1)   fill_table[seq_len(row_labels),] <- head_fill

        # 处理cell 填充
        rs <- seq(headrows + head_anno + 1, nr - col_annots)
        cs <- seq(row_labels + 1, nc - row_annots)
        probe <- is.na(fill_table)
        unique_cell <- unique(table[probe])
    }

    {
        # cell fill 是命名颜色向量
        if (! is.null(names(cell_fill))){
            if (! all(unique_cell %in% names(cell_fill))) {
                ecp <- unique_cell[! unique_cell %in% names(cell_fill)]
                stop(ecp, "not in names of 'cell_fill'")
            }
        }

        # cell_fill 是function
        if (is.function(cell_fill)){
            cell_fill <- cell_fill(length(unique_cell))
            names(cell_fill) <- unique_cell
        } else {
            if (length(unique_cell) > length(cell_fill))
                cell_fill <- rep(cell_fill, length(unique_cell) %/% length(cell_fill) + 1)
            names(cell_fill)[seq_len(length(unique_cell))] <- unique_cell
        }
        fill_table[probe] <- cell_fill[as.matrix(table[probe])]
    }

    {
        # 设置viewport
        # GFF非空
        # ROOT---TOP---model
        #         \
        #          ----table
        # GFF为空
        # ROOT---table
        if (equal_col_width){
           tws <- rep(1, nc)
        } else if (TRUE){
            # 获取每一列最宽值
            tws <- c()
            for(i in seq_len(nc)){
                w <- 0
                for(j in seq_len(nr)) w <- max(w, nchar(table[j,i]), na.rm = TRUE)
                tws <- c(tws, w)
            }
        } else {
            tws <- rep(1, nc)
        }
        vpslayout <- grid.layout(nrow = nr + 1, ncol = nc,
                                 widths = tws,
                                 heights = rep(1, nr + 1),
                                 default.units = "null")
        if (! is.null(gff)){
            if (gene_model_height >= 1) stop('gene_model_height should less than 1')
            if (space_height >= 1) stop('space_height should less than 1')
            if (gene_model_height + space_height >= 1){
                warning('gene_model_height and space_height is set to default')
                gene_model_height <- 0.2
                space_height <- 0.1
            }
            if (is.null(table_height)) table_height <- 1 - gene_model_height - space_height

            grid.newpage()
            gmlayout <- grid.layout(nrow = 3,
                                    heights = c(gene_model_height, space_height, table_height),
                                    default.units = "npc")
            vpTOP <- viewport(x = 0.5, y = 0.5, width = 0.9, height = 0.9,
                              default.units = "npc",
                              layout = gmlayout,
                              gp = style,
                              name = "TOP")
            pushViewport(vpTOP)
            vpm <- viewport(layout.pos.row = 1, name = "model", gp = gpar(fill = gene_model_fill))
            pushViewport(vpm)
            # geneModel.plot()
            # 表格内容
            seekViewport("TOP")
            vps <- viewport(x = 0.5, y = 0.5, width = 1, height = 1,
                            default.units = "npc",
                            layout = vpslayout,
                            layout.pos.row = 3,
                            name = 'table')
            pushViewport(vps)
        } else {
            grid.newpage()

            # 表格内容
            vps <- viewport(x = 0.5, y = 0.5, width = 0.9, height = 0.9,
                            default.units = "npc",
                            layout = vpslayout,
                            gp = style,
                            name = 'table')
            pushViewport(vps)
        }
    }


    # 绘制基因model 和 变异位点
    if(! is.null(gff)){
        # 获取染色体位置信息
        opts <- attr(hapSummary, "options")
        if (is.null(Chr)) Chr <- opts["CHROM"]
        pos <- as.numeric(unlist(strsplit(opts["POS"],"-")))
        tmp <- c(start, end)
        if (! is.numeric(start)) start <- min(pos) else start <- min(tmp, na.rm = TRUE)
        if (! is.numeric(end)) end <- max(pos) else end <- max(tmp, na.rm = TRUE)

        seq_names <- unique(as.vector(GenomicRanges::seqnames(gff)))
        if (! Chr %in% seq_names) stop(Chr, " in 'hapSummary' was not found in GFF") else {
            gff <- gff[as.vector(gff@seqnames) %in% Chr]
            gff_s <- min(GenomicRanges::start(gff))
            gff_e <- max(GenomicRanges::end(gff))
            if (start > gff_e | end < gff_s) stop("There was no overlap between hapSummary: ", start, "-", end,
                                                  " and gff: ", gff_s,"-",gff_e, ".")
        }

        # 获取目标区间的注释信息
        gene <- GenomicRanges::GRanges(Chr,
                                       IRanges::IRanges(
                                           start = start,
                                           end = end))

        over <- try(gff[gff %over% gene], silent = TRUE)
        if (length(over) == 0) stop("There is no overlap between GFF and hapSummary")
        probe <- stringr::str_detect(tolower(over$type), "utr|cds")
        over <- over[probe]
        if (length(over) == 0) stop("There is no 'utr' or 'cds' annotation at gene range in GFF")
        over$Parent <- unlist(over$Parent)
        ntranscript <- length(unique(over$Parent))

        # 基础绘图参数(gene model)
        seekViewport("model")

        p_ <- over
        np <- length(unique(over$Parent))
        y_ <- 1/(1 + np)

        # 绘制变异位点
        {
            seekViewport("TOP")
            ps_ <- as.numeric(as.vector(as.matrix(table[1, seq_len(sites(hapSummary)) + 1])))
            ps <- rep(ps_,3)
            ps <- ps[order(ps)]
            ps <- (ps - start)/abs(start - end)
            grid.lines(x = ps,
                       y = rep(c(1 - gene_model_height,
                                 1 - gene_model_height + (y_) * gene_model_height,
                                 NA), length(ps)/2))
            psm <- ps

            # 变异位点与表头连线
            # 表格中每一列的宽度(包含位置)
            tps_r <- cumsum(tws)/sum(tws) # 第n列的右边界
            tps_l <- c(0, tps_r[-length(tps_r)]) # 第n列的左边界
            tps_m <- (tps_r + tps_l)/2 # 第n列的中间位置

            psm[3*seq_len(sites(hapSummary)) - 2] <- tps_m[seq_len(sites(hapSummary)) + 1]
            grid.lines(x = psm,
                       y = rep(c(table_height,
                                 1 - gene_model_height,
                                 NA),
                               length(ps)/2))
        }

        # 绘制基因model
        {
            # Genome Lines
            xs_g <- rep(c(0.01,0.99,NA), np)
            ys_g <- rep(y_ * seq_len(np), 3)
            ys_g <- ys_g[order(ys_g)]
            grid.lines(x = xs_g, y = ys_g, vp = vpm)

            # gene structure
            is_CDS <- stringr::str_detect(tolower(p_$type),"cds")
            is_utr <- stringr::str_detect(tolower(p_$type),"utr")
            hs <- (1 * is_CDS + 0.5 * is_utr) * CDS_height
            ws <- p_@ranges@width/abs(start - end)
            xs <- (p_@ranges@start + p_@ranges@width/2 - start)/abs(start - end)
            ys <- y_ * as.numeric(as.factor(p_$Parent))
            grid.rect(x = xs, y = ys, width = ws, height = hs, default.units = "npc", vp = vpm)

            # grid.text(unique(p_$Parent), x = rep(0.1, np), y = y_ * seq_len(np) + CDS_height, vjust = 0, vp = vpm)
        }

    }

    # 绘制表格
    seekViewport('table')
    for(i in seq_len(nr)){
        for(j in seq_len(nc)){
            vp_name_ij <- paste("vp",i,j,sep = ".")
            rect_name_ij <- paste("rect",i,j,sep = ".")
            txt_name_ij <- paste("txt",i,j,sep = ".")
            assign(vp_name_ij,
                   viewport(layout.pos.row = i,
                            layout.pos.col = j,
                            name = vp_name_ij))

            grid.rect(x = 0.5, y = 0.5, width = 1, height = 1,
                      default.units = "npc",just = "centre",
                      gp = gpar(fill = fill_table[i,j], lty = 1, col = "white"),
                      vp = get(vp_name_ij))
            grid.text(x = 0.5, y = 0.5, just = "centre",table[i,j], vp = get(vp_name_ij))

        }
    }

    if (! is.null(footbar)){
        seekViewport(name = "table")
        current.parent()
        current.viewport()
        i <- nr + 1
        j <- nc
        vp_name_ij <- paste("vp",i,j,sep = ".")
        assign(vp_name_ij,
               viewport(layout.pos.row = i,
                        layout.pos.col = j,
                        name = vp_name_ij))

        grid.text(paste(ft_i, ft_t, footbar, sep = "\n"),x = 1, y = 0.9, default.units = "npc",
                  hjust = 1, vjust = 1,
                  vp = get(vp_name_ij))
    }
    upViewport(0)
    grid.text(title,x = 0.5, y = 0.99, just = "top", gp = gpar(fontsize = 20))
}

