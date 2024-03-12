# 汉字转unicode 编码http://tools.bugscaner.com/unicode/
#' @title plotHapTable2
#' @description
#' plot the hapResult in table like style using grid system.
#' This function is under development and may not stable.
#' Some parameters may deleted or renamed in future.
#' @examples
#' #
#' data(geneHapR_test)
#' plotHapTable2(hapResult)
#' plotHapTable2(hapResult, gff = gff)
#' @param hapSummary the hapSummary or hapResult object
#' @param gff,filter_gff_by_geneID gff or bed annotation, filter_gff_by_geneID: filter gff annotation by geneID if provided.
#' @param title title of plot
#' @param CDS_height a numeric vector specified the height of CDS, and the height of utr is half of that,
#'  only useful when gff is provided,
#' @param equal_col_width a bool or numeric vector specified whether column with should be equal
#' @param annot_fill,head_fill,labels_fill the fill color of annotation, head and label row or columns
#' @param table_line_col the line color in genotype table
#' @param INFO_split,INFO_tag,geneID,tag_field used to set annotation in haplotype table.
#'  And the geneID was used to fileter annotation in INFO field.
#' @param footbar the foot notes
#' @param cell_fill a color vector or function or named vector specified cell fill color
#' @param replaceMultiAllele replace multi-allele title by 'T1,T2,...' or not
#' @param col_annots,head_anno,headrows,row_annots,row_labels the column or row number of annotation or labels or heads
#' @param show_INFO show annotation field or not, default as `FALSE`
#' @param show_chr_name show chromosome name at left-top cell or not
#' @param show_indel_size the Indel length longer will be replaced by "i1,i2,i3,..."
#' @param style see help(gpar)
#' @param Chr,start,end which range should be plotted in gene model
#' @param gene_model_height,table_height,space_height the plotting range height of gene model, table and spacer
#' @param model_rect_col,model_rect_fill,model_line_col a string specified the color for line/rectangle in gene model
#' @param model_anno_pos,model_anno_adj,model_anno_cex,model_anno_col,model_anno_txt
#' the position (x,y), just (hjust, vjust), color, size and content of annotation text in gene model
#' @param table_txt_col,table_txt_cex controls the color and size of texts in genotype table
#' @param link_line_type the type of link lines for mutations in gene model and genotype table
#' @param angle the angle of number positions
#' @param annot_for_each_transcrips mark the strand and trancripts name for each gene modle
#' @param ... param not used
#' @export
plotHapTable2 <- function(hapSummary, show_indel_size = 1, replaceMultiAllele = TRUE, angle = 0,
                          show_INFO = FALSE, INFO_split = c(";", ",","\\|"),INFO_tag = "ANN",geneID = NA, tag_field = -1,
                          title = "", filter_gff_by_geneID = TRUE,
                          gff = NULL,show_chr_name = TRUE,
                          Chr = NULL, start = NULL, end = NULL,
                          model_rect_col = "black", model_rect_fill = "grey80", model_line_col = "black",
                          model_anno_txt = NULL, model_anno_col = "black", model_anno_cex = 1,
                          table_txt_col = "black", table_txt_cex = 1,
                          model_anno_pos = c(-1, -1), model_anno_adj = c(0,1),
                          gene_model_height = 0.2, space_height = 0.1, table_height = NULL, CDS_height = 0.3,
                          link_line_type = 3,
                          headrows = 1, equal_col_width = FALSE,
                          head_anno = 1, col_annots = 0,
                          row_labels = 1, row_annots = 1,
                          table_line_col = "white", annot_for_each_transcrips = TRUE,
                          labels_fill = "white", annot_fill = "grey90", head_fill = NULL,
                          cell_fill = NULL, style = gpar(fontfamily = "sans", fontface = 1, cex = 0.7),
                          footbar = "", ...){
    # function start
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

        # 处理Indel
        ft_i <- c()
        for (i in seq_len(length(als))){
            al_i <- als[i]
            table[table == al_i] <- paste0("i",i )
            ft_i <- c(ft_i, paste0("i",i))
        }


        # 处理表头中的Indel序列
        names(als) <- ft_i
        ft_i <- paste(ft_i, als, sep = ":", collapse = ";")
        als <- als[order(nchar(als), decreasing = TRUE)]
        # print(als)
        # print(ft_i)
        # print(ALLELE)
        for(i in seq_len(length(als))){
            ALLELE <- stringr::str_replace(ALLELE, als[i], names(als)[i])
        }


        # 处理多等位位点
        if(replaceMultiAllele & sum(probe_mula) > 0){
            ft_t <- paste("T",seq_len(sum(probe_mula)),":",ALLELE[probe_mula], sep = "", collapse = ";")
            ALLELE[probe_mula] <- paste0("T",seq_len(sum(probe_mula)))
        } else ft_t <- NULL
        table[table$Hap == "ALLELE",] <- ALLELE
    }

    # 显示INFO信息处理
    if(show_INFO) {
        INFO <- as.vector(as.matrix(table[table$Hap == "INFO", ]))
        for (i in seq_len(length(INFO))){
            INFO_i <- INFO[[i]]
            if (is.na(INFO_i)){
                INFO[[i]] <- NA
                next
            }

            # 第一次分离出tags
            INFO_i <- unlist(strsplit(INFO_i, INFO_split[1]))
            p_INFO_i <- stringr::str_detect(INFO_i,paste0(INFO_tag, "="))
            if(any(p_INFO_i)){
                INFO_i <- stringr::str_remove(INFO_i,paste0(INFO_tag, "="))
                INFO_i <- INFO_i[p_INFO_i]
            } else {
                INFO[[i]] <- NA
                next
            }


            # tags
            if (! is.na(INFO_i) & ! is.na(INFO_split[2]) & ! is.na(geneID)){
                INFO_i <- unlist(strsplit(INFO_i, INFO_split[2]))
                INFO_i <- INFO_i[stringr::str_detect(INFO_i, geneID)]
            }


            # tag field
            if (tag_field > 0 & ! any(is.na(INFO_i))){
                tmp <- c()
                for (j in INFO_i){
                    splt <- INFO_split[length(INFO_split)]
                    tmp <- c(tmp, unlist(strsplit(j, splt))[tag_field])

                }
                INFO_i <- paste(unique(tmp), sep = ",",collapse = "\n")
            }
            # print(INFO_i)
            INFO[[i]] <- INFO_i
        }
        table[table$Hap == "INFO",] <- INFO
    }



    {
        # 生成填充表
        if(isFALSE(cell_fill)) {
            cell_fill <- 0
            head_fill <- 0
            annot_fill <- 0
            labels_fill <- 0
            if(table_line_col == "white") table_line_col <- "grey10"
        } else {
            if(is.null(cell_fill)) cell_fill <- OPTS$pie.colors.function
            head_fill <- if (is.null(head_fill)) labels_fill else head_fill
        }

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

        # 获取染色体位置信息
        opts <- attr(hapSummary, "options")
        if (is.null(Chr)) Chr <- opts["CHROM"]
        pos <- as.numeric(unlist(strsplit(opts["POS"],"-")))
        nwchar_fstcol <- max(nchar(table[1,1]) * cospi(angle/180),
                             nchar(table[nrow(table) - 2,1]),
                             nchar(table[2,1])) * table_txt_cex
        nwchar_lstcol <- max(nchar(table[1,ncol(table)]) * cospi(angle/180),
                             nchar(table[ifelse(show_INFO, 4,3),ncol(table)]))  * table_txt_cex
        tws_old <- tws
        tws <- unit.c(unit(nwchar_fstcol + 1, "strwidth", "H"),
                      unit(tws[2:(length(tws) - 1)],"null"),
                      unit(nwchar_lstcol + 1, "strwidth", "5"))

        vpslayout <- grid.layout(nrow = nr + 1, ncol = nc,
                                 widths = tws,
                                 heights = unit.c(
                                     # unit(sinpi(angle / 180) + 0.1, "strwidth", max(pos)),
                                     unit.pmax(unit(1, "strwidth", table[1,])) * table_txt_cex,
                                     unit(rep(1, nr), "null")),
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
            vpm <- viewport(layout.pos.row = 1, name = "model")
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

        # 转录本数量
        par <- c()
        for (i in seq_len(length(over))){
            par <- c(par, unlist(over$Parent[i])[1])
        }
        over$Parent <- par
        if (!is.na(geneID) & filter_gff_by_geneID){
            lower_par <- tolower(par)
            over <- over[stringr::str_detect(lower_par, tolower(geneID))]
        }
        ntranscript <- length(unique(over$Parent))
        # print(over)
        # print(ntranscript)
        # print(over$Parent)






        # 基础绘图参数(gene model)
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
            psm <- unit(ps, "npc")

            # 变异位点与表头连线
            # 表格中每一列的宽度(包含位置)
            # 第n列的中间位置
            # 基因型列的每列宽度
            tws_cumsum <- cumsum(tws_old)
            # 基因型列的每字符宽度
            tws_sum <- sum(tws_old[-c(1, length(tws_old))])
            gtc_wid_per_s <- (unit(1, "npc") -  tws[1] - tws[length(tws)])/tws_sum


            tps_m <- unit(1,"npc")
            for (ci in 1:sites(hapSummary)){
                # 第一列宽度加平均每字符宽度
                tps_m <- unit.c(tps_m, tws[1] + gtc_wid_per_s * (tws_cumsum[ci] + tws_old[ci + 1]/2 - tws_old[1]))
            }

            #
            psm[3*seq_len(sites(hapSummary)) - 2] <- unit(tps_m[seq_len(sites(hapSummary)) + 1], "npc")

            grid.lines(x = psm,
                       y = rep(c(table_height,
                                 1 - gene_model_height,
                                 NA),
                               length(ps)/2),
                       gp = gpar(lty = link_line_type))

        }

        # 绘制基因model
        {
            # Genome Lines
            xs_g <- rep(c(0, 1, NA), np)
            ys_g <- rep(y_ * seq_len(np), 3)
            ys_g <- ys_g[order(ys_g)]
            grid.lines(x = xs_g, y = ys_g,
                       gp = gpar(col = model_line_col),
                       vp = vpm)

            # gene structure
            is_CDS <- stringr::str_detect(tolower(p_$type),"cds")
            is_utr <- stringr::str_detect(tolower(p_$type),"utr")
            hs <- (1 * is_CDS + 0.5 * is_utr) * CDS_height / ntranscript
            ws <- p_@ranges@width/abs(start - end)
            xs <- (p_@ranges@start + p_@ranges@width/2 - start)/abs(start - end)
            ys <- y_ * as.numeric(as.factor(p_$Parent))
            grid.rect(x = xs, y = unit(1 - ys, "npc"), width = ws, height = hs,
                      default.units = "npc",
                      gp = gpar(fill = model_rect_fill, col = model_rect_col), vp = vpm)

            # 注明每条转录本的方向
            strands <- c()

            for (i in unique(over$Parent)){
                strands <- c(strands, as.vector(over@strand[over$Parent == i]@values))
            }
            if(annot_for_each_transcrips){
                # print(ys)
                # print(np)
                # print(paste(ifelse(unique(strands) == "+", "5'->3'", "3'<-5'"), unique(p_$Parent), sep = ";"))
                strands[is.na(strands)] = ""
                grid.text(paste(ifelse(strands == "+", "5'->3'", "3'<-5'"),
                                unique(p_$Parent), sep = ";"),
                          x = unit(rep(0, np), "npc"),
                          # y = unit(unique(ys) + max(hs), "npc") + unit(0.3, "strheight", "53"),
                          y = unit(1, "npc") - unit(unique(ys) - max(hs)/2, "npc") + unit(0.15, "strheight","1"),
                          vjust = 0, # 垂直0:下对齐1:上对齐
                          hjust = 0, # 水平
                          gp = gpar(col = model_anno_col, cex = model_anno_cex),
                          vp = vpm)
            } else {

                # 注明方向
                seekViewport("model")
                if (is.null(model_anno_txt)){
                    if (length(unique(over@strand)) > 1) {
                        warning("strand is ambiguous in annotation, Please spcified the strand")
                    } else {
                        strand <- as.vector(unique(over@strand)[1])
                        model_anno_txt <- ifelse(strand == "+", "5' -> 3'", "3' <- 5'")
                    }
                }


                # model annotation
                if (! is.null(model_anno_txt)){
                    if(model_anno_pos[1] != -1){
                    model_anno_x <- model_anno_pos[1]
                    model_anno_y <- model_anno_pos[2]
                    } else {
                        model_anno_x = 0
                        model_anno_y = 1
                    }
                    model_anno_h <- model_anno_adj[1]
                    model_anno_v <- model_anno_adj[2]
                    # message("model_anno_txt: ", model_anno_txt)
                    # message("model_anno_x: ", model_anno_x)
                    # message("model_anno_y: ", model_anno_y)

                    grid.text(model_anno_txt,
                              x = model_anno_x, y = model_anno_y,
                              hjust = model_anno_h, vjust = model_anno_v,
                              gp = gpar(col = model_anno_col, cex = model_anno_cex),
                              vp = vpm)

                }
            }
        }


    }

    # 填充表格文字
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

            if (i == 1){
                grid.rect(x = 0.5, y = 0.5, width = 1, height = 1,
                          default.units = "npc",just = "centre",
                          gp = gpar(fill = fill_table[i,j], lty = 1,
                                    col = table_line_col,
                                    alpha = ifelse(head_fill == 0, 1, 0)),
                          vp = get(vp_name_ij))

                grid.text(x = 0.5, y = 0.5, just = "centre",table[i,j],
                          gp = gpar(col = table_txt_col, cex = table_txt_cex),
                          vp = get(vp_name_ij), rot = angle)

            }else{
                grid.rect(x = 0.5, y = 0.5, width = 1, height = 1,
                          default.units = "npc",just = "centre",
                          gp = gpar(fill = fill_table[i,j], lty = 1, col = table_line_col),
                          vp = get(vp_name_ij))

                grid.text(x = 0.5, y = 0.5, just = "centre",table[i,j],
                          gp = gpar(col = table_txt_col, cex = table_txt_cex),
                          vp = get(vp_name_ij))
            }
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
        if (!is.null(ft_t))  ft_i <- paste(ft_i, ft_t, sep = "\n")
        grid.text(paste(ft_i, footbar, sep = "\n"),x = 1, y = 0.9, default.units = "npc",
                  hjust = 1, vjust = 1,
                  vp = get(vp_name_ij))
    }
    upViewport(0)
    grid.text(title,x = 0.5, y = 0.99, just = "top", gp = gpar(fontsize = 20))
}

