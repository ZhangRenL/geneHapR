#' @name addPromoter
#' @title add promoter to annotation
#' @param anno anotation, imported gff/bed
#' @param PromoterLength the length of promoter region, default as 1500
#' @param bedFile the output bed file name
#' @examples
#' data("geneHapR_test")
#' bed <- addPromoter(gff)
#' @export
addPromoter <-
    function(anno,
             PromoterLength = 1500,
             bedFile = NULL) {
        bed <- gff <- anno
        gff <- gff[gff$type == "CDS"]
        gff$Parent <- unlist(gff$Parent)

        df <- data.frame()
        dfn <- data.frame()
        n <- 0
        for (i in unique(gff$Parent)) {
            n <- n + 1
            gffi <- gff[gff$Parent == i]
            rangesi <- data.frame(gffi@ranges)
            strandi <- unique(gffi@strand) %>% as.vector()
            chri <- gffi@seqnames %>% unique() %>% as.vector()
            if (length(strandi) > 1)
                warning("Can't find strands of ", i)
            if (strandi[1] == "+") {
                endi <- min(rangesi$start) - 1
                starti <- endi - PromoterLength
            } else{
                starti <- max(rangesi$end) - 1
                endi <- starti + PromoterLength
            }
            dfi <-
                c(chri, starti, endi, paste(i, "Promoter"), ".", strandi)
            dfn <- rbind(dfn, dfi)
            if ((n %% 500 == 1) & (n > 1)) {
                message("Processed ", n,
                        ", ",
                        round(n / length(unique(
                            gff$Parent
                        )), 2) * 100,
                        "%")
                names(dfn) <-
                    c("Chr", "start", "end", "name", "X", "strand")
                df <- rbind(df, dfn)
                dfn <- data.frame()
            }
        }
        if (nrow(dfn) != 0) {
            names(dfn) <- c("Chr", "start", "end", "name", "X", "strand")
            df <- rbind(df, dfn)
        }

        tempf <- tempfile()
        on.exit(unlink(tempf))
        write.table(
            df,
            tempf,
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t"
        )
        if (is.null(bedFile))
            return(c(bed, import_bed(tempf)))
        else{
            gff2bed6(c(bed, import_bed(tempf)), bedFile = bedFile)
        }
    }


gff2bed6 <- function(gff,
                     gffFile = NULL,
                     bedFile = NULL) {
    if (!is.null(gffFile))
        gff <- import_gff(gffFile)
    df <- data.frame(gff)
    tmp = c()
    n = 1
    for (i in df$Parent) {
        tmp[n] = if (length(i) == 0)
            df$Name[n]
        else
            i
        n = n + 1
    }
    df$Parent <- tmp

    bed <- data.frame(
        Chr = df$seqnames,
        start = df$start - 1,
        end = df$end,
        name = paste(df$Parent, df$type),
        score = ".",
        strand = df$strand
    )

    if (is.null(bedFile)) {
        tempf <- tempfile()
        on.exit(unlink(tempf))
        write.table(
            bed,
            tempf,
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t"
        )
        import_bed(tempf)
    } else
        write.table(
            bed,
            bedFile,
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t"
        )
}
