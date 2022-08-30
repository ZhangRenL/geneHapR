root = "../example/snp3kvars-CHR8-25947258-25951166-plink"
#' @title import_plink.pedmap
#' @name import_plink.pedmap
#' @description used for import small p.link file stored in map and ped format
#' @param root this function only support p.link file
#'   format stored in "map" and "ped" format, the file names after removed suffix
#'   should be same with each other.
#' @param sep_ped a character indicate the separation of ped file, default as "\t"
#' @param sep_map a character indicate the separation of map file, default as "\t"
#' @param pedfile,mapfile if `root` is missing then `pedfile` and `mapfile` are needed
#' @return list, contains map information stored in data.frame and ped
#'   information stored in data.frame
#' @usage
#'   import_plink.pedmap(root = root,
#'                       sep_ped = " ", sep_map = "\t",
#'                       pedfile = pedfile, mapfile = mapfile)
#' @inherit plink.pedmap2hap examples
#' @export
import_plink.pedmap <- function(root = root,
                                sep_ped = " ", sep_map = "\t",
                                pedfile = pedfile, mapfile = mapfile){
    if(!missing(root)){
        pedfile <- paste0(root, ".ped")
        mapfile <- paste0(root, ".map")
    }
    ped <- read.table(file <- pedfile, sep = sep_ped, # nrows = 100,
                      header = FALSE)
    map <- read.table(file = mapfile, header = FALSE, sep = sep_map)
    if((ncol(ped) - 6) / nrow(map) != 2)
        warning("ped or map file may corrupted")
    return(list(map = map, ped = ped))
}


#' @name filter_plink.pedmap
#' @title filter_plink.pedmap
#' @description used for filtration of p.link
#' @usage
#'   filter_plink.pedmap(x,
#'                       mode = c("POS", "type", "both"),
#'                       Chr = Chr, start = start, end = end,
#'                       gff = gff, type = type)
#' @param x a list stored the p.link information
#' @param mode filtration mode, one of c("POS", "type", "both")
#' @param Chr the chromosome name, need if mode set as POS or both
#' @param start,end numeric, the range of filtration, and the start should smaller than end
#' @param gff the imported gff object
#' @param type should be in `unique(gff$type)`, usually as "CDS", "genome".
#' @param cusTyp if `type` set as custom, then `cusTyp` is needed
#' @inherit plink.pedmap2hap examples
#' @return list, similar with `x`, but filtered
#' @export
filter_plink.pedmap <- function(x,
                                mode = c("POS", "type", "both"),
                                Chr = Chr, start = start, end = end,
                                gff = gff, type = type, cusTyp = cusTyp){
    on.exit(gc(verbose = FALSE))
    if(missing(Chr)) stop("Chr is missing")
    map <- x$map
    probe <- rep(TRUE, nrow(map))


    if(mode %in% c("POS", "both")){
        if(missing(start) | missing(end)) stop("start or end is missing")
        POS <- c(start, end)
        probe <- probe & as.character(map[, 1]) %in% as.character(Chr)
        probe <- probe & map[, 4] > min(POS) & map[, 4] < max(POS)
    }


    if(mode %in% c("type", "both")){
        if (missing(gff))
            stop("gff is missing!")
        p <- type %in% unique(gff$type)
        m <- paste(unique(gff$type), collapse = "','")
        if (FALSE %in% p)
            stop("type should in c('",m,"')")
        if (type == "custom")
            type <- cusTyp
        if ("genome" %in% type) {
            gff <- gff
        } else {
            gff <- gff[gff$type %in% type]
        }
        POSRange <- POS2GRanges(Chr = "scaffold_1", POS = map[, 4])
        probe <- probe & POSRange %over% gff
    }
    if(! TRUE %in% probe)
        warning("There is no overlaps")

    probe_ped <- c(1:6, which(probe) * 2 + 6, which(probe) * 2 + 5)
    probe_ped <- probe_ped[order(probe_ped)]
    return(list(map = x$map[which(probe), ],
                ped = x$ped[, probe_ped]))
}


#' @name plink.pedmap2hap
#' @title plink.pedmap2hap
#' @description convert p.link format data into hapResult
#' @param p.link list contains p.link information
#' @param hapPrefix prefix of haplotype names
#' @param hyb_remove whether remove accessions contains hyb-sites
#' @param na.drop whether drop accessions contains missing data ("N", NA)
#' @usage
#'   plink.pedmap2hap(p.link,
#'                    hapPrefix = "H",
#'                    hyb_remove = TRUE,
#'                    na.drop = TRUE)
#' @examples
#' \donttest{
#'    pedfile <- system.file()
#'    mapfile <- system.file()
#'    p.link <- import_plink.pedmap(pedfile = pedfile, mapfile = mapfile,
#'                                  sep_map = "\t", sep_ped = "\t")
#'    p.link <- filter_plink.pedmap(p.link, mode = "POS",
#'                                  Chr = "chr08", start = 25948004,
#'                                  end = 25949944)
#'    hapResult <- plink.pedmap2hap(p.link, hapPrefix = "H",
#'                                  hyb_remove = TRUE,
#'                                  na.drop = TRUE)
#' }
#' @return object of hapSummary class
#' @export
plink.pedmap2hap <- function(p.link,
                             hapPrefix = "H",
                             hyb_remove = TRUE,
                             na.drop = TRUE) {
    map <- p.link$map

    # set options
    options <- c(hapPrefix = hapPrefix)
    options <- c(options, CHROM = unique(map[,1]))
    options <- c(options, POS = paste0(range(map[,4]), collapse = "-"))

    POS <- map[,4]
    CHR <- map[,1]
    INFO <- map[,2]
    hapData <- plink.pedmap2hapdata(p.link, POS = POS, Chr = Chr)
    hap <- hapData$hap
    allS_new <- hapData$allS_new
    ALLELE <- hapData$alleles

    # Drop hyb or N
    if (hyb_remove) {
        hap[grepl("|", hap, fixed = T)] <- NA
        hap <- na.omit(hap)
        options <- c(options, hyb_remove = "YES")
    } else
        options <- c(options, hyb_remove = "NO")

    if (na.drop) {
        hap[hap == "N"] <- NA
        hap <- na.omit(hap)
        options <- c(options, NA_remove = "YES")
    } else
        options <- c(options, NA_remove = "NO")

    hap <- assign_hapID(hap, hapPrefix)

    # add infos
    meta <- rbind(c("CHR", CHR, ""),
                  c("POS", POS, ""),
                  c("INFO", INFO, ""),
                  c("ALLELE", ALLELE, "")) %>% data.frame()
    colnames(meta) <- colnames(hap)
    hap <- rbind(meta, hap)
    rownames(hap) <- seq_len(nrow(hap))



    # set attributes
    hap <- remove_redundancy_col(hap)
    class(hap) <- unique(c('hapResult', "data.frame"))
    accAll <- p.link$ped[, 1]
    attr(hap, "AccAll") <- accAll
    accRemain <- hap$Accession[hap$Accession != ""]
    attr(hap, "AccRemain") <- accRemain
    attr(hap, "AccRemoved") <- accAll[!accAll %in% accRemain]
    attr(hap, "options") <- options
    hap2acc <- hap$Accession[-c(1:4)]
    names(hap2acc) <- hap$Hap[-c(1:4)]
    attr(hap, "hap2acc") <- hap2acc
    attr(hap, "freq") <- table(hap$Hap[-c(1:4)])

    return(hap)
}


plink.pedmap2hapdata <- function(p.link, POS = POS, Chr = Chr){
    allS_new <- allS

    # prepare for hapdata
    map <- p.link$map
    ped <- as.matrix(p.link$ped)
    ped.cls <- ncol(ped)
    ped[ped == 0] <- "N"
    hap <- data.frame(row.names = ped[,1])
    alleles <- c()
    for(i in seq_len(ped.cls / 2 - 3)){
        al.i <- ped[, 5 + 2 * i]
        al.i1 <- ped[, 6 + 2 * i]

        # prepare genotype
        col.i <- p <- al.i == al.i1
        col.i[p] <- al.i[p]
        col.i[!p] <- paste(al.i[!p], al.i1[!p], sep = "|")

        # hapdata column
        hap <- cbind(hap, col.i)

        # prepare alleles
        alele.i <- c(al.i, al.i1)
        alele.i[alele.i == "N"] <- NA
        alele.i <- na.omit(alele.i) %>% unique()
        if(length(alele.i) > 2){
            alele.i[2] <- paste0(alele.i[1:2], collapse = "/")
            alele.i <- paste(alele.i[2:length(alele.i)], collapse = ",")
        } else if(length(alele.i) == 2){
            alele.i <- paste(alele.i, collapse = "/")
        } else if(length(alele.i) == 1){
            alele.i <- paste(alele.i, alele.i, sep = "/")
        }
        alleles <- c(alleles, alele.i)
    }
    colnames(hap) <- POS
    hap <- as.matrix(hap)
    allS_new <- table(hap) %>% names()
    allS_new
    return(list(hap = hap, allS_new = allS_new, alleles = alleles))
}


<<<<<<< HEAD
# p.link <- import_plink.pedmap(root, sep_ped = "\t", sep_map = "\t")
# x = p.link
# Chr = x$map[1,1]
# POS = range(x$map[,4])
# p.link = filter_plink.pedmap(x, mode = "POS", Chr = Chr, POS = POS)
=======
p.link <- import_plink.pedmap(root, sep_ped = "\t", sep_map = "\t")
x = p.link
Chr = x$map[1,1]
POS = range(x$map[,4])
p.link = filter_plink.pedmap(x, mode = "POS", Chr = Chr, POS = POS)
>>>>>>> 80ac43a1d1649cf55788d2652d2c01878c835b32
