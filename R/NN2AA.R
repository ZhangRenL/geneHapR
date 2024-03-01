condon <- c("TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG",
            "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG",
            "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG",
            "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG",
            "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG",
            "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG",
            "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG",
            "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG")
condon_rev <- c("AAA", "AAG", "AAT", "AAC", "AGA", "AGG", "AGT", "AGC",
                "ATA", "ATG", "ATT", "ATC", "ACA", "ACG", "ACT", "ACC",
                "GAA", "GAG", "GAT", "GAC", "GGA", "GGG", "GGT", "GGC",
                "GTA", "GTG", "GTT", "GTC", "GCA", "GCG", "GCT", "GCC",
                "TAA", "TAG", "TAT", "TAC", "TGA", "TGG", "TGT", "TGC",
                "TTA", "TTG", "TTT", "TTC", "TCA", "TCG", "TCT", "TCC",
                "CAA", "CAG", "CAT", "CAC", "CGA", "CGG", "CGT", "CGC",
                "CTA", "CTG", "CTT", "CTC", "CCA", "CCG", "CCT", "CCC")
condon2p3 <- c("Phe", "Phe", "Leu", "Leu", "Ser", "Ser", "Ser", "Ser",
               "Tyr", "Tyr", "***", "***", "Cys", "Cys", "***", "Trp",
               "Leu", "Leu", "Leu", "Leu", "Pro", "Pro", "Pro", "Pro",
               "His", "His", "Gln", "Gln", "Arg", "Arg", "Arg", "Arg",
               "Ile", "Ile", "Ile", "Met", "Thr", "Thr", "Thr", "Thr",
               "Asn", "Asn", "Lys", "Lys", "Ser", "Ser", "Arg", "Arg",
               "Val", "Val", "Val", "Val", "Ala", "Ala", "Ala", "Ala",
               "Asp", "Asp", "Glu", "Glu", "Gly", "Gly", "Gly", "Gly")
condon2p1 <- c("F", "F", "L", "L", "S", "S", "S", "S",
               "Y", "Y", "*", "*", "C", "C", "*", "W",
               "L", "L", "L", "L", "P", "P", "P", "P",
               "H", "H", "Q", "Q", "R", "R", "R", "R",
               "I", "I", "I", "M", "T", "T", "T", "T",
               "N", "N", "K", "K", "S", "S", "R", "R",
               "V", "V", "V", "V", "A", "A", "A", "A",
               "D", "D", "E", "E", "G", "G", "G", "G")
comnucl <- c("A","T","C","G")
names(comnucl)<- c("T","A","G","C")
names(condon2p3) <- condon
names(condon2p1) <- condon
p3_2p1 <- condon2p1; names(p3_2p1) <- condon2p3
p1_2p3 <- condon2p3; names(p1_2p3) <- condon2p1

genome2CDS <- function(ref_seq, gff, ref_start = 1, geneID = ""){
    gene_gff <- gff[1]
    for(i in seq_len(length(gff))){
        if (stringr::str_detect(gff[i]$Parent[[1]], geneID))
            if ((stringr::str_detect(gff[i]$type, "cds")))
                gene_gff <- c(gene_gff, gff[i])
    }
    gene_gff <- gene_gff[-1]

    return(NA)
}

genomePos2CDS <- function(pos, gff){
    return(NA)
}
condon2aa <- function(nucl,aaStyle){
    if (aaStyle == 1) {
        return(condon2p1[nucl])
    } else if (aaStyle == 3) {
        return(condon2p3[nucl])
    } else {
        stop("aaStyle should be 1 or 3")
    }
}

nucl2aa <- function(nucl, aaStyle = 1){
    if (nchar(nucl) %% 3 != 0) stop("Error In Seq Length")
    if (nchar(nucl) == 3) {
        return(condon2aa(nucl,aaStyle))
    } else {
        cons <- c()
        n <- 1
        l <- nchar(nucl)
        while(n < nucl){
            cons <- c(cons, substr(nucl, n, n + 2))
            n <- n + 3
        }
        return(paste0(condon2aa(nucl, aaStyle)))
    }
}
#' @param hap hapResult or hapSummary object
#' @description
#' annotat variants effect in hapResult
#' @details
#' variant effect list:
#'   氨基酸替换
#'   移码变异
#'   可变剪接
#' @param ref_seq reference sequence
#'
hapSnpEff <- function(hap, ref_seq, ref_start = 1, gff, replace_INFO = TRUE, ref_type = c("CDS", "genome"),
                      aaStyle = 3, geneID){
    # TODO:
    # 判断可变剪接
    # 判断基因方向
    if (! any(inherits(hap, "hapSummary"), inherits(hap, "hapResult")))
        stop("hap should be 'hapResult' or 'hapSummary'")
    # 参考序列变为大写字母，检查参考序列
    ref_seq <- as.character(ref_seq[1])
    ref_seq <- toupper(ref_seq)
    cref <- unique(unlist(strsplit(ref_seq, "|")))
    if (! all(cref %in% c("A", "T", "C","G"))) {
        stop("Non-nucleotide letters found: ", cref[! cref %in% c("A", "T", "C","G")])
    }

    # 判断基因方向
    if (tolower(ref_type[1]) == "cds") {
        strand <- "+"
    } else {
        gff_gene <- gff[stringr::str_detect(gff$Parent, geneID)]
        gff_gene_cds <- gff_gene[tolower(gff_gene$type) == "cds"]
        strand <- as.vector(gff_gene@strand@values)[1]
    }

    # CDS gff
    if(tolower(ref_type[1]) != "cds"){
        # cds sequences
        gff_gene_cds
        gff_gene_cds$Parent
        cds_s <- gff_gene_cds@ranges@start
        cds_e <- cds_s + gff_gene_cds@ranges@width - 1
        p <- order(cds_s)
        cds_s <- cds_s[p]
        cds_e <- cds_e[p]
        cds_seq <- c()
        for (i in seq_along(cds_s)){
            cds_seq <- c(cds_seq, substr(ref_seq, cds_s[i] - ref_start, cds_e[i] - ref_start))
        }
        ref_seq <- paste0(cds_seq, collapse = "")
        print("CDS Length: ", nchar(ref_seq))
        if (strand == "-"){
            ref_seq <- rev(unlist(strsplit(ref_seq, "|")))
            ref_seq <- comnucl[ref_seq]
            ref_seq <- paste0(ref_seq, collapse = "")
        }

    }
    # check start condon
    atgStart <- 1 - ref_start
    startCondon <- substr(ref_seq, atgStart, atgStart + 2)
    if(startCondon != "ATG") warning("The start condon should be 'ATG' rather than ", startCondon,"!")
    if (ref_start >= 1) warning("ref_start should be start from 1 when the ref_type is CDS")


    # 处理多态性位点
    for (ci in seq_len(ncol(hap))){
        pos <- suppressMessages(as.numeric(hap[hap$Hap == "POS", ci]))
        if (!is.na(pos)){
            if (tolower(ref_type[1]) == "cds") {
                cdsPos <- pos
                strand <- "+"
            } else {
                cdsPos <- genomePos2CDS(pos, gff)
            }

            allele <- hap[hap$Hap == "ALLELE", ci]
            als <- unlist(strsplit(allele, "[|/,]"))

            # 密码子变化
            ref_len <- nchar(als[1])
            scon <- ceiling(cdsPos/3)
            econ <- ceiling((cdsPos + ref_len - 1)/3)
            condon_ci <- substr(ref_seq, scon * 3 - 2, econ * 3)
            condon_ci_p <- (cdsPos - atgStart) %% ifelse(strand == "+", 3, -3)

            # 密码子变化
            if ((condon_ci_p == 1 & strand == "+") | (condon_ci_p == 0 & strand == "-")) {
                condon_ci = paste0(als, substr(condon_ci, ref_len + 1, ref_len + 2))
            } else if (condon_ci_p == 2 & strand == "+"| condon_ci_p == 2  & strand == "-"){
                condon_ci = paste0(substr(condon_ci, 1, 1),
                                   als,
                                   substr(condon_ci, ref_len + 2, ref_len + 2))
            } else {
                condon_ci <- paste0(substr(condon_ci, 1, 2), als)
            }

            # 差异是3的倍数
            p <- nchar(condon_ci) %% 3 == 0
            aa_ci <- condon_ci
            aa_ci[p] <- nucl2aa(condon_ci[p], aaStyle)
            # 差异不是3的倍数
            aa_ci[!p] <- "Shift"


            # INFO 变化
            nninfo <- paste0("CDS=",
                             condon_ci[1],
                             cdsPos - atgStart + 1,
                             paste(condon_ci[seq_along(condon_ci)[-1]],
                                   sep =  ","))
            aainfo <- paste0("Protein=", aa_ci[1],
                             scon,
                             paste(aa_ci[seq_along(aa_ci)[-1]],
                                   sep =  ","))

            if(replace_INFO) {
                hap[hap$Hap == "INFO", ci] <- paste(nninfo,aainfo, sep = ";")
            } else {
                hap[hap$Hap == "INFO", ci] <- paste(hap[hap$Hap == "INFO", ci], nninfo, aainfo, sep = ";")
            }

        }
    }

    return(hap)


    # # genome2CDS
    # if (tolower(ref_type[1]) != "cds") {
    #     gene_gff <- gff[1]
    #     for(i in seq_len(length(gff))){
    #         if (stringr::str_detect(gff[i]$Parent[[1]], geneID))
    #             if ((stringr::str_detect(gff[i]$type, "cds")))
    #                 gene_gff <- c(gene_gff, gff[i])
    #     }
    #     gene_gff <- gene_gff[-1]
    #
    #     if(length(gene_gff) == 0){
    #         stop("No CDS region for ", geneID, " in Provided annotation,\n",
    #              "please check your 'geneID' and annotation file.")
    #     } else {
    #         cds <- c()
    #         for (i in seq_len(length(gene_gff))){
    #             cds <- c(cds, substr(ref_seq, ref_start, ref_start))
    #         }
    #         cds <- paste(cds, sep = "", collapse = "")
    #     }
    # }
    # ref_seq <- toupper(ref_seq)

}
