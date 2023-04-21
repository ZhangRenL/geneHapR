#' not used yet
merge_Vcf <- function(...){
    vcfs <- list(...)
    for(i in vcfs){
        if(!inherits(i, "vcfR")) stop("Only VCFs supported")
    }

    if(length(vcfs) == 1){
        return(vcfs[[1]])
    } else if(length(vcfs) >= 2){
        vcf <- vcfs[[1]]
        for(i in 2 : length(vcfs)){
            vcfi <- vcfs[[i]]
            if(ncol(vcf@gt) == ncol(vcfi@gt)){
                if(all(colnames(vcf@gt) %in% colnames(vcfi@gt)) &
                   all(colnames(vcfi@gt) %in% colnames(vcf@gt))){
                    vcf@fix <- rbind(vcf@fix, vcfi@fix)
                    vcf@gt <- rbind(vcf@gt, vcfi@gt)
                } else {
                    stop("Samples in all files should be same")
                }
            } else {
                stop("Samples in all files should be same")
            }
        }
    }

    return(vcf)
}

#' not used yet
merge_plink <- function(...){
    files <- list(...)
    for(i in files){
        if(!inherits(i, "list")){stop("Were all inputs p.link ?")}
    }
    if(length(files) == 1){
        return(files[[1]])
    } else if(length(files) >= 2){
        f <- files[[1]]
        for(i in 2 : length(files)){
            fi <- files[[i]]
            if(ncol(vcf@gt) == ncol(vcfi@gt)){
                if(all(colnames(vcf@gt) %in% colnames(vcfi@gt)) &
                   all(colnames(vcfi@gt) %in% colnames(vcf@gt))){
                    vcf@fix <- rbind(vcf@fix, vcfi@fix)
                    vcf@gt <- rbind(vcf@gt, vcfi@gt)
                } else {
                    stop("Samples in all files should be same")
                }
            } else {
                stop("Samples in all files should be same")
            }
        }
    }

    return(vcf)
}
