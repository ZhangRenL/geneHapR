.onLoad <- function(libname, pkgname) {
    message("")
    message("lastupdate: 2022.6.5")
    message("Author: Zhang Ren Liang\n(Mail: zhang_renliang@163.com)")
    message("")
    options(stringsAsFactors = FALSE)
}

packageStartupMessage("
lastupdate: 2022.5.30
Author: Zhang Ren Liang\n(Mail: zhang_renliang@163.com)
")

homo <- c("A", "C", "G", "T")
names(homo) <- paste(homo, homo, sep = "/")
unknown <- c(".", "N")
names(unknown) <- paste(unknown, unknown, sep = "/")
hetero <- c("A|C","A|G","A|T",
            "C|A","C|T","C|G",
            "G|A","G|C","G|T",
            "T|A","T|C","T|G")
names(hetero) <- stringr::str_replace(hetero,"[|]", "/")

allS <- list(homo = homo,
             unknown = unknown,
             hetero = hetero,
             all = c(homo, hetero))

# B: C G T
# D: A G T
# H: A C T
# V：A C G
# N: A C G T

update_allS <- function(allS_new, REF = REF, ALT = ALT){
    ALT <- stringr::str_split(ALT, pattern = ",")
    ALT <- unlist(ALT)
    # update homo
    for(l in c(REF, ALT)){
        if(!l %in% allS_new$homo){
            homonew <- l
            allS_new$homo <- c(allS_new$homo, homonew)
        }
    }
    allS_new$homo <- unique(allS_new$homo)
    names(allS_new$homo) <- paste(allS_new$homo,allS_new$homo, sep = "/")

    # update hetero
    for(i in c(REF,ALT)){
        for(j in c(REF,ALT)){
            if(i != j){
                hetero_ij <- paste(i,j,sep = "|")
                allS_new$hetero <- c(allS_new$hetero, hetero_ij)
            }
        }
    }
    allS_new$hetero <- unique(allS_new$hetero)
    names(allS_new$hetero) <- stringr::str_replace(allS_new$hetero, "[|]", "/")

    # update allS_new$all
    allS_new$all <- c(allS_new$homo, allS_new$hetero, allS_new$unknown)

    return(allS_new)
}
