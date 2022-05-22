.onLoad <- function(libname, pkgname) {
    options(stringsAsFactors = FALSE)
}

homo <- c("A", "C", "G", "T", "-", "+")
names(homo) <- c("AA", "CC", "GG", "TT", "--", "++")
unknown <- c(".", "N")
hetero <- c("W", "W", "S", "S", "R","R","Y","Y","K","K","M","M","+-","-+")
names(hetero) <- c("AT", "TA", "CG", "GC", "AG", "GA",
                   "CT","TC","GT","TG","AC","CA","+-","-+")
allS <- c(homo, hetero)
# B: C G T
# D: A G T
# H: A C T
# V：A C G
# N: A C G T
