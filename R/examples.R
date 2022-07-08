#### 1. I/O ####
# 1.1 import VCF
vcf <- import_vcf("my_vcf_file.vcf")
vcf <- import_vcf("my_vcf_file.vcf.gz")

# 1.2 import gff
gff <- import_gff("my_gff_file.gff")

# 1.3 import seqs (fasta format)
seqs <- import_seqs("my_seqs_file.fa", format = "fasta")

# 1.4 import/export haps
# import
hapSummary <- import_hap("hapSummary.txt")
hapResult <- import_hap("hapResult.txt")

# export
write.hap(hapResult, file = "hapResult.txt")
write.hap(hapSummary, file = "hapSummary.txt")

# 1.5 convert into haplotypes (pegas)
haplo <- as.haplotype(hap = hapSummary)




#### 2. VCF and sequences manipulation ####
# 2.1 VCF manipulation
# 2.1.1 big vcf preprocess
# new VCF file will be saved to disk
# extract a single gene/range from a large vcf
filterLargeVCF(VCFin = "Ori.vcf.gz",
               VCFout = "filtered.vcf.gz",
               Chr = "scaffold_8",
               POS = c(19802,24501),
               override = TRUE)

# extract multi genes/ranges from large vcf
filterLargeVCF(VCFin = "Ori.vcf.gz",          # surfix should be .vcf.gz or .vcf
               VCFout = c("filtered1.vcf.gz", # surfix should be .vcf.gz or .vcf
                          "filtered2.vcf.gz", # length of VCFout should be equal with Chr and POS
                          "filtered3.vcf.gz"),
               Chr = c("scaffold_8",
                       "scaffold_8",
                       "scaffold_7"),
               POS = list(c(19802,24501),
                          c(27341,28949),
                          c(38469,40344)),
               override = TRUE)               # if TRUE, existed file will be override without warning

# 2.1.2 VCF filter
# filter VCF by position
vcf_f1 <- filter_vcf(vcf, mode = "POS",
                     Chr = "scaffold_1",
                     start = 4300, end = 5890)

# filter VCF by annotation
vcf_f2 <- filter_vcf(vcf, mode = "type",
                     gff = gff,
                     type = "CDS")

# filter VCF by position and annotation
vcf_f3 <- filter_vcf(vcf, mode = "both",
                     Chr = "scaffold_1",
                     start = 4300, end = 5890,
                     gff = gff,
                     type = "CDS")

# 2.2 sequences manipulation
# 2.2.1 sequences alinment
seqs <- allignSeqs(seqs)

# 2.2.2 sequences trim
seqs <- trimSeqs(seqs, minFlankFraction = 0.1)




#### 3. calculate haplotypes ####
# 3.1 from VCF
hapResult <- vcf2hap(vcf,
                     hapPrefix = "H",
                     filter_Chr = FALSE, Chr = "scaffold_8",
                     filter_POS = FALSE, startPOS = 603658, endPOS = 604956,
                     hyb_remove = TRUE,
                     na.drop = TRUE)

# 3.2 from sequences
hapResult <- seqs2hap(seqs,
                      Ref = names(seqs)[1],
                      hyb_remove = TRUE,
                      na.drop = TRUE,
                      maxGapsPerSeq = 0.25,
                      hapPrefix = "H")



#### 4. set position of ATG As 0 for visualization (optional) ####
# Note that: after modification the position of ATG are 0, 1 and 2 separately
# 4.1 set ATG position as zero in gff
newgff <- gffSetATGas0(gff = gff, hap = hapResult,
                       geneID = "test1G0387",
                       Chr = "scaffold_1",
                       POS = c(4300, 7910))

# set position of ATG as zero in hapResult/hapSummary
newhap <- hapSetATGas0(gff = gff, hap = hapResult,
                       geneID = "test1G0387",
                       Chr = "scaffold_1",
                       POS = c(4300, 7910))



#### 5. add annotations ####
# add annotations to INFO field
# length of values must be equal to sites number
hapSummary <-addINFO(hapSummary,
                     tag = "PrChange",
                     values = rep(c("C->D", "V->R", "G->N"),3),
                     replace = FALSE, sep = ";")
hapSummary <-addINFO(hapSummary1,
                     tag = "CDSChange",
                     values = rep(c("C->A", "T->C", "G->T"),3),
                     replace = FALSE, sep = ";")


#### 5. Summary hap results and visualization ####
# 5.1 summary hap result
hapSummary <- hap_summary(hapResult)

plotHapTable(hapSummary = hapSummary)
plotHapTable(hapSummary,
             hapPrefix = "H",
             INFO_tag = c("CDSChange", "PrChange"),
             displayIndelSize = 1, angle = 45,
             replaceMultiAllele = TRUE,
             ALLELE.color = "grey90")




