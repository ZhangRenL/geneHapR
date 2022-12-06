# `geneHapR`:

## Installation
```r
# a few non-CRAN packages may need install manually
if(! require(BiocManager)) 
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "GenomicRanges", "muscle", "IRanges", "rtracklayer", "trackViewer"))
install.packages("geneHapR")

devtools::install_git("https://gitee.com/zhangrenl/genehapr")  # need devtools package and git software
```


## Vignettes 

[Introduction of geneHapR](https://gitee.com/zhangrenl/genehapr/wikis/Introduction.md)

[An example in *Setaria italica*](https://gitee.com/zhangrenl/genehapr/wikis/An_example_in_Setaria_italica)

[An example with lazy data](https://gitee.com/zhangrenl/genehapr/wikis/An%20example%20using%20lazy%20data%20in%20geneHapR)

## Issues

If you find any bugs or suggesstions, please add a new issue at [**zhangrenl/genehapr**](https://gitee.com/zhangrenl/genehapr/issues)