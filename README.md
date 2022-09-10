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

[Introduction of geneHapR](https://gitee.com/zhangrenl/genehapr/wikis/Introduction.html)

[The workflow of geneHapR](https://gitee.com/zhangrenl/genehapr/wikis/workflow.html)

[DATA Format of geneHapR](https://gitee.com/zhangrenl/genehapr/wikis/data_format.html)

[An example in *Setaria italica*](https://gitee.com/zhangrenl/genehapr/wikis/An_example)

## Issues

If you find any bugs or suggesstions, please add a new issue at [**zhangrenl/genehapr**](https://gitee.com/zhangrenl/genehapr/issues)