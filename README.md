# `geneHapR`:

## Installation

Windows, Linux and mac users:

```r
# a few non-CRAN packages may need install manually
if(! require(BiocManager)) 
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "GenomicRanges", "muscle", "IRanges", "rtracklayer", "trackViewer"))
install.packages("geneHapR")
if(! require(devtools)) 
    install.packages("devtools")
devtools::install_git("https://gitee.com/zhangrenl/genehapr")  # need devtools package and git software
```

For linux system with docker:

```bash
# pull the docker from hub.docker.com
docker pull zhangrenl/genehapr:latest

# run the image
docker run -ti --rm zhangrenl/genehapr bash
```



## Vignettes 

[Introduction of geneHapR](https://gitee.com/zhangrenl/genehapr/wikis/Introduction.md)

[An example in *Setaria italica*](https://gitee.com/zhangrenl/genehapr/wikis/An_example_in_Setaria_italica)

[An example with lazy data](https://gitee.com/zhangrenl/genehapr/wikis/An%20example%20using%20lazy%20data%20in%20geneHapR)

## Issues

If you find any bugs or suggestions, please add a new issue at [**zhangrenl/genehapr**](https://gitee.com/zhangrenl/genehapr/issues)

## Questions

Please join us using **QQ Group** by ID 255742992,

Or contact me using **QQ** by ID 1205654509

Or contact me using **weChat** by ID suiyuanhongzhu

Or send an email to "1205654509@qq.com"

## Citation

If you using geneHapR in your work, please cite the following paper.

"Zhang, R., Jia, G. & Diao, X. geneHapR: an R package for gene haplotypic statistics and visualization. BMC Bioinformatics 24, 199 (2023). https://doi.org/10.1186/s12859-023-05318-9"

## Acknowledgement

Thanks to Dr. Li Dongdong, Dr. He Qiang and Dr. Zhang Linlin, Pang Jianzhou for their valuable discussions and comments.

Thanks to the users of geneHapR for their feedback suggestions, especially users from Diao's Lab,
H.L., Li Zecong, Dr. Aashish Gyawali, et al.
