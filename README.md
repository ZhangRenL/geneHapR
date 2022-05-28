---
Title:  quickHapR
Author: ZhangRenL
Date:   2022.05.27
lastUpdate:   2022.05.20
---

# quickHapR

[TOC]

## 1. 介绍

**通过vcf文件对基因进行单倍型分析**

1.  基于二代测序得到的大量群体的基因型信息进行单倍型分析，

2.  结合表型数据筛选出优异单倍型进行后续研究 [quickHapR](https://gitee.com/zhangrenl/quickHapR)

## 2. 基本逻辑/步骤

1.  导入数据（vcf数据，GFF注释，表型数据）
2.  结合gff文件，或直接指定目标区域，筛选vcf文件，获取目标区域变异信息
3.  通过vcf文件计算单倍型，导出单倍型数据
4.  （可选）筛选单倍型，将ATG起始位置设为0
5.  单倍型数据结合注释信息将变异位点标注在基因模式图上
6.  展示各单倍型之间的进化关系（网络图）
7.  比较不同单倍型之间的表型差异，筛选优势单倍型

## 3 输入输出数据格式及在R中的类型

### 3.1 输入文件

输入文件都属于文本文件，都可以使用文本编辑软件打开查看

#### 3.1.1 vcf

**文件格式：** [VCF](https://www.jianshu.com/p/b2b30b23c866)

**读取方法：**`import_vcf()`或`vcfR::read.vcfR()`均可用于该文件的读取,推荐使用 `import_vcf()`

**R数据类型：**vcfR

#### 3.1.2 gff: 基因组注释文件

**文件格式：**[GFF3](www.gmod.org/wiki/GFF3)

**读取方法：**`import_gff()`或`rtracklayer::import()`均可用于该文件的读取,推荐使用`import_gff()`

**R数据类型：**GenomicRanges

#### 3.1.3 phenos

**注意：**这里用到的表型都应为数量性状

**文件格式：**至少两列，第一列为材料名称，与vcf文件中的individuals对应，之后各列为不同表型；第一行为表型名称，表型名称如果包括除表型名外的其他信息如时间、地点等可以用‘.’与表型名分隔开不同元素间可用‘\_’分隔，如：`plantHeight.2021_BeiJing`，`Weight.female`


**读取方法：**`import_pheno()`或`read.table(file = "", header = TRUE, row.names = 1, check.names = FALSE)`, `read.delim()`均可用于该文件的读取，推荐使用\`import_pheno()


**R数据类型：**data.frame，材料名称作为row.names， 表型名称作为col.names

#### 3.1.4 Accession type


**文件格式：**至少两列，第一列为材料名称；随后不同列为对应的材料类别

    Acc     col     class    
    In1     red     Landrace
    In2     red     Landrace
    In3     blue    Cultivar
    In4     blue    Cultivar
    In5     red     Cultivar

**读取方法：**`read.table(file = "", header = TRUE, row.names = 1)` R数据类型：data.frame，材料名称作为row.names， 分类方式作为col.names

### 3.2 结果文件

#### 3.2.1 hap:

#### 3.2.2 hapResult:

#### 3.2.3 hapNet:

#### 3.2.4 hapTree:

## 4. 安装教程

**注意:**

为避免不必要的安装问题，

1.  Windows10及以上系统推荐使用最新版的[R](https://cran.r-project.org/bin/windows/base/)和[Rtools4.2](https://cran.r-project.org/bin/windows/Rtools/rtools42/rtools.html)

2.  Windows7 系统推荐使用 [R ( = 4.1.3)](https://cran.r-project.org/bin/windows/base/old/4.1.3/)和[Rtools4.0](https://cran.r-project.org/bin/windows/Rtools/rtools40.html) (否则有一定概率会出现依赖的Packages无法正常安装)

### 4.1 **安装准备** 

-   [ ] 安装[Rtools](https://cran.r-project.org/bin/windows/Rtools/)软件，点击[这里](https://cran.r-project.org/bin/windows/Rtools/)下载安装

-   [ ] 安装[git](https://git-scm.com/downloads)软件，点击[这里](https://git-scm.com/downloads)下载安装

-   [ ] 安装R packages: `[devtools](https://devtools.r-lib.org/)`，`[BiocManager](https://www.bioconductor.org/install/)`，安装命令`install.packages(c("devtools", "BiocManager"))`

### 4.2 quickHapR 自动安装方式

``` r
# 自动安装命令
# 1. 选择从gitee安装（国内代码托管平台）
# 注意：该方式需要先安装Rtools和Git两个软件
devtools::install_git("https://gitee.com/zhangrenl/quickHapr")

# 2. 选择从github安装
devtools::install_github("zhangrenl/quickHapr")
```

### 4.3 quickHapR手动安装方式

``` r
# 如果上述命令安装失败可前往Gitee下载预编译R包，选择本地安装
# 注意：该方式安装不保证是最新版
# 本地安装quickHapR前还需手动安装依赖的R packages：
install.packages("BiocManager")
library(BiocManager)
install(c("ggpubr", "vcfR", "tidyverse", "stringr", "resHape2", "randomcoloR",
          "rtracklayer", "trackViewer", "GenomicRanges", "IRanges"))
```

## 5. 使用方法

### 5.1 软件测试

``` r
library(quickHapR)
data("quickHap_test")
hap <- get_hap(vcf, hyb_remove = TRUE, na.drop = TRUE)
hapVsPheno(hap = hap,pheno = pheno, phenoName = "GrainWeight.2021",minAcc = 3)
hapResult <- hap_result(hap)
plotHapTable(hapResult)
phenoResult <- hapVsPheno(hap = hap,
                      pheno = phenos,
                      phenoName = "GrainWeight.2021",
                      minAcc = 3,
                      mergeFigs = TRUE)
plot(phenoResult$figs)
hapNet = get_hapNet(hapResult, 
                    accGroup = accGroup)
plot(hapNet)
plotHapNet(hapNet)

```

### 5.2 软件使用

``` r
# 加载quickHapR
library(quickHapR)
data("quickHap_test") # 加载测试数据,处理自己的数据时不必执行该行

# 设定工作目录
setwd("/your/working/directory")


# 基本数据设定
geneID <- "Seita.1G001600"
Chr <- "scaffold_1"
strand <- "-"
start <- 136756
end <- 144094
hapPreFix <- "H"

# 导入数据
vcf = import_vcf("vcf/rawvcf/Seita.1G001600_136756_144094_-_3k.vcf.gz")
gff = import_gff("gff/Sitalica.gff3")
phenos = import_pheno("phenos.txt")
accGroup <- read.table("accgroup.txt", header = TRUE, check.names = FALSE, row.names = 1)


# 根据GFF注释信息对外显子上的变异位点进行筛选
vcf <- filter_vcf(vcf,                # import_vcf() 导入的vcfR
                  mode = "type",      # 筛选模式：POS/type/both之一
                  Chr = Chr, # 筛选模式为POS或both时必须，染色体名称
                  start = start,      # 开始位置，numeric，筛选模式为POS或both时必须，
                  end = end,          # 结束位置，numeric，筛选模式为POS或both时必须，
                  gff = gff,          # gff，筛选模式为type或both时必须
                  type = "CDS")       # 筛选模式为type或both时必须，CDS/exon/gene/genome之一

# 计算并输出单倍型结果
# hap, data.frame:第一列与最后一列分别固定为Hap和Accession，中间列为位置及对应的基因型
# 前四行为注释信息分别是：CHR，POS，ALLELE,INFO
hap = get_hap(vcf,                   # import_vcf() 导入的vcfR
              hapPrefix = hapPreFix, # 单倍型名称前缀
              filter_Chr = FALSE,    # 筛选染色体选项
              Chr = Chr,             # 通过染色体对vcf信息进行筛选
              filter_POS = TRUE,     # 通过位置进行筛选
              startPOS = start,      # Numeric, 起始位置，通过位置对vcf信息进行筛选
              endPOS = end)          # Numeric, 终止位置，通过位置对vcf信息进行筛选

# hapResult, data.frame: 第一列固定为Hap，最后两列分别固定为Accession和freq，中间列为位置及对应的基因型
# 前四行为注释信息分别是：CHR, POS, ALLELE, INFO
hapResult = hap_result(hap,                    # hap 结果
                       hapPrefix = hapPreFix,  # 单倍型名称前缀
                       out = FALSE,            # 是否输出文件，如果为TRUE， 必须指定输出路径file
                       file = "results/Seita.1G001600_HapResult.txt")  # 输出文件路径（tab分隔的表格）


# 可视化单倍型结果
plotGeneStructure(gff,                # gff注释信息
                  hapResult,          # 单倍型结果
                  Chr = Chr, # 基因所在染色体
                  startPOS = start,  # 基因结构示意图的起始位点
                  endPOS = end,    # 基因结构示意图的终止位置
                  type = "pin",       # SNP类型
                  cex = 1,            # circle大小
                  CDS_h = 0.05,       # 不同基因结构的高度
                  fiveUTR_h = 0.02, 
                  threeUTR_h = 0.01) 
# 单倍型表格
plotHapTable(hapResult,               # 单倍型结果
             hapPrefix = hapPrefix,   # 单倍型前缀（阿拉伯数字前的字母）
             geneID = geneID,         # 基因ID， 作为图表Title
             title.color = "grey90")  # 表头底色

# 单倍型与表型的关联分析
phenoResult = hapVsPheno(hap,         # data.frame:第一列与最后一列分别固定为Hap和Accession，中间列为位置及对应的基因型
                         phenos,      # data.frame: 第一列固定为Accession，随后各列为表型数据，phenoName作为colnames
                         phenoName = "yourPhenoName", # 本次分析中使用的表型名称
                         hapPrefix = hapPrefix,   # 单倍型编号的前缀
                         geneID = geneID,         # 基因ID， 作为表头信息
                         mergeFigs = TRUE,        # 是否将两图融合
                         minAcc = 5)              # 需要分析的单倍型包含的数据量最小值

# plot(phenoResult$fig_pvalue)
# plot(phenoResult$fig_Violin)

plot(phenoResult$figs)


# 获取单倍型进化关系
hapNet = get_hapNet(hapResult,            # 单倍型结果
                    accGroup = accGroup,  # 数据框，品系分组
                    groupName = colnames(accGroup)[1])  # 分组名称，accGroup的列名称之一

# 单倍型网络
plot(hapNet)
plotHapNet(hapNet,
           size = "freq", scale.ratio = 1, cex = 0.8, # circle的大小
           col.link = 1, link.width = 1, lwd = 1,     # 连接线颜色、粗细
           show.mutation = 1,                         # 标注突变位点的方式，（0,1,2,3）之一
           pieCol = pieCol,                           # 颜色向量，饼图颜色，与分组数量长度相同
           addLegend = TRUE)                          # 是否标注图注


```
