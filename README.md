---
Title:  quickHapR
Author: ZhangRenL
Date:   2022.05.20
lastUpdate:   2022.05.14
---

# quickHapR

[TOC]

## 1. 介绍

**通过vcf文件对基因进行单倍型分析**

1.  基于二代测序得到的大量群体的基因型信息进行单倍型分析，

2.  结合表型数据筛选出优异单倍型进行后续研究 [quickHapR](https://gitee.com/zhangrenl/quickHapR)

## 2. 基本逻辑

1.  导入数据（vcf数据，GFF注释，表型数据）
2.  通过vcf文件计算单倍型，导出单倍型数据
3.  单倍型数据结合注释信息将变异位点标注在基因示意图上；展示各单倍型之间的进化关系
4.  比较不同单倍型之间的表型差异，筛选优势单倍型

## 3. 安装教程

**注意:**

为避免不必要的安装问题，

1.  建议Windows10及以上系统使用最新版的R ( \>= 4.2.1)和Rtools4.2

2.  建议Windows7 系统使用 R ( = 4.1.3)和Rtools4.0 (因系统差异，可能会出现依赖的Packages无法正常安装)

### 3.1 **安装准备**

-   [ ] 安装Rtools软件

-   [ ] 安装git软件

-   [ ] 安装R packages: `devtools`，`BiocManager`

### 3.2 quickHapR 自动安装方式

``` r
# 自动安装命令
# 1. 选择从gitee安装（国内代码托管平台）
# 注意：该方式需要先安装Rtools和Git两个软件
devtools::install_git("https://gitee.com/zhangrenl/quickhapr")

# 2. 选择从github安装
devtools::install_github("zhangrenl/quickhapr")
```

### 3.3 quickHapR手动安装方式

``` r
# 如果上述命令安装失败可前往Gitee下载预编译R包，选择本地安装
# 注意：该方式安装不保证是最新版
# 本地安装quickHapR前还需手动安装依赖的R packages：
install.packages("BiocManager")
library(BiocManager)
install(c("ggpubr", "vcfR", "tidyverse", "stringr", "reshape2", "randomcoloR",
          "rtracklayer", "trackViewer", "GenomicRanges", "IRanges"))
```

## 4. 使用方法

### 4.1 软件测试

``` r
library(quickHapR)
data("quickHap_test")
hap <- get_hap(vcf,hyb_remove = TRUE, na.drop = TRUE)
hapVsPheno(hap = hap,pheno = pheno,phenoName = "GrainWeight.2021",minAcc = 3)
hapResult <- hap_result(hap)
plotHapTable(hapResult)
plotHapTable(hapResult)
phenoResult <- hapVsPheno(hap = hap,
                      pheno = pheno,
                      phenoName = "GrainWeight.2021",
                      minAcc = 3,
                      mergeFigs = TRUE)
plot(phenoResult$figs)
```

### 4.2 软件使用

``` r
# 加载quickhapR
library(quickHapR)
data("quickHap_test") # 加载测试数据,处理自己的数据时不必执行该行

# 设定工作目录
setwd("/your/working/directory")

# 导入数据
vcf = import_vcf("Seita.1G001600_136756_144094_-_3k_final.vcf.gz")
gff = import_gff("Yugu1.gff3")
phenos = import_pheno("allPheno.txt")

vcf <- filter_vcf(vcf,                # import_vcf() 导入的vcfR
                  mode = "type",      # 筛选模式：POS/type/both之一
                  Chr = "scaffold_1", # 筛选模式为POS或both时必须，染色体名称
                  start = 136756,     # 筛选模式为POS或both时必须，开始位置
                  end = 144094,       # 筛选模式为POS或both时必须，结束位置
                  gff = gff,          # 筛选模式为type或both时必须，gff
                  type = "CDS") # 筛选模式为type或both时必须，CDS/exon/gene/genome之一
                  
# 计算并输出单倍型结果
# hap,data.frame:第一列与最后一列分别固定为HAP和Accession，中间列为位置及对应的基因型
# 前四行为注释信息分别是：CHR，POS，ALLELE,INFO
hap = get_hap(vcf,                 # import_vcf() 导入的vcfR
              filter_Chr = FALSE,  # 筛选染色体选项
              Chr = "scaffold_1",  # 通过染色体对vcf信息进行筛选
              filter_POS = TRUE,   # 通过位置进行筛选
              startPOS = 136756,   # Numeric, 起始位置，通过位置对vcf信息进行筛选
              endPOS = 144094)     # Numeric, 终止位置，通过位置对vcf信息进行筛选

# hapResult, data.frame: 第一列固定为HAP，最后两列分别固定为Accession和freq，中间列为位置及对应的基因型
# 前四行为注释信息分别是：CHR，POS，ALLELE,INFO
hapResult = hap_result(hap,        # hap 结果
                       hap_prefix = "H",  # 单倍型前缀
                       out = FALSE,# 是否输出文件，如果为TRUE， 必须指定输出路径file
                       file = "results/Seita.1G001600_hapResult.txt")  # 输出文件路径（tab分隔的表格）


# 可视化单倍型结果
plotGeneStructure(gff,                # gff注释信息
                  hapResult,          # 单倍型结果
                  Chr = "scaffold_1", # 基因所在染色体
                  startPOS = 136756,  # 基因结构示意图的起始位点
                  endPOS = 144094,    # 基因结构示意图的终止位置
                  type = "pin",       # SNP类型
                  cex = 1,            # circle大小
                  CDS_h = 0.05,       # 不同基因结构的高度
                  fiveUTR_h = 0.02, 
                  threeUTR_h = 0.01) 
plotHapTable(hapResult,               # 单倍型结果
             hapPrefix = "H",         # 单倍型前缀（阿拉伯数字前的字母）
             geneID = "",             # 基因ID， 作为图表Title
             title.color = "grey90")  # 表头底色

# 单倍型与表型的关联分析
phenoResult = hapVsPheno(hap,        # data.frame:第一列与最后一列分别固定为HAP和Accession，中间列为位置及对应的基因型
                 phenos,      # data.frame: 第一列固定为Accession，随后各列为表型数据，phenoName作为colnames
                 phenoName = "yourPhenoName", # 本次分析中使用的表型名称
                 hapPrefix = "H",             # 单倍型编号的前缀
                 geneID = "Seita.1G000000",   # 基因ID， 作为表头信息
                 mergeFigs = TRUE,    # 是否将两图融合
                 minAcc = 5)          # 需要分析的单倍型包含的数据量最小值
                 
# plot(phenoResult$fig_pvalue)
# plot(phenoResult$fig_Violin)

plot(phenoResult$figs)
```
