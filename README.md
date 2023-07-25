# scCancer2

## Introduction

We updated our R toolkit, `scCancer`, based on massive single-cell transcriptome and spatial transcriptome data.

1. `Cell subtype annotation` and `cross-dataset label similarity`: Our analysis mainly focused on cell subtype annotation by training multiple lightweight machine-learning models on scRNA-seq data. We proposed a method for quantitatively evaluating the similarity of cell subtype labels originating from different published datasets. We fully preserved the original labeling in cell atlases and analyzed the relationship between cell subtypes across datasets.

2. `Malignant cell identification`: We constructed a reference dataset combining scRNA-seq and bulk RNA-seq data across multiple cancer types to identify the malignant cell in TME. We trained a model to identify malignant cells with high generalization ability and computational efficiency. 

3. `Spatial transcriptome analysis`: Finally, we integrated a spatial transcriptome analysis pipeline. It enables us to analyze TME from a spatial perspective systematically and automatically.

With `scCancer2`, researchers can understand the composition of the TME more accurately from multiple dimensions.

For old version of `scCancer`, see https://github.com/wguo-research/scCancer.

## Overview of scCancer2
![image](https://github.com/czythu/scCancer/blob/master/inst/Overview.png)

## System Requirements

R version: >= 3.5.0

My R version: R 4.0.5

## Package Installation and Quick Start for scRNA-seq analysis

Quick start of scCancer2:

1. Dependency installation

Some dependent packages for scCancer2 (old version of scCancer, edgeR, garnett, xgboost, and org.Hs.eg.db) may not be installed automatically, so you can install them from the following steps. After installing them successfully, you can run the demos.

```R
checkPkg <- function(pkg){
    return(requireNamespace(pkg, quietly = TRUE))
}

# Some frequently used packages
if(!checkPkg("BiocManager")) install.packages("BiocManager")
if(!checkPkg("devtools")) install.packages("devtools")
if(!checkPkg("Seurat")) BiocManager::install("Seurat")
if(!checkPkg("Biobase")) BiocManager::install("Biobase")
if(!checkPkg("knitr")) BiocManager::install("knitr")
if(!checkPkg("GSVA")) BiocManager::install("GSVA")
if(!checkPkg("pheatmap")) BiocManager::install("pheatmap")
if(!checkPkg("ComplexHeatmap")) BiocManager::install("ComplexHeatmap")
if(!checkPkg("markdown")) install.packages("markdown")
if(!checkPkg("R.utils")) install.packages("R.utils")
if(!checkPkg("RcppArmadillo")) install.packages("RcppArmadillo")
if(!checkPkg("RcppProgress")) install.packages("RcppProgress")

# Algorithms integrated in scCancer1.0
if(!checkPkg("NNLM")) devtools::install_github("linxihui/NNLM")
# if(!checkPkg("harmony")) devtools::install_github("immunogenomics/harmony")
# if(!checkPkg("liger")) devtools::install_github("MacoskoLab/liger")
if(!checkPkg("harmony")) install.packages("harmony")
if(!checkPkg("liger")) install.packages("rliger")

# Algorithms integrated in scCancer2.0
if(!checkPkg("monocle")) BiocManager::install(c("monocle"))
if(!checkPkg("edgeR")) BiocManager::install(c("edgeR"))
if(!checkPkg(c('DelayedArray','DelayedMatrixStats','org.Hs.eg.db','org.Mm.eg.db'))) BiocManager::install(c('DelayedArray','DelayedMatrixStats','org.Hs.eg.db','org.Mm.eg.db'))
if(!checkPkg("garnett")) devtools::install_github("cole-trapnell-lab/garnett")
if(!checkPkg("xgboost")) install.packages("xgboost")
```

if errors occur when installing "NNLM", "edgeR" or "harmony", you may install them from the .tar.gz file:

https://cran.r-project.org/src/contrib/Archive/NNLM/

https://bioconductor.org/packages/release/bioc/html/edgeR.html

We noticed that directly installing harmony from CRAN might meet this bug: [harmony/issues/](https://github.com/immunogenomics/harmony/issues/159),
you may download the source package from https://github.com/immunogenomics/harmony/releases/tag/0.1. to run scCombination.

2. If you have already installed the above dependencies, you have 2 ways to run scCancer2.0:

(a) Recommended: if you want to completely update scCancer to the next version:

```R
# install scCancer2.0
devtools::install_github("czythu/scCancer")
library(scCancer)
```

See [scCancer2.rmd](https://github.com/czythu/scCancer/blob/master/vignettes/) for demos.

(b) Download .zip file of R package. Install [scCancer](https://github.com/wguo-research/scCancer) and run temporary installation in [scCancer2.rmd](https://github.com/czythu/scCancer/blob/master/vignettes/) in the scCancer folder.

```R
# install scCancer1.0
devtools::install_github("wguo-research/scCancer")
# Load all files in the folder
suppressWarnings(load_all())
suppressWarnings(document())
library(scCancer)
```

(c) We provide 5 recommended data sets, including 2 large-scale published datasets (multi-sample) and 3 unpublished data (single-sample).

[`CRC-example-immune (Source: GSE146771)`](https://cloud.tsinghua.edu.cn/f/dc6178e9a37746cf9f11/?dl=1)

[`PAC-example-tumor (Source: CRA001160)`](https://cloud.tsinghua.edu.cn/f/a7b70953a42048ccb231/?dl=1)

[`KC-example-tumor`](https://cloud.tsinghua.edu.cn/f/6b29aab86fc94340832e/?dl=1)

[`PAC-example-normal`](https://cloud.tsinghua.edu.cn/f/3f4715952407477b8b3a/?dl=1)

[`Organoid-example-epithelial`](https://cloud.tsinghua.edu.cn/f/5519909386244a058255/?dl=1)

## A Python module for malignant cell identification

If you are only interested in identifying malignant cells in your own samples or want to reproduce Figure4 in our manuscript, we highly recommend using the pipeline [Figure4.ipynb](https://github.com/czythu/scCancer_MicroEnv/tree/master/MalignantCellIdentification)

The constructed reference data set, several query data sets and the trained model have been uploaded to: ......

The basic processing steps are relied on package `scanpy`, `sklearn` and `xgboost`.

## Report Generation for scRNA-seq analysis

1. The results of cell subtype annotation are stored in folder: cellSubtypeAnno/

2. The results of malignant cell identification by machine learning method are directly insert into original report (report-scAnno.html).

## Package Installation and Quick Start for spatial transcriptome analysis

See https://github.com/Miaoyx323/stCancer for details.

Most of the dependencies are the same, only copyKAT need to be installed.

```R
checkPkg <- function(pkg){
    return(requireNamespace(pkg, quietly = TRUE))
}
if(!checkPkg("copykat")) devtools::install_github("navinlabcode/copykat")
devtools::install_github("Miaoyx323/stCancer")
```

See [stCancer.rmd](https://github.com/czythu/scCancer/blob/master/vignettes/) for demos.

## Citation