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

## Package Installation and Quick Start for scRNA-seq analysis

### System Requirements

R version: >= 3.5.0

We have test scCancer2 on R 4.0.5 (Recommended), R 4.1.1 (Recommended) and R 4.2.0 (Workable in most samples)

R tools need to be previously installed: https://cran.r-project.org/bin/windows/Rtools/

To avoid the version conflicts of R packages, we highly recommend that you install a brand new R environment and switch it in anaconda or RStudio.

### Quick start of scCancer2

#### Dependency installation

You can install the dependencies from the following steps. After installing them successfully, you can run the demos.

```R
checkPkg <- function(pkg){
    return(requireNamespace(pkg, quietly = TRUE))
}

# Some frequently used packages
if(!checkPkg("BiocManager")) install.packages("BiocManager")
if(!checkPkg("devtools")) install.packages("devtools")
# Our Seurat version: 4.0.2
if(!checkPkg("Seurat")) BiocManager::install("Seurat")
if(!checkPkg("Biobase")) BiocManager::install("Biobase")
if(!checkPkg("GSVA")) BiocManager::install("GSVA")
if(!checkPkg("pheatmap")) BiocManager::install("pheatmap")
if(!checkPkg("ComplexHeatmap")) BiocManager::install("ComplexHeatmap")
if(!checkPkg("markdown")) install.packages("markdown")
if(!checkPkg("R.utils")) install.packages("R.utils")

# Algorithms integrated in scCancer
if(!checkPkg("NNLM")) devtools::install_github("linxihui/NNLM")
if(!checkPkg("harmony")) install.packages("harmony")
if(!checkPkg("liger")) install.packages("rliger")

# Algorithms newly integrated in scCancer2
if(!checkPkg("xgboost")) install.packages("xgboost")
if(!checkPkg("DESeq2")) BiocManager::install("DESeq2")
# https://cole-trapnell-lab.github.io/garnett/docs/
if(!checkPkg("monocle")) BiocManager::install(c("monocle"))
if(!checkPkg("edgeR")) BiocManager::install(c("edgeR"))
BiocManager::install(c('DelayedArray', 'DelayedMatrixStats', 'org.Hs.eg.db', 'org.Mm.eg.db'))
if(!checkPkg("garnett")) devtools::install_github("cole-trapnell-lab/garnett")
```

Notice: if errors occur when installing "NNLM", "edgeR" or "harmony", you may install them from the .tar.gz file:

https://cran.r-project.org/src/contrib/Archive/NNLM/

https://bioconductor.org/packages/release/bioc/html/edgeR.html

We also noticed that directly installing harmony from CRAN might meet this bug: [github/harmony/issues/](https://github.com/immunogenomics/harmony/issues/159),
you may download and install the source package from https://github.com/immunogenomics/harmony/releases/tag/0.1. to run `scCombination` with harmony method smoothly.

#### Run scCancer2

If you have already installed the above dependencies, you have 2 ways to run scCancer2:

(a) Recommended: if you want to completely update scCancer to the next version:

```R
# install scCancer2
devtools::install_github("czythu/scCancer")
library(scCancer)
# Run demos
# ...... See vignettes/scCancer2.rmd
```

See [scCancer2.rmd](https://github.com/czythu/scCancer/blob/master/vignettes/) for demos.

(b) Download .zip file of R package. Open scCancer.rproj and run temporary installation in [scCancer2.rmd](https://github.com/czythu/scCancer/blob/master/vignettes/).

```R
# Check the dependencies
suppressMessages(library(Seurat))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(garnett))
suppressMessages(library(xgboost))
suppressMessages(library(ggplot2))
suppressMessages(library(ggsci))
suppressMessages(library(cowplot))
suppressMessages(library(viridis))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(edgeR))
library(devtools)
# Load all files in the folder
suppressWarnings(load_all())
suppressWarnings(document())
library(scCancer)
# Run demos
# ...... See vignettes/scCancer2.rmd
```

(c) We have uploaded 5 recommended data sets, including 3 unpublished data (single-sample) and 2 large-scale published datasets (multi-sample).
The first three data sets are recommended for reproducing the whole pipeline because there are fewer samples for faster operation. 
The last two data sets are recommended for cell subtype annotation task because there are more cell types and richer cell numbers.

[`KC-example-tumor`](https://cloud.tsinghua.edu.cn/f/6b29aab86fc94340832e/?dl=1)

[`PAC-example-normal`](https://cloud.tsinghua.edu.cn/f/3f4715952407477b8b3a/?dl=1)

[`Organoid-example-epithelial`](https://cloud.tsinghua.edu.cn/f/5519909386244a058255/?dl=1)

[`CRC-example-immune (Source: GSE146771)`](https://cloud.tsinghua.edu.cn/f/dc6178e9a37746cf9f11/?dl=1)

[`PAC-example-tumor (Source: CRA001160)`](https://cloud.tsinghua.edu.cn/f/a7b70953a42048ccb231/?dl=1)

### R module for newly implemented scRNA-seq analysis

If you have a processed dataset (matrix or Seurat object), you can use cell subtype annotation and malignant cell identification module alone.

`scStatistics` and `scAnnotation` are not needed. See [cellSubtypeAnno.Rmd](https://github.com/czythu/scCancer/blob/master/vignettes/) and [malignantCellIden.Rmd](https://github.com/czythu/scCancer/blob/master/vignettes/) in vignettes folder for tutorials.


### A Python module for malignant cell identification

If you are only interested in identifying malignant cells in your own samples or want to reproduce Figure4 in our manuscript, we highly recommend using the pipeline [Figure4.ipynb](https://github.com/czythu/scCancer_MicroEnv/tree/master/MalignantCellIdentification)

The constructed reference data set will be uploaded later.

The 5 recommended data sets above can be served as query data sets.

The trained model: [`sc_xgboost_alldata.model`](https://github.com/czythu/scCancer_MicroEnv/tree/master/MalignantCellIdentification/model). It has been integrated into the R package.

The basic processing steps are relied on package `scanpy`, `sklearn` and `xgboost`.

Due to the differences between Seurat and scanpy and the parameters setting at the preprocessing steps, the results of malignant cell identification are slightly different in R and Python.

### Report Generation for scRNA-seq analysis

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

[1] Wenbo Guo, Dongfang Wang, Shicheng Wang, Yiran Shan, Changyi Liu, Jin Gu, scCancer: a package for automated processing of single-cell RNA-seq data in cancer, Briefings in Bioinformatics, bbaa127, https://doi.org/10.1093/bib/bbaa127

