# scCancer

## Introduction

The `scCancer` package focuses on processing and analyzing droplet-based scRNA-seq data for cancer research. Except basic data processing steps, this package takes several special considerations for cancer-specific features.

The workflow of  `scCancer` mainly consists of three modules: `scStatistics`, `scAnnotation`, and `scCombination`.
* The `scStatistics` performs basic statistical analyses of raw data and quality control.
* The `scAnnotation` performs functional data analyses and visualizations, such as low dimensional representation, clustering, cell type classification, cell malignancy estimation, cellular phenotype analyses, gene signature analyses, cell-cell interaction analyses, etc.
* The `scCombination` perform multiple samples data integration, batch effect correction and analyses visualization.

After the computational analyses, detailed and graphical reports were generated in user-friendly HTML format.

<img src="http://lifeome.net/software/sccancer/scCancer-workflow.png" width="70%" alt="scCancer-workflow" align=center>

([Click to view larger workflow picture](http://lifeome.net/software/sccancer/scCancer-workflow.png))

## System Requirements
* R version: >= 3.5.0 (**suggest:** R 3.6, **not 4.0**)
* **Hint:  For R (version>=4.0) under Windows system**, the Rtools needs to be updated to version 4.0 from https://cran.r-project.org/bin/windows/Rtools/. So, if you are not familiar with R environment configuration, we **don't** suggest to use R (>=4.0).

## Current version

* scCancer 2.2.1 (update at 2021.03.02)
* [All version log](https://github.com/wguo-research/scCancer/wiki/Version-Log)

## Installation

The detailed installation instruction can be found in the project [wiki]( https://github.com/wguo-research/scCancer/wiki/2.-Installation).


## Usage

The vignette of `scCancer` can be found in the project [wiki]( https://github.com/wguo-research/scCancer/wiki).

* [Quick start](https://github.com/wguo-research/scCancer/wiki/3.-Quick-start)
* [Step by step introduction](https://github.com/wguo-research/scCancer/wiki/4.-Step-by-step-introduction)
* [Other personalized settings](https://github.com/wguo-research/scCancer/wiki/5.-Other-personalized-settings)

We provide an [example data](http://lifeome.net/software/sccancer/KC-example.tar.gz) of kidney cancer from 10X Genomics, and following are the generated HTML reports:

* [`report-scStat.html`](http://lifeome.net/software/sccancer/KC-example-report-scStat.html)
* [`report-scAnno.html`](http://lifeome.net/software/sccancer/KC-example-report-scAnno.html)

For multi-datasets, following is a generated HTML report for three kidney cancer samples integration analysis:

* [`report-scAnnoComb.html`](http://lifeome.net/software/sccancer/KC123-report-scAnnoComb.html)


## Citation
Please use the following citation:

Wenbo Guo, Dongfang Wang, Shicheng Wang, Yiran Shan, Changyi Liu, Jin Gu, scCancer: a package for automated processing of single-cell RNA-seq data in cancer, _Briefings in Bioinformatics_,  bbaa127, [https://doi.org/10.1093/bib/bbaa127](https://doi.org/10.1093/bib/bbaa127)

## License
GPL-3


# scCancer2

## Introduction

We updated our R toolkit, scCancer, based on massive single-cell transcriptome and spatial transcriptome data.

1. Cell subtype annotation and cross-dataset label similarity: Our analysis mainly focused on cell subtype annotation by training multiple lightweight machine-learning models on scRNA-seq data. We proposed a method for quantitatively evaluating the similarity of cell subtype labels originating from different published datasets. We fully preserved the original labeling in cell atlases and analyzed the relationship between cell subtypes across datasets.

2. Malignant cell identification: We constructed a reference dataset combining scRNA-seq and bulk RNA-seq data across multiple cancer types to identify the malignant cell in TME. We trained a model to identify malignant cells with high generalization ability and computational efficiency. 

3. Spatial transcriptome analysis: Finally, we integrated a spatial transcriptome analysis pipeline. It enables us to analyze TME from a spatial perspective systematically and automatically.

With scCancer2, researchers can understand the composition of the TME more accurately from multiple dimensions.

## Overview of scCancer2
![image](https://github.com/czythu/scCancer_MicroEnv/blob/master/CellSubtypeAnnotation/vignettes/Figure1.png)

## System Requirements

R version: >= 3.5.0

## Package Installation and Quick Start for scRNA-seq analysis

Quick start of scCancer2:

1. Dependency installation

Some dependent packages for scCancer2 (old version of scCancer, edgeR, garnett, xgboost, and org.Hs.eg.db) may not be installed automatically, so you can install them from the following steps. After installing them successfully, then you can run the vignettes/scCancer2.rmd

```R
checkPkg <- function(pkg){
    return(requireNamespace(pkg, quietly = TRUE))
}

# Some frequently used packages 
if(!checkPkg("Seurat")) BiocManager::install("Seurat")
if(!checkPkg("knitr")) BiocManager::install("knitr")
if(!checkPkg("GSVA")) BiocManager::install("GSVA")
if(!checkPkg("pheatmap")) BiocManager::install("pheatmap")
if(!checkPkg("ComplexHeatmap")) BiocManager::install("ComplexHeatmap")
if(!checkPkg("markdown")) install.packages("markdown")
if(!checkPkg("R.utils")) install.packages("R.utils")
if(!checkPkg("RcppArmadillo")) install.packages("RcppArmadillo")
if(!checkPkg("RcppProgress")) install.packages("RcppProgress")

# Algorithms integrated in scCancer1.0
if(!checkPkg("harmony")) devtools::install_github("immunogenomics/harmony")
if(!checkPkg("NNLM")) devtools::install_github("linxihui/NNLM")
if(!checkPkg("liger")) devtools::install_github("MacoskoLab/liger")

# Algorithms integrated in scCancer2.0
if(!checkPkg("monocle")) BiocManager::install(c("monocle"))
if(!checkPkg("edgeR")) BiocManager::install(c("edgeR"))
if(!checkPkg(c('DelayedArray','DelayedMatrixStats','org.Hs.eg.db','org.Mm.eg.db'))) BiocManager::install(c('DelayedArray','DelayedMatrixStats','org.Hs.eg.db','org.Mm.eg.db'))
if(!checkPkg("garnett")) devtools::install_github("cole-trapnell-lab/garnett")
if(!checkPkg("xgboost")) install.packages("xgboost")
```

if errors occur when installing "NNLM" or "edgeR", you may install them from the .tar.gz file:

https://cran.r-project.org/src/contrib/Archive/NNLM/

https://bioconductor.org/packages/release/bioc/html/edgeR.html

2. If you have already installed the above dependencies, you have 2 ways to run scCancer2.0:

(a) Recommended: if you want to completely update scCancer to the next version:

```R
# install scCancer2.0
devtools::install_github("czythu/scCancer")
```

(b) Temporary installation

```R
# install scCancer1.0
# devtools::install_github("wguo-research/scCancer")
```

Download .zip file of scCancer2 package. Run vignettes/scCancer2.rmd in the scCancer folder. 

See vignettes/scCancer2.rmd for temporary installation and demos.

Recommended demo: CRC-example (Source: GSE146771)

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
```

## Citation
