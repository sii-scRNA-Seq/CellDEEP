
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CellDEEP

<!-- badges: start -->

<!-- badges: end -->

## Introduction to CellDEEP

Single-cell RNA sequencing (scRNA-seq) allows us to explore gene
expression at an unprecedented resolution, but it faces significant
challenges like high dropout rates and data sparsity. CellDEEP was
developed to bridge the gap between robust but coarse pseudobulk methods
and sensitive but potentially biased single-cell methods.

Check our paper for details:
<https://gla-my.sharepoint.com/:w:/g/personal/domenico_somma_glasgow_ac_uk/IQClckLR7tE6SLF5Raw7J39PAWbMgDQf-qdc73jdZYNDUHE?e=j0le0G>

## Installation

You can install the development version of CellDEEP via github:

``` r
devtools::install_github("sii-scRNA-Seq/CellDEEP")
```

## Example

Before using, don’t forget library it:

``` r
library(CellDEEP)
```

Vignette:

``` r
vignette("CellDEEP_vignette", package = "CellDEEP")
```

To quickly run CellDEEP, pass your metadata column names directly into
`FindMarker.CellDEEP`:

``` r
data("sim")

# Pool defaults to TRUE
de.test <- FindMarker.CellDEEP(sim, 
                          group_id = "Status", 
                          sample_id = "DonorID", 
                          cluster_id = "cluster_id",
                          Pool = TRUE,
                          test.use = "wilcox", 
                          n_cells = 3, 
                          min_cells_per_subgroup = 1,
                          cell_selection = "random", 
                          readcounts = "sum", 
                          logfc.threshold = 0.25, 
                          ident.1 = "Case", 
                          ident.2 = "Control")
```

This section introduce what is updated.

For version 1.0.1:  
1. FindMarker.CellDEEP Pool default should be TRUE.  
2. Not clear what is cell_cutoff, replaced with new parameter  
3. Change “pool_way” to cell_selection  
4. Vignette easy to access/read  
5. Change toy data to simulated data, will generate DE result now.  

For publish version(1.0.0):  
Delete code used for experiment, keep only CellDEEP function code.  
Delete all the comments, clean the code.  
Rename pooling function as CellDEEP.Kmean and CellDEEP.Random.  
