
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

You can install the development version of CellDEEP like so:

Please download the package from github directly, then load the Rproject
file.

``` r
#This is not working for now
devtools::install_github("sii-scRNA-Seq/CellDEEP")
```

## Example

Before using, don’t forget library it:

``` r
library(CellDEEP)
```

To quickly run CellDEEP, just need to prepare the data, then run
`FindMarker.CellDEEP`:

``` r
#This is a test to see if the code is working, will reutrn 0 DE genes.

data("pDC")
pdc <- prepare_data(pDC, sample_id = "sample_id", group_id = "Worst_Clinical_Status", cluster_id = "initial_clustering")

#K-mean pooling:
de.test <- FindMarker.CellDEEP(pdc, Pool = TRUE, ident.1 = "Healthy", ident.2 = "Severe", group.by = "group_id", n_cells = 4, cell_cutoff = 4, assay = "RNA", pool_way = "kmean")
#Randomly pooling:
de.test <- FindMarker.CellDEEP(pdc, Pool = TRUE, ident.1 = "Healthy", ident.2 = "Severe", group.by = "group_id", n_cells = 4, cell_cutoff = 4, assay = "RNA", pool_way = "random")
```

## Update

This section introduce what is updated.

For this publish version: Delete code used for experiment, keep only
CellDEEP function code. Delete all the comments, clean the code. Rename
pooling function as CellDEEP.Kmean and CellDEEP.Random.

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
