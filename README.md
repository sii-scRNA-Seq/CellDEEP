
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

<figure>
<img src="man/figures/CellDEEP_overview.png" alt="CellDEEP overview" />
<figcaption aria-hidden="true">CellDEEP overview</figcaption>
</figure>

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

### Vignette:

- [Tutorial
  vignette](https://sii-scrna-seq.github.io/CellDEEP/articles/CellDEEP_vignette.html)

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

de.test
```

The output here would be a gene list. If you want the whole DE output
list include p-value and logFC, simply change the parameter `name_only`
to FALSE:

``` r
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
                          ident.2 = "Control",
                          name.only = FALSE)
#> Start Pooling.....
#> Pooling...
#> Warning: Data is of class matrix. Coercing to dgCMatrix.
#> Pooling summary (random):
#> Input cells: 200
#> Cells kept in pooled pseudocells: 180
#> Cells not kept (approx): 20
#> Skipped empty groups: 0
#> Skipped empty clusters: 0
#> Skipped empty samples: 20
#> Skipped subgroups (<= min_cells_per_subgroup): 0
#> Dropped remainder cells (< n_cells) after random pooling: 20
#> FindMarker running.....
#> 1st ident is:
#> Case
#> 2nd ident is:
#> Control
#> group by:
#> group_id
#> Normalizing layer: counts
#> Finding variable features for layer counts
#> Centering and scaling data matrix
#> For a (much!) faster implementation of the Wilcoxon Rank Sum Test,
#> (default method for FindMarkers) please install the presto package
#> --------------------------------------------
#> install.packages('devtools')
#> devtools::install_github('immunogenomics/presto')
#> --------------------------------------------
#> After installation of presto, Seurat will automatically use the more 
#> efficient implementation (no further action necessary).
#> This message will be shown once per session

head(de.test)
#>                 p_val avg_log2FC pct.1 pct.2    p_val_adj
#> Gene1728 2.398196e-11   3.493846   1.0 0.567 4.796393e-08
#> Gene1992 5.772013e-11 -11.401061   0.0 0.900 1.154403e-07
#> Gene1626 8.774766e-11  -3.879532   0.2 1.000 1.754953e-07
#> Gene1864 1.933222e-10 -10.188828   0.0 0.867 3.866443e-07
#> Gene1715 1.934584e-10  -9.855602   0.0 0.867 3.869168e-07
#> Gene1807 6.247026e-10  -9.667764   0.0 0.833 1.249405e-06
```

And if you want full list without any filtering for analysis like GSEA,
simply set `full_list` to TRUE:

``` r
de.test.full <- FindMarker.CellDEEP(sim, 
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
                          ident.2 = "Control",
                          name.only = FALSE,
                          full_list = TRUE)
#> Start Pooling.....
#> Pooling...
#> Warning: Data is of class matrix. Coercing to dgCMatrix.
#> Pooling summary (random):
#> Input cells: 200
#> Cells kept in pooled pseudocells: 180
#> Cells not kept (approx): 20
#> Skipped empty groups: 0
#> Skipped empty clusters: 0
#> Skipped empty samples: 20
#> Skipped subgroups (<= min_cells_per_subgroup): 0
#> Dropped remainder cells (< n_cells) after random pooling: 20
#> FindMarker running.....
#> 1st ident is:
#> Case
#> 2nd ident is:
#> Control
#> group by:
#> group_id
#> Normalizing layer: counts
#> Finding variable features for layer counts
#> Centering and scaling data matrix

head(de.test.full)
#>                 p_val avg_log2FC pct.1 pct.2    p_val_adj
#> Gene1728 2.398196e-11   3.493846   1.0 0.567 4.796393e-08
#> Gene1992 5.772013e-11 -11.401061   0.0 0.900 1.154403e-07
#> Gene1626 8.774766e-11  -3.879532   0.2 1.000 1.754953e-07
#> Gene1864 1.933222e-10 -10.188828   0.0 0.867 3.866443e-07
#> Gene1715 1.934584e-10  -9.855602   0.0 0.867 3.869168e-07
#> Gene1807 6.247026e-10  -9.667764   0.0 0.833 1.249405e-06

#see how different they are
length(rownames(de.test))
#> [1] 20
length(rownames(de.test.full))
#> [1] 1054
```

## This section introduce what is updated.

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
