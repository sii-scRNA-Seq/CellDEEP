
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
                          test.use = "MAST", 
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
                          test.use = "MAST", 
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
#> 
#> Done!
#> Combining coefficients and standard errors
#> Calculating log-fold changes
#> Calculating likelihood ratio tests
#> Refitting on reduced model...
#> 
#> Done!

head(de.test)
#>                 p_val avg_log2FC pct.1 pct.2    p_val_adj
#> Gene1728 2.135401e-18   3.493846   1.0 0.567 4.270801e-15
#> Gene1992 5.520147e-15 -11.401061   0.0 0.900 1.104029e-11
#> Gene1864 5.149952e-14 -10.188828   0.0 0.867 1.029990e-10
#> Gene1715 5.149952e-14  -9.855602   0.0 0.867 1.029990e-10
#> Gene1807 3.904984e-13  -9.667764   0.0 0.833 7.809967e-10
#> Gene1626 7.039584e-13  -3.879532   0.2 1.000 1.407917e-09
```

And if you want full list without any filtering for analysis like GSEA,
simply set `full_list` to TRUE:

``` r
de.test.full <- FindMarker.CellDEEP(sim, 
                          group_id = "Status", 
                          sample_id = "DonorID", 
                          cluster_id = "cluster_id",
                          Pool = TRUE,
                          test.use = "MAST", 
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
#> 
#> Done!
#> Combining coefficients and standard errors
#> Calculating log-fold changes
#> Calculating likelihood ratio tests
#> Refitting on reduced model...
#> 
#> Done!

head(de.test.full)
#>                 p_val avg_log2FC pct.1 pct.2    p_val_adj
#> Gene1728 2.135401e-18   3.493846   1.0 0.567 4.270801e-15
#> Gene1992 5.520147e-15 -11.401061   0.0 0.900 1.104029e-11
#> Gene1864 5.149952e-14 -10.188828   0.0 0.867 1.029990e-10
#> Gene1715 5.149952e-14  -9.855602   0.0 0.867 1.029990e-10
#> Gene1807 3.904984e-13  -9.667764   0.0 0.833 7.809967e-10
#> Gene1626 7.039584e-13  -3.879532   0.2 1.000 1.407917e-09

#see how different they are
length(rownames(de.test))
#> [1] 27
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
