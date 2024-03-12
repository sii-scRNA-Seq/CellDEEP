
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CellDEEP

<!-- badges: start -->
<!-- badges: end -->

The goal of CellDEEP is to …

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
pdc <- prepare_data(pDC, sample_id = "sample_id", group_id = "Worst_Clinical_Status", cluster_id = "initial_clustering", assay = "covid")

de.test <- FindMarker.CellDEEP(pdc, Pool = TRUE, ident.1 = "Healthy", ident.2 = "Severe", group.by = "group_id", n_cells = 6, cell_cutoff = 10, assay = "RNA")
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
