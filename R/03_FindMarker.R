#' @title Perform Differential Expression and Filter Results
#'
#' @rdname return.DE
#'
#' @description
#' A wrapper for \code{Seurat::FindMarkers} that simplifies the extraction of
#' Differentially Expressed (DE) genes. It allows for automatic p-value
#' filtering and provides the option to return either gene names or a full results table.
#'
#' @param dataset A Seurat object.
#' @param test.use Character. Denotes which test to use. Default is "wilcox".
#' @param DE.ident.1 Identifier(s) for the first group of cells.
#' @param DE.ident.2 Identifier(s) for the second group of cells.
#' @param DE.group Character. The metadata column to group by.
#' @param assay Character. The assay to use (default "RNA").
#' @param p_cutoff Numeric. Adjusted p-value threshold for filtering (default 0.05).
#' @param name.only Logical. If TRUE (default), returns a vector of gene names; if FALSE, returns a data.frame.
#' @param logfc.threshold Numeric. Minimum log fold change required (default 0.1).
#' @param min.pct Numeric. Minimum percentage of cells expressing the gene (default 0.01).
#' @param full_list Logical. If TRUE, overrides p-value filtering and returns all results.
#' @param ... Extra parameters passed to \code{Seurat::FindMarkers}.
#'
#' @return A character vector of gene names or a data.frame of markers depending on \code{name.only}.
#'
#' @export
#'
#' @examples
#' # Get gene names only
#' # genes <- return.DE(my_seurat, DE.ident.1 = "Control", DE.group = "condition")
#'
#' # Get full table of significant markers
#' # marker_table <- return.DE(my_seurat, DE.ident.1 = "B_Cell", DE.group = "celltype", name.only = FALSE)
return.DE <- function(
    dataset,
    test.use       = "wilcox",
    DE.ident.1,
    DE.ident.2,
    DE.group,
    assay          = "RNA",
    p_cutoff       = 0.05,
    name.only      = TRUE,
    logfc.threshold = 0.1,
    min.pct        = 0.01,
    full_list      = FALSE,
    ...
) {

  markers <- Seurat::FindMarkers(
    dataset,
    ident.1         = DE.ident.1,
    ident.2         = DE.ident.2,
    group.by        = DE.group,
    test.use        = test.use,
    assay           = assay,
    logfc.threshold = logfc.threshold,
    min.pct         = min.pct,
    ...
  )

  ## If user wants full list, do NOT filter by p-value
  if (isTRUE(full_list)) {
    if (isTRUE(name.only)) {
      DE <- rownames(markers)
      print(length(DE))
      print(head(DE))
      return(DE)
    } else {
      return(markers)
    }
  }

  ## Otherwise, apply p-value cutoff (using adjusted p-values)
  if (!"p_val_adj" %in% colnames(markers)) {
    stop("Column 'p_val_adj' not found in markers; check FindMarkers output.")
  }

  keep <- !is.na(markers$p_val_adj) & markers$p_val_adj < p_cutoff

  if (isTRUE(name.only)) {
    DE <- rownames(markers)[keep]
    print(length(DE))
    print(head(DE))
    return(DE)
  } else {
    markers_filt <- markers[keep, , drop = FALSE]
    return(markers_filt)
  }
}



#' @title Differential Expression with Optional Cell Pooling
#'
#' @rdname FindMarker.CellDEEP
#'
#' @description
#' A high-level wrapper that performs Differential Expression (DE) analysis.
#' It can run standard Seurat DE or first aggregate cells into pseudocells
#' using K-means or Random pooling before performing the analysis.
#'
#' @param object A Seurat object.
#' @param ident.1 Character. The first identity group to compare.
#' @param ident.2 Character. The second identity group to compare (optional).
#' @param group.by Character. The metadata column used for grouping cells.
#' @param test.use Character. DE test to use (default "wilcox").
#' @param Pool Logical. Whether to perform CellDEEP pooling before DE (default FALSE).
#' @param readcounts Character. Aggregation method for pools: "sum", "mean", or "10X".
#' @param n_cells Integer. Target number of cells per pool.
#' @param assay Character. Assay to use (default "RNA").
#' @param cell_cutoff Integer. Minimum cells per sample/cluster for pooling logic.
#' @param pool_way Character. Choice of "kmean" or "random" pooling.
#' @param name.only Logical. If TRUE, returns gene names only; otherwise, returns a data.frame.
#' @param logfc.threshold Numeric. Minimum log fold-change for DE filtering.
#' @param min.pct Numeric. Minimum detection rate for genes to be considered.
#' @param p_cutoff Numeric. P-value threshold for filtering (adjusted p-value).
#' @param full_list Logical. If TRUE, returns all genes regardless of p-value.
#' @param ... Additional arguments passed to \code{Seurat::FindMarkers}.
#'
#' @return A vector of gene names or a data.frame containing DE statistics.
#'
#' @export
#'
#' @examples
#' # Run standard Wilcoxon DE
#' # res <- FindMarker.CellDEEP(my_obj, ident.1 = "A", group.by = "group")
#'
#' # Run DE with K-means pooling
#' # res <- FindMarker.CellDEEP(my_obj, ident.1 = "A", group.by = "group", Pool = TRUE, pool_way = "kmean")
FindMarker.CellDEEP <- function(
    object,
    ident.1         = NULL,
    ident.2         = NULL,
    group.by        = NULL,
    test.use        = "wilcox",
    Pool            = FALSE,
    readcounts      = "sum",
    n_cells         = 10,
    assay           = "RNA",
    cell_cutoff     = 25,
    pool_way        = "kmean",
    name.only       = TRUE,
    logfc.threshold = 0.1,
    min.pct         = 0.01,
    p_cutoff        = 0.05,
    full_list       = FALSE,
    ...
) {

  # no pooling: run directly
  if (!isTRUE(Pool)) {
    de.markers <- return.DE(
      dataset        = object,
      DE.ident.1     = ident.1,
      DE.ident.2     = ident.2,
      DE.group       = group.by,
      assay          = assay,
      test.use       = test.use,
      name.only      = name.only,
      logfc.threshold = logfc.threshold,
      min.pct        = min.pct,
      p_cutoff       = p_cutoff,
      full_list      = full_list,
      ...
    )

  } else {  # Pool = TRUE : run CellDEEP

    print("Start Pooling.....")

    if (pool_way == "kmean") {

      # Pooling the cells
      pooled.object <- CellDEEP.Kmean(
        object,
        readcounts = readcounts,
        n_cells    = n_cells,
        assay_name = assay,
        cell_cutoff = cell_cutoff
      )

    } else if (pool_way == "random") {

      # Pooling the cells
      pooled.object <- CellDEEP.Random (
        object,
        readcounts = readcounts,
        n_cells    = n_cells,
        assay_name = assay
      )

    } else {
      stop("Wrong pooling way, choose 'kmean' or 'random'")
    }

    print("FindMarker running.....")
    print("1st ident is:"); print(ident.1)
    print("2nd ident is:"); print(ident.2)
    print("group by:");    print(group.by)

    pooled.object <- Seurat::NormalizeData(pooled.object)
    pooled.object <- Seurat::FindVariableFeatures(
      pooled.object,
      selection.method = "vst",
      nfeatures        = 2000
    )
    all.genes <- rownames(pooled.object)
    pooled.object <- Seurat::ScaleData(pooled.object, features = all.genes)

    # Run FindMarker and return DE gene list
    de.markers <- return.DE(
      dataset        = pooled.object,
      DE.ident.1     = ident.1,
      DE.ident.2     = ident.2,
      DE.group       = group.by,
      assay          = assay,
      test.use       = test.use,
      name.only      = name.only,
      logfc.threshold = logfc.threshold,
      min.pct        = min.pct,
      p_cutoff       = p_cutoff,
      full_list      = full_list,
      ...
    )
  }

  return(de.markers)
}
