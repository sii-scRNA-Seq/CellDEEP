#' @title Perform Differential Expression and Filter Results
#'
#' @rdname return.DE
#'
#' @description
#' A wrapper for \code{Seurat::FindMarkers} that simplifies the extraction of
#' Differentially Expressed (DE) genes. It supports p-value filtering and can
#' return either gene names or a full results table.
#'
#' @param dataset A Seurat object.
#' @param test.use Character. DE test to use (default \code{"wilcox"}).
#' @param DE.ident.1 Identifier(s) for the first group of cells.
#' @param DE.ident.2 Identifier(s) for the second group of cells.
#' @param DE.group Character. Metadata column to group by.
#' @param assay Character. Assay to use (default \code{"RNA"}).
#' @param p_cutoff Numeric. Adjusted p-value threshold (default 0.05).
#' @param name.only Logical. If TRUE, return gene names only.
#' @param logfc.threshold Numeric. Minimum log fold change (default 0.1).
#' @param min.pct Numeric. Minimum fraction of cells expressing a gene.
#' @param full_list Logical. If TRUE, return all genes and skip p-value filter.
#' @param ... Extra arguments passed to \code{Seurat::FindMarkers}.
#'
#' @importFrom Seurat FindMarkers
#' @return A character vector of genes or a marker data.frame.
#' @export
return.DE <- function(
    dataset,
    test.use        = "wilcox",
    DE.ident.1,
    DE.ident.2,
    DE.group,
    assay           = "RNA",
    p_cutoff        = 0.05,
    name.only       = TRUE,
    logfc.threshold = 0.25,
    min.pct         = 0.01,
    full_list       = FALSE,
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

  if (isTRUE(full_list)) {
    if (isTRUE(name.only)) {
      DE <- rownames(markers)
      message(length(DE))
      message(head(DE))
      return(DE)
    }
    return(markers)
  }

  if (!"p_val_adj" %in% colnames(markers)) {
    stop("Column 'p_val_adj' not found in markers; check FindMarkers output.")
  }

  keep <- !is.na(markers$p_val_adj) & markers$p_val_adj < p_cutoff

  if (isTRUE(name.only)) {
    DE <- rownames(markers)[keep]
    message(length(DE))
    message(head(DE))
    return(DE)
  }

  markers[keep, , drop = FALSE]
}


#' @title Standardize Seurat Metadata for CellDEEP
#'
#' @rdname prepare_data
#'
#' @description
#' Standardizes metadata columns to \code{sample_id}, \code{group_id}, and
#' \code{cluster_id} so CellDEEP functions can run consistently.
#'
#' @param Subset.Seurat A Seurat object.
#' @param assay Character. Assay to use (default \code{"RNA"}).
#' @param sample_id Character. Metadata column name for sample IDs.
#' @param group_id Character. Metadata column name for group IDs.
#' @param cluster_id Character. Metadata column name for cluster IDs.
#' @param file_path Character. Reserved for compatibility.
#'
#' @import Seurat
#' @return A Seurat object with standardized metadata fields.
#' @export
prepare_data <- function(Subset.Seurat, assay = "RNA",
                         sample_id, group_id, cluster_id,
                         file_path = NULL) {
  if (is.null(Subset.Seurat@meta.data[[group_id]])) {
    stop("group_id not present")
  }
  if (is.null(Subset.Seurat@meta.data[[sample_id]])) {
    stop("sample_id not present")
  }

  if (is.null(Subset.Seurat@meta.data[[cluster_id]])) {
    message("No cluster info, so assuming all cells belong to 1 cluster")
    Subset.Seurat@meta.data[["cluster_id"]] <- "Cluster 0"
  } else {
    Subset.Seurat@meta.data[["cluster_id"]] <- Subset.Seurat@meta.data[[cluster_id]]
  }

  Subset.Seurat@meta.data[["group_id"]] <- Subset.Seurat@meta.data[[group_id]]
  Subset.Seurat@meta.data[["sample_id"]] <- Subset.Seurat@meta.data[[sample_id]]

  Subset.Seurat@meta.data[["group_id"]] <- droplevels(as.factor(Subset.Seurat@meta.data[["group_id"]]))
  Subset.Seurat@meta.data[["cluster_id"]] <- droplevels(as.factor(Subset.Seurat@meta.data[["cluster_id"]]))
  Subset.Seurat@meta.data[["sample_id"]] <- droplevels(as.factor(Subset.Seurat@meta.data[["sample_id"]]))

  Subset.Seurat
}


#' @title Differential Expression with Optional Cell Pooling
#'
#' @rdname FindMarker.CellDEEP
#'
#' @description
#' It can run Seurat DE directly or first aggregate cells
#' into metacells using CellDEEP pooling.
#'
#' @param object A Seurat object.
#' @param ident.1 Character. First identity group to compare.
#' @param ident.2 Character. Second identity group to compare.
#' @param group.by Character. Metadata column used for grouping (default \code{"group_id"}).
#' @param sample_id Character. Input metadata column for sample IDs.
#' @param group_id Character. Input metadata column for group IDs.
#' @param cluster_id Character. Input metadata column for cluster IDs.
#' @param prepare Logical. If TRUE, run \code{prepare_data} first.
#' @param test.use Character. DE test to use.
#' @param Pool Logical. If TRUE, perform CellDEEP pooling before DE (default TRUE).
#' @param readcounts Character. Pool aggregation method: \code{"sum"}, \code{"mean"}, or \code{"10X"}.
#' @param n_cells Integer. Target number of cells per pool.
#' @param assay Character. Assay to use (default \code{"RNA"}).
#' @param min_cells_per_subgroup Integer. Minimum cells in each sample-cluster subgroup required for pooling.
#' @param cell_selection Character. Pooling strategy: \code{"kmean"} or \code{"random"}.
#' @param name.only Logical. If TRUE, return gene names only.
#' @param logfc.threshold Numeric. Minimum log fold-change.
#' @param min.pct Numeric. Minimum detection rate.
#' @param p_cutoff Numeric. Adjusted p-value threshold.
#' @param full_list Logical. If TRUE, return all genes regardless of p-value.
#' @param ... Additional arguments passed to \code{Seurat::FindMarkers}.
#'
#' @import Seurat
#' @return A vector of gene names or a DE data.frame.
#' @export
FindMarker.CellDEEP <- function(
    object,
    ident.1                 = NULL,
    ident.2                 = NULL,
    group.by                = "group_id",
    sample_id               = NULL,
    group_id                = NULL,
    cluster_id              = NULL,
    prepare                 = TRUE,
    test.use                = "wilcox",
    Pool                    = TRUE,
    readcounts              = "sum",
    n_cells                 = 10,
    assay                   = "RNA",
    min_cells_per_subgroup  = 25,
    cell_selection          = "kmean",
    name.only               = TRUE,
    logfc.threshold         = 0.25,
    min.pct                 = 0.01,
    p_cutoff                = 0.05,
    full_list               = FALSE,
    ...
) {
  if (isTRUE(prepare)) {
    if (is.null(sample_id) || is.null(group_id) || is.null(cluster_id)) {
      stop("When prepare = TRUE, provide sample_id, group_id, and cluster_id.")
    }

    object <- prepare_data(
      Subset.Seurat = object,
      assay = assay,
      sample_id = sample_id,
      group_id = group_id,
      cluster_id = cluster_id
    )
  }

  if (!isTRUE(Pool)) {
    return(return.DE(
      dataset         = object,
      DE.ident.1      = ident.1,
      DE.ident.2      = ident.2,
      DE.group        = group.by,
      assay           = assay,
      test.use        = test.use,
      name.only       = name.only,
      logfc.threshold = logfc.threshold,
      min.pct         = min.pct,
      p_cutoff        = p_cutoff,
      full_list       = full_list,
      ...
    ))
  }

  message("Start Pooling.....")

  if (cell_selection == "kmean") {
    pooled.object <- CellDEEP.Kmean(
      object,
      readcounts              = readcounts,
      n_cells                 = n_cells,
      assay_name              = assay,
      min_cells_per_subgroup  = min_cells_per_subgroup
    )
  } else if (cell_selection == "random") {
    pooled.object <- CellDEEP.Random(
      object,
      readcounts              = readcounts,
      n_cells                 = n_cells,
      assay_name              = assay,
      min_cells_per_subgroup  = min_cells_per_subgroup
    )
  } else {
    stop("Wrong cell selection method, choose 'kmean' or 'random'")
  }

  message("FindMarker running.....")
  message("1st ident is:"); message(ident.1)
  message("2nd ident is:"); message(ident.2)
  message("group by:"); message(group.by)

  pooled.object <- Seurat::NormalizeData(pooled.object)
  pooled.object <- Seurat::FindVariableFeatures(
    pooled.object,
    selection.method = "vst",
    nfeatures = 2000
  )
  pooled.object <- Seurat::ScaleData(pooled.object, features = rownames(pooled.object))

  return.DE(
    dataset         = pooled.object,
    DE.ident.1      = ident.1,
    DE.ident.2      = ident.2,
    DE.group        = group.by,
    assay           = assay,
    test.use        = test.use,
    name.only       = name.only,
    logfc.threshold = logfc.threshold,
    min.pct         = min.pct,
    p_cutoff        = p_cutoff,
    full_list       = full_list,
    ...
  )
}
