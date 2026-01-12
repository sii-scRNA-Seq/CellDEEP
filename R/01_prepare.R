# Here is the function to help user prepare the data format CellDEEP required.
# Input shall be an seurat object, output shall be an object with *sample_id*, *group_id* and *cluster_id*.
# User input not only the ids it requires, but also does it need normalization/scale/UMAP generation
# Output object will be a seurat object with sample_id, group_id and cluster_id.

#' @title Standardize and Subset Seurat Objects
#'
#' @rdname prepare_data
#'
#' @description
#' Standardizes metadata column names to "sample_id", "group_id", and "cluster_id"
#' while handling optional data subsetting and ground truth extraction.
#'
#' @param Subset.Seurat A Seurat object.
#' @param assay Character. The assay to use (default "RNA").
#' @param sample_id Character. Metadata column name for samples.
#' @param group_id Character. Metadata column name for groups.
#' @param cluster_id Character. Metadata column name for clusters.
#' @param file_path Character. Path to save ground truth reports.
#'
#' @import Seurat
#'
#' @return A Seurat object with standardized metadata and updated factor levels.
#'
#' @export
#'
#' @examples
#' # Simple usage example:
#' # prepared_data <- prepare_data(
#' #   Subset.Seurat = my_obj,
#' #   sample_id = "orig.ident",
#' #   group_id = "condition",
#' #   cluster_id = "res.0.5"
#' # )
prepare_data <- function(Subset.Seurat, assay = "RNA",
                         sample_id, group_id, cluster_id,
                         file_path = NULL) {


  #Check metadata needed
  if (is.null(Subset.Seurat@meta.data[[cluster_id]]) == TRUE) {
    message("No cluster info, so assuming all cells belong to 1 cluster")
    Subset.Seurat$cluster_id <- "Cluster 0"
  }
  if (is.null(Subset.Seurat@meta.data[[group_id]]) == TRUE) {
    stop("group_id not present")
  }
  if (is.null(Subset.Seurat@meta.data[[sample_id]]) == TRUE) {
    stop("sample_id not present")
  }

  table(Subset.Seurat@meta.data[[sample_id]])
  table(Subset.Seurat@meta.data[[group_id]])
  table(Subset.Seurat@meta.data[[cluster_id]])

  Subset.Seurat@meta.data[["group_id"]] <- Subset.Seurat@meta.data[[group_id]]
  Subset.Seurat@meta.data[["sample_id"]] <- Subset.Seurat@meta.data[[sample_id]]
  Subset.Seurat@meta.data[["cluster_id"]] <- Subset.Seurat@meta.data[[cluster_id]]

  Subset.Seurat@meta.data[["group_id"]] <- droplevels(as.factor(Subset.Seurat@meta.data[["group_id"]]))
  Subset.Seurat@meta.data[["cluster_id"]] <- droplevels(as.factor(Subset.Seurat@meta.data[["cluster_id"]]))
  Subset.Seurat@meta.data[["sample_id"]] <- droplevels(as.factor(Subset.Seurat@meta.data[["sample_id"]]))

  return(Subset.Seurat)

}



