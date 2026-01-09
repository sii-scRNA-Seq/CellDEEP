# Here is the function to help user prepare the data format CellDEEP required.
# Input shall be an seurat object, output shall be an object with *sample_id*, *group_id* and *cluster_id*.
# User input not only the ids it requires, but also does it need normalization/scale/UMAP generation
# Output object will be a seurat object with sample_id, group_id and cluster_id.


#' @title Check DE genes against a ground-truth table
#' @rdname check.ground.truth
#' @description Compare a list of DE genes to a ground-truth reference and report overlaps.
#'
#' @param ground.truth.table Ground-truth object/list containing genes.
#' @param genes.DE Character vector of DE genes to check.
#' @param title Character; header label for output.
#' @param file_path Character; path to append results.
#' @param verbose Logical; print extra information.
#'
#' @return NULL (writes summary to file and prints to console).
#' @export
#'
#' @examples
#' # prepared_obj <- prepare_data(
#' #   Original.Seurat = my_obj,
#' #   sample_id = "orig.ident",
#' #   group_id = "condition",
#' #   cluster_id = "seurat_clusters"
#' # )
check.ground.truth <- function(ground.truth.table, genes.DE, title,file_path, verbose=FALSE) {

  title = paste0("### ", title, " ###")
  print(title)
  print(length(genes.DE))
  to_print = c()
  write(title, file_path, append = TRUE)
  write(length(genes.DE), file_path, append = TRUE)

  if (exists("Ground.truth")==TRUE){
    if(Ground.truth$type=="Muscat") {
      if (verbose) {
        print(head(ground.truth.table$genes))
        print(table(ground.truth.table$genes$category))
        cat('\n')
      }

      FP_EE=intersect(ground.truth.table$genes[ground.truth.table$genes$category=="ee", ]$gene, genes.DE)
      print(paste0("False Positive EE: ",length(FP_EE)))
      to_print <- append(to_print, length(FP_EE))

      FP_EP=intersect(ground.truth.table$genes[ground.truth.table$genes$category=="ep", ]$gene, genes.DE)
      print(paste0("False Positive EP: ",length(FP_EP)))
      to_print <- append(to_print, length(FP_EP))

      TP_DE=intersect(ground.truth.table$genes[ground.truth.table$genes$category=="de", ]$gene, genes.DE)
      print(paste0("True Positive DE: ", length(TP_DE), " on ", length(ground.truth.table$genes[ground.truth.table$genes$category=="de", ]$gene)))
      to_print <- append(to_print, length(TP_DE))

      TP_DP=intersect(ground.truth.table$genes[ground.truth.table$genes$category=="dp", ]$gene, genes.DE)
      print(paste0("True Positive DP: ", length(TP_DP), " on ", length(ground.truth.table$genes[ground.truth.table$genes$category=="dp", ]$gene)))
      to_print <- append(to_print, length(TP_DP))

      TP_DM=intersect(ground.truth.table$genes[ground.truth.table$genes$category=="dm", ]$gene, genes.DE)
      print(paste0("True Positive DM: ", length(TP_DM), " on ", length(ground.truth.table$genes[ground.truth.table$genes$category=="dm", ]$gene)))
      to_print <- append(to_print, length(TP_DM))

      TP_DB=intersect(ground.truth.table$genes[ground.truth.table$genes$category=="db", ]$gene, genes.DE)
      print(paste0("True Positive DB: ", length(TP_DB), " on ", length(ground.truth.table$genes[ground.truth.table$genes$category=="db", ]$gene)))
      to_print <- append(to_print, length(TP_DB))

    } else if (Ground.truth$type=="Other") {
      TP=intersect(ground.truth.table$genes, genes.DE)
      print(paste0("TRUE Positive: ",length(TP)))
      TP=paste0("TRUE Positive: ",length(TP))
      to_print <- append(to_print, TP)
    }
    write(to_print, file_path, append = TRUE, sep = "\t", ncolumns = 6)
    write("\n", file_path, append = TRUE)
  }
}

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
#' @param Ground.truth Logical. If TRUE, extracts ground truth genes.
#' @param file_path Character. Path to save ground truth reports.
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
                         Ground.truth = FALSE,
                         file_path = NULL) {


  #Check metadata needed
  if (is.null(Subset.Seurat@meta.data[[cluster_id]]) == TRUE) {
    print("No cluster info, so assuming all cells belong to 1 cluster")
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

  # --- Ground Truth Extraction Logic ---
  if (Ground.truth == TRUE) {
    # Warning: Check if file_path is provided
    if (is.null(file_path)) {
      stop("Argument 'file_path' is missing. A file path is required when Ground.truth is TRUE.")
    }
    if (!is.null(Subset.Seurat@misc$muscat_ground_truth)) {
      # Muscat ground truth
      GT_obj <- list() # Creating a local list to hold truth data
      GT_obj$type <- "Muscat"
      GT_obj$genes <- Subset.Seurat@misc$muscat_ground_truth
      write("Ground Truth:", file_path)
      write(names(table(GT_obj$genes$category)), file_path, append = TRUE, sep = "\t", ncolumns = 6)
      write(table(GT_obj$genes$category), file_path, append = TRUE, sep = "\t", ncolumns = 6)
      write("\n", file_path, append = TRUE)
    } else if (!is.null(Subset.Seurat@misc$DE.genes)) {
      # Other DE.genes as ground truth
      GT_obj <- list()
      GT_obj$type <- "Other"
      GT_obj$genes <- Subset.Seurat@misc$DE.genes
      write("Ground Truth:", file_path)
      write(length(GT_obj$genes), file_path, append = TRUE)
      write("\n", file_path, append = TRUE)

    } else {
      stop("Ground truth not found in @misc slot.")
    }
  }


  return(Subset.Seurat)

}



