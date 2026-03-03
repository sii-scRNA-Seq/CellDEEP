
#' @title K-means Based Cell Pooling for Seurat Objects
#'
#' @rdname CellDEEP.Kmean
#'
#' @description
#' Pools cells into "pseudocells" by applying k-means clustering to PCA embeddings.
#' This reduces data sparsity while maintaining the biological grouping of
#' sample, cluster, and condition.
#'
#' @param dataset A Seurat object. Must have PCA reductions calculated.
#' @param n_cells Integer. Target number of cells to pool into each pseudocell.
#' @param nstart Integer. Number of random sets to start with in \code{kmeans}.
#' @param assay_name Character. The assay to pull counts from (default "RNA").
#' @param readcounts Character. Aggregation method: "mean" (rounded average),
#' "sum", "10X" (mean * 10).
#' @param min_cells_per_subgroup Integer. Minimum cells required in each
#' sample-cluster subgroup to perform pooling (default 25).
#'
#' @import Seurat
#'
#' @return A new Seurat object where each "cell" is a pooled group of original cells.
#'
#' @note
#' This function requires that PCA has already been run on the input \code{dataset},
#' as it uses the "pca" reduction for clustering.
#'
#' @export
#'
#' @examples
#' # pooled_obj <- CellDEEP.Kmean(dataset = my_seurat, n_cells = 10, readcounts = "mean")

CellDEEP.Kmean <- function(dataset, n_cells= 10, nstart=100, assay_name="RNA",
                           readcounts = "mean", min_cells_per_subgroup = 25){

  pseudo_cell_mtx <- matrix(, nrow=length(dataset[[assay_name]]$counts@Dimnames[[1]]), ncol=0)

  #Here: filter cluster(after all splitting) whose cell number < 25
  meta_data = c()
  group_id = c()
  cluster_id = c()
  sample_id = c()
  ktable = data.frame(row.names = rownames(dataset))
  drop_out_counter = 0
  total_input_cells <- ncol(dataset)
  empty_group_skips <- 0
  empty_cluster_skips <- 0
  empty_sample_skips <- 0
  below_min_subgroup_skips <- 0
  pooled_input_cells <- 0
  singleton_dropped_cells <- 0

  message("Pooling...")
  for (x in levels(as.factor(dataset$group_id))){ # This is important
    group_cells <- rownames(dataset@meta.data)[dataset@meta.data$group_id == x]
    if (length(group_cells) == 0) {
      empty_group_skips <- empty_group_skips + 1
      next
    }
    group_subset <- subset(dataset, cells = group_cells)

    #for each group...
    for (z in levels(as.factor(dataset$cluster_id))) {
      cluster_cells <- rownames(group_subset@meta.data)[group_subset@meta.data$cluster_id == z]
      if (length(cluster_cells) == 0) {
        empty_cluster_skips <- empty_cluster_skips + 1
        next
      }
      cluster_subset <- subset(group_subset, cells = cluster_cells)

        for(y in levels(as.factor(dataset$sample_id))){
          sample_cells <- rownames(cluster_subset@meta.data)[cluster_subset@meta.data$sample_id == y]
          if (length(sample_cells) == 0) {
            empty_sample_skips <- empty_sample_skips + 1
            next
          }
          counter = 0
          sample_subset <- subset(cluster_subset, cells = sample_cells)
          if(as.integer(length(colnames(sample_subset))) > min_cells_per_subgroup){

              k = as.integer(length(colnames(sample_subset))/n_cells)+1

              if (k<2){
                stop("Error: k < 2 for at least one sample. Change n_cells")
              }
              sample_subset@meta.data$kmeans <- kmeans(x = sample_subset@reductions[["pca"]]@cell.embeddings,centers = k, nstart = nstart)$cluster


              for(h in levels(as.factor(sample_subset@meta.data$kmeans))){


                k.clusters <- subset(sample_subset,subset=kmeans==h)

                #pool cells
                cells <- rownames(k.clusters@meta.data) #get cell rownames for kcluster
                cell_number <- length(cells)
                if(cell_number > 1){
                  pooled_input_cells <- pooled_input_cells + cell_number

                 #get cell number that would be pooled
                  pool <- k.clusters[[assay_name]]$counts[,cells] #get the cell information inside kcluster
                  exp_mtx <- as.matrix(pool) #make a matrix of cell information inside kcluster
                  sum_total <- rowSums(exp_mtx)


                  if (readcounts == "mean") {
                    #mean_total <- round(sum_total/n_cells)
                    mean_total <- round(sum_total/cell_number)
                    mean_total <- data.frame(mean_total) #make a dataframe
                    pseudo_cell_mtx <- cbind(pseudo_cell_mtx, mean_total$mean_total)
                  } else if (readcounts == "sum") {
                    sum_total <- data.frame(sum_total)
                    pseudo_cell_mtx <- cbind(pseudo_cell_mtx, sum_total$sum_total)
                  } else if (readcounts == "10X") {
                    #mean_total <- round(10*(sum_total/n_cells))
                    mean_total <- round(10*(sum_total/cell_number))
                    mean_total <- data.frame(mean_total) #make a dataframe
                    pseudo_cell_mtx <- cbind(pseudo_cell_mtx, mean_total$mean_total)
                  }
                  else {
                    stop("Error: readcounts parameter not known")
                  }

                  #increase counters
                  counter = counter + 1
                  meta_data <- append(meta_data, paste(y,"_",counter)) # New cell name
                  sample_id <- append(sample_id, paste(y))
                  group_id <- append(group_id, paste(x))
                  cluster_id <- append(cluster_id,paste(z))


                  ktable.cells <- data.frame(row.names = cells, pooled_cells=rep(paste(y,h,sep = "_"), length(cells)))
                  ktable <- rbind(ktable,ktable.cells)

                }else{
                  drop_out_counter = drop_out_counter + 1
                  singleton_dropped_cells <- singleton_dropped_cells + cell_number
                }
              }
          } else {
            below_min_subgroup_skips <- below_min_subgroup_skips + 1
          }
        }
    }
  }

  if (ncol(pseudo_cell_mtx) == 0) {
    stop("No pseudocells were generated. Check group/sample/cluster IDs or lower min_cells_per_subgroup.")
  }

  #Create Seurat object
  row.names(pseudo_cell_mtx) <- dataset[[assay_name]]$counts@Dimnames[[1]]
  colnames(pseudo_cell_mtx) <- meta_data
  pseudo_cell_seurat <- Seurat::CreateSeuratObject(counts = pseudo_cell_mtx)
  pseudo_cell_seurat$group_id <- group_id
  pseudo_cell_seurat$sample_id <- sample_id
  pseudo_cell_seurat$cluster_id <- cluster_id


  dataset <- AddMetaData(dataset, metadata = ktable, col.name = "Pooled_kmeans_cells")
  Idents(dataset) <- "Pooled_kmeans_cells"

  table(pseudo_cell_seurat@meta.data$sample_id)

  message("Drop out cell number during kmean pooling is:")
  message(drop_out_counter)
  message("Pooling summary (kmean):")
  message(paste0("Input cells: ", total_input_cells))
  message(paste0("Cells kept in pooled pseudocells: ", pooled_input_cells))
  message(paste0("Cells not kept (approx): ", total_input_cells - pooled_input_cells))
  message(paste0("Skipped empty groups: ", empty_group_skips))
  message(paste0("Skipped empty clusters: ", empty_cluster_skips))
  message(paste0("Skipped empty samples: ", empty_sample_skips))
  message(paste0("Skipped subgroups (<= min_cells_per_subgroup): ", below_min_subgroup_skips))
  message(paste0("Dropped singleton cells after kmeans split: ", singleton_dropped_cells))

  return(pseudo_cell_seurat)
}


#' @title Random Cell Pooling for Seurat Objects
#'
#' @rdname CellDEEP.Random
#'
#' @description
#' Pools cells into pseudocells by random selection within biological groups.
#' Includes a minimum threshold filter of 25 cells per subgroup to ensure
#' pooling quality.
#'
#' @param dataset A Seurat object.
#' @param n_cells Integer. The number of cells to pool into each pseudocell.
#' @param assay_name Character. The assay to use for counts (default "RNA").
#' @param readcounts Character. Method to aggregate counts: "sum" or "mean".
#' @param min_cells_per_subgroup Integer. Minimum cells required in each
#' sample-cluster subgroup to perform pooling (default 25).
#'
#' @import Seurat
#'
#' @return A new Seurat object containing the aggregated pseudocells.
#'
#' @note
#' Subgroups (sample-cluster combinations) with fewer than 25 cells are
#' automatically skipped. The function also generates a DimPlot to visualize
#' the random pooling across samples.
#'
#' @export
#'
#' @examples
#' # random_pooled_obj <- CellDEEP.Random(dataset = my_seurat, n_cells = 10, readcounts = "mean")
CellDEEP.Random <- function(dataset, n_cells= 10, assay_name="RNA",
                            min_cells_per_subgroup = 25, readcounts = "mean"){

  pseudo_cell_mtx <- matrix(, nrow=length(dataset[[assay_name]]$counts@Dimnames[[1]]), ncol=0)

  meta_data = c()
  group_id = c()
  cluster_id = c()
  sample_id = c()
  rtable = data.frame(row.names = rownames(dataset))
  total_input_cells <- ncol(dataset)
  empty_group_skips <- 0
  empty_cluster_skips <- 0
  empty_sample_skips <- 0
  below_min_subgroup_skips <- 0
  pooled_input_cells <- 0
  remainder_dropped_cells <- 0

  message("Pooling...")
  for (x in levels(as.factor(dataset$group_id))){
    group_cells <- rownames(dataset@meta.data)[dataset@meta.data$group_id == x]
    if (length(group_cells) == 0) {
      empty_group_skips <- empty_group_skips + 1
      next
    }
    group_subset <- subset(dataset, cells = group_cells)

    #for each group...
    for (z in levels(as.factor(group_subset$cluster_id))) {
      cluster_cells <- rownames(group_subset@meta.data)[group_subset@meta.data$cluster_id == z]
      if (length(cluster_cells) == 0) {
        empty_cluster_skips <- empty_cluster_skips + 1
        next
      }
      cluster_subset <- subset(group_subset, cells = cluster_cells)

        #For each sample/patient/replicate
        for(y in levels(as.factor(dataset$sample_id))){
          counter = 0

          sample_cells <- rownames(cluster_subset@meta.data)[cluster_subset@meta.data$sample_id == y]
          if (length(sample_cells) == 0) {
            empty_sample_skips <- empty_sample_skips + 1
            next
          }
          sample_subset <- subset(cluster_subset, cells = sample_cells)  # subsets according to the sample/replicate

          if(as.integer(length(colnames(sample_subset))) > min_cells_per_subgroup){

              real.cells <- rownames(sample_subset@meta.data) #get cell rowname to pool
              cluster_counter = 0

              while (length(real.cells) >=n_cells){  #when there are more than n cells in the cluster
                pool<- sample(real.cells, n_cells, replace = FALSE) #randomly pool n cells from the subset
                pooled_input_cells <- pooled_input_cells + n_cells

                cell_id_to_pool <- pool

                real.cells <- subset(real.cells, !(real.cells %in% pool)) #delete those n cells from the subset
                pool <- cluster_subset[[assay_name]]$counts[,pool] #get the n cells readcounts(before was only names)
                exp_mtx <- as.matrix(pool) #make a matrix and the mean
                sum_total <- rowSums(exp_mtx)

                if (readcounts == "mean") {
                  mean_total <- round(sum_total/n_cells)
                  mean_total <- data.frame(mean_total) #make a dataframe
                  pseudo_cell_mtx <- cbind(pseudo_cell_mtx, mean_total$mean_total)
                } else if (readcounts == "sum") {
                  sum_total <- data.frame(sum_total)
                  pseudo_cell_mtx <- cbind(pseudo_cell_mtx, sum_total$sum_total)
                } else {
                  stop("Error: readcounts parameter not known")
                }

                #increase counters
                counter = counter + 1
                cluster_counter= cluster_counter+1

                meta_data <- append(meta_data, paste(y,"_",counter)) # New cell name
                sample_id <- append(sample_id, paste(y))
                group_id <- append(group_id, paste(x))
                cluster_id <- append(cluster_id,paste(z))

                rtable.cells <- data.frame(row.names = cell_id_to_pool, pooled_cells=rep(paste(y,cluster_counter,sep = "_"), length(cell_id_to_pool)))
                rtable <- rbind(rtable,rtable.cells)
              }
              remainder_dropped_cells <- remainder_dropped_cells + length(real.cells)
          } else {
            below_min_subgroup_skips <- below_min_subgroup_skips + 1
          }
        }
    }
  }

  if (ncol(pseudo_cell_mtx) == 0) {
    stop("No pseudocells were generated. Check group/sample/cluster IDs or lower min_cells_per_subgroup.")
  }

  #Create Seurat object
  row.names(pseudo_cell_mtx) <- dataset[[assay_name]]$counts@Dimnames[[1]]
  colnames(pseudo_cell_mtx) <- meta_data
  pseudo_cell_seurat <- CreateSeuratObject(counts = pseudo_cell_mtx)
  pseudo_cell_seurat$group_id <- group_id
  pseudo_cell_seurat$sample_id <- sample_id
  pseudo_cell_seurat$cluster_id <- cluster_id

  #split UMAP
  dataset <- AddMetaData(dataset, metadata = rtable, col.name = "Pooled_randomly_cells")
  Idents(dataset) <- "Pooled_randomly_cells"

  table(pseudo_cell_seurat@meta.data$sample_id)
  message("Pooling summary (random):")
  message(paste0("Input cells: ", total_input_cells))
  message(paste0("Cells kept in pooled pseudocells: ", pooled_input_cells))
  message(paste0("Cells not kept (approx): ", total_input_cells - pooled_input_cells))
  message(paste0("Skipped empty groups: ", empty_group_skips))
  message(paste0("Skipped empty clusters: ", empty_cluster_skips))
  message(paste0("Skipped empty samples: ", empty_sample_skips))
  message(paste0("Skipped subgroups (<= min_cells_per_subgroup): ", below_min_subgroup_skips))
  message(paste0("Dropped remainder cells (< n_cells) after random pooling: ", remainder_dropped_cells))

  return(pseudo_cell_seurat)
}
