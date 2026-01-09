
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
#' "sum", "10X" (mean * 10), or "normalized.sum".
#' @param cell_cutoff Integer. Minimum cells required in a sample-cluster group
#' to perform pooling (default 25).
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

CellDEEP.Kmean <- function(dataset, n_cells= 10, nstart=100, assay_name="RNA", readcounts = "mean", cell_cutoff = 25){

  pseudo_cell_mtx <- matrix(, nrow=length(dataset[[assay_name]]@counts@Dimnames[[1]]), ncol=0)

  #Here: filter cluster(after all splitting) whose cell number < 25
  meta_data = c()
  group_id = c()
  cluster_id = c()
  sample_id = c()
  ktable = data.frame(row.names = rownames(dataset))

  print("Pooling...")
  for (x in levels(as.factor(dataset$group_id))){ # This is important
    group_subset <- subset(dataset, subset= group_id ==x)

    #for each group...
    for (z in levels(as.factor(dataset$cluster_id))) {
      cluster_subset <- subset(group_subset,subset=cluster_id==z)

      if(z %in% cluster_subset$cluster_id ){

        for(y in levels(as.factor(dataset$sample_id))){

          if(y %in% as.factor(cluster_subset$sample_id)){
            counter = 0
            sample_subset <- subset(cluster_subset, subset=sample_id==y)
            if(as.integer(length(colnames(sample_subset))) > cell_cutoff){
              k = as.integer(length(colnames(sample_subset))/n_cells)+1

              if (k<2){
                stop("Error: k < 2 for at least one sample. Change n_cells")
              }
              sample_subset@meta.data$kmeans <- kmeans(x = sample_subset@reductions[["pca"]]@cell.embeddings,centers = k, nstart = nstart)$cluster


              for(h in levels(as.factor(sample_subset@meta.data$kmeans))){


                k.clusters <- subset(sample_subset,subset=kmeans==h)

                #pool cells
                cells <- rownames(k.clusters@meta.data) #get cell rownames for kcluster


                cell_number <- length(cells) #get cell number that would be pooled
                pool <- k.clusters[[assay_name]]@counts[,cells] #get the cell information inside kcluster
                exp_mtx <- as.matrix(pool) #make a matrix of cell information inside kcluster
                sum_total <- rowSums(exp_mtx)

                #prepare read counts for normalized sum
                k.clusters[[assay_name]]@data <- as.matrix(k.clusters[[assay_name]]@data)
                pool.nor <- k.clusters[[assay_name]]@data[,cells] #get the cell information inside kcluster
                exp_mtx.nor <- as.matrix(pool.nor) #make a matrix of cell information inside kcluster
                sum_total.nor <- rowSums(exp_mtx.nor)

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
                }else if (readcounts == "normalized.sum") {
                  sum_total.nor <- data.frame(sum_total.nor)
                  pseudo_cell_mtx <- cbind(pseudo_cell_mtx, sum_total.nor$sum_total.nor)
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
              }
            }
          }
        }
      }
    }
  }
  #Create Seurat object
  row.names(pseudo_cell_mtx) <- dataset[[assay_name]]@counts@Dimnames[[1]]
  colnames(pseudo_cell_mtx) <- meta_data
  pseudo_cell_seurat <- Seurat::CreateSeuratObject(counts = pseudo_cell_mtx)
  pseudo_cell_seurat$group_id <- group_id
  pseudo_cell_seurat$sample_id <- sample_id
  pseudo_cell_seurat$cluster_id <- cluster_id


  dataset <- AddMetaData(dataset, metadata = ktable, col.name = "Pooled_kmeans_cells")
  Idents(dataset) <- "Pooled_kmeans_cells"

  table(pseudo_cell_seurat@meta.data$sample_id)
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
#' @param cell_cutoff Integer. Minimum cells required in a sample-cluster group
#' to perform pooling (default 25).
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
CellDEEP.Random <- function(dataset, n_cells= 10, assay_name="RNA", cell_cutoff = 25, readcounts = "mean"){
  pseudo_cell_mtx <- matrix(, nrow=length(dataset[[assay_name]]@counts@Dimnames[[1]]), ncol=0)

  meta_data = c()
  group_id = c()
  cluster_id = c()
  sample_id = c()
  rtable = data.frame(row.names = rownames(dataset))

  print("Pooling...")
  for (x in levels(as.factor(dataset$group_id))){
    print(x)
    group_subset <- subset(dataset, subset= group_id ==x)

    #for each group...
    for (z in levels(as.factor(group_subset$cluster_id))) {
      cluster_subset <- subset(group_subset,subset=cluster_id==z)

      #if the group has the cluster (?! Not sure why this step, to handle multiple clusters?)
      if(z %in% cluster_subset$cluster_id ){
        #print("IF Z...")

        #For each sample/patient/replicate
        for(y in levels(as.factor(dataset$sample_id))){
          counter = 0

          #If the patient belong to the same group z
          if(y %in% as.factor(cluster_subset$sample_id)){
            sample_subset <- subset(cluster_subset, subset= sample_id ==y)  # subsets according to the sample/replicate

            if(as.integer(length(colnames(sample_subset))) > cell_cutoff){

              real.cells <- rownames(sample_subset@meta.data) #get cell rowname to pool
              cluster_counter = 0

              while (length(real.cells) >=n_cells){  #when there are more than n cells in the cluster
                pool<- sample(real.cells, n_cells, replace = FALSE) #randomly pool n cells from the subset

                cell_id_to_pool <- pool

                real.cells <- subset(real.cells, !(real.cells %in% pool)) #delete those n cells from the subset
                pool <- cluster_subset[[assay_name]]@counts[,pool] #get the n cells readcounts(before was only names)
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
            }
          }
        }
      }
    }
  }
  #Create Seurat object
  row.names(pseudo_cell_mtx) <- dataset[[assay_name]]@counts@Dimnames[[1]]
  colnames(pseudo_cell_mtx) <- meta_data
  pseudo_cell_seurat <- CreateSeuratObject(counts = pseudo_cell_mtx)
  pseudo_cell_seurat$group_id <- group_id
  pseudo_cell_seurat$sample_id <- sample_id
  pseudo_cell_seurat$cluster_id <- cluster_id

  #split UMAP
  dataset <- AddMetaData(dataset, metadata = rtable, col.name = "Pooled_randomly_cells")
  Idents(dataset) <- "Pooled_randomly_cells"
  print(Seurat::DimPlot(dataset, split.by = "sample_id",ncol = 4))

  table(pseudo_cell_seurat@meta.data$sample_id)

  return(pseudo_cell_seurat)
}
