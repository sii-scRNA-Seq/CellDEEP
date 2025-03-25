#' @title Subset Seurat object
#'
#' @rdname Subset.seurat.object
#'
#' @description This function take a seurat object as input, and subset it as user input.
#'
#' @param Seurat.object A Seurat object user want to subset
#' @param Ident.to.subset The ident to subset, eg: cluster_id
#' @param Specific.Idents.to.subset The ident want to susbet, eg: 5
#'
#' @import Seurat
#' @import SeuratObject
#'
#' @return Subset Seurat object
#' @export
#'
#' @examples
Subset.seurat.object <- function(Seurat.object, Ident.to.subset, Specific.Idents.to.subset) {
  Moment.Seurat.object = Seurat.object

  print(Moment.Seurat.object)
  for (i in 1:length(Ident.to.subset)) {
    Seurat::Idents(Moment.Seurat.object) <- Ident.to.subset[[i]]
    Moment.Seurat.object <- subset(Moment.Seurat.object, idents= Specific.Idents.to.subset[[i]])
    Moment.Seurat.object <- Seurat::ScaleData(Moment.Seurat.object)
  }
  print(Moment.Seurat.object)
  return(Moment.Seurat.object)
}


#' @title Adjust input data format, Normalize/Scale/UMAP if needed, make the object ready for further steps.
#'
#' @rdname prepare_data
#'
#' @description This function adjust input data format to make it ready for further steps. Add "sample_id",
#' "cluster_id" and "group_id". Can Normalize/Scale/RunUMAP if needed.
#'
#' @param obj A Seurat object, should have metadata for sample, cluster and group info.
#' @param assay RNA as default. If assay = "covid", create a new seurat obejct based on RAW read counts of input
#' data, and go with it.
#' @param sample_id The metadata has sample information.
#' @param group_id The metadata has group information.
#' @param cluster_id The metadata has cluster information.
#' @param Ident.to.subset The idents to subset(if user want to subset the data).
#' @param Need.to.Norm_Transf If the object need normalizarion or not.
#' @param Need.to.Scale If the object need scale or not
#' @param Need.UMAP.generation If the object need generate an UMAP or not.
#' @param Ground.truth Does this dataset have ground truth or not.
#' @param Specific.Idents.to.subset If need subset, which element do you want to choose.
#' @param file_path The file path to store ground truth
#'
#' @import Seurat
#' @import SeuratObject
#' @import data.table
#'
#' @returns
#' @export
#'
#' @examples
prepare_data <- function(obj, assay = "RNA",
                         sample_id, group_id, cluster_id,
                         Ident.to.subset = NULL,
                         Need.to.Norm_Transf = FALSE,
                         Need.to.Scale = FALSE,
                         Need.UMAP.generation = FALSE,
                         Ground.truth = FALSE,
                         Specific.Idents.to.subset = NULL,
                         file_path = NULL) {

  #get the data
  Original.Seurat <- obj

  ##For covid PBMC only, create another object with raw read counts only
  if (assay == "covid"){
    matrix <- as.matrix(Original.Seurat[["raw"]]@counts)
    Original.Seurat.test <- Seurat::CreateSeuratObject(matrix)
    Original.Seurat.test@meta.data <- Original.Seurat@meta.data
    Original.Seurat.test@reductions <- Original.Seurat@reductions
    Original.Seurat <- Original.Seurat.test
  }


  #Check metadata needed
  if (is.null(Original.Seurat@meta.data[[cluster_id]]) == TRUE) {
    print("No cluster info, so assuming all cells belong to 1 cluster")
    Original.Seurat$cluster_id <- "Cluster 0"
  }
  if (is.null(Original.Seurat@meta.data[[group_id]]) == TRUE) {
    stop("group_id not present")
  }
  if (is.null(Original.Seurat@meta.data[[sample_id]]) == TRUE) {
    stop("sample_id not present")
  }

  if (is.null(Ident.to.subset) == TRUE) {
    Subset.Seurat = Original.Seurat
  } else {
    Subset.Seurat = Subset.seurat.object(Original.Seurat, Ident.to.subset,Specific.Idents.to.subset)
    Subset.Seurat@meta.data[[group_id]] <- droplevels(as.factor(Subset.Seurat@meta.data[[group_id]]))
    Subset.Seurat@meta.data[[cluster_id]] <- droplevels(as.factor(Subset.Seurat@meta.data[[cluster_id]]))
    Subset.Seurat@meta.data[[sample_id]] <- droplevels(as.factor(Subset.Seurat@meta.data[[sample_id]]))
  }

  table(Subset.Seurat@meta.data[[sample_id]])
  table(Subset.Seurat@meta.data[[group_id]])
  table(Subset.Seurat@meta.data[[cluster_id]])

  if (Need.to.Norm_Transf == TRUE) {
    Subset.Seurat <- Seurat::NormalizeData(Subset.Seurat)
    Subset.Seurat <- Seurat::FindVariableFeatures(Subset.Seurat, selection.method = "vst", nfeatures = 2000)
    Subset.Seurat <- Seurat::ScaleData(Subset.Seurat, features = rownames(Subset.Seurat))
  }

  if (Need.to.Scale == TRUE) {
    Subset.Seurat <- Seurat::ScaleData(Subset.Seurat, features = rownames(Subset.Seurat))
  }

  if (Need.UMAP.generation == TRUE) {
    Subset.Seurat <- Seurat::RunPCA(object = Subset.Seurat, npcs = 30)
    Seurat::ElbowPlot(object = Subset.Seurat, ndims  = 30)
    Seurat::DimPlot(Subset.Seurat, reduction = "pca")
    Subset.Seurat <- Seurat::FindNeighbors(Subset.Seurat, dims = 1:20)
    Subset.Seurat <- Seurat::FindClusters(Subset.Seurat,resolution = 0.1)
    Subset.Seurat <- Seurat::RunUMAP(Subset.Seurat, dims = 1:20)
  }

  #Subset.Seurat <- subset(x = emma, subset = seurat_clusters == "2")

  # #SaveRDS
  # if (is.null(Ident.to.subset) == FALSE) {
  #   saveRDS(Subset.Seurat, paste(folder_middle_path, Subset.filename.tosave, sep = "/"))
  # }

  rm(Original.Seurat) #You don't need it and we save memory

  Subset.Seurat@meta.data[["group_id"]] <- Subset.Seurat@meta.data[[group_id]]
  Subset.Seurat@meta.data[["sample_id"]] <- Subset.Seurat@meta.data[[sample_id]]
  Subset.Seurat@meta.data[["cluster_id"]] <- Subset.Seurat@meta.data[[cluster_id]]

  # drop extra levels in Subset.Seurat object
  Subset.Seurat@meta.data[["group_id"]] <- droplevels(as.factor(Subset.Seurat@meta.data[["group_id"]]))
  Subset.Seurat@meta.data[["cluster_id"]] <- droplevels(as.factor(Subset.Seurat@meta.data[["cluster_id"]]))
  Subset.Seurat@meta.data[["sample_id"]] <- droplevels(as.factor(Subset.Seurat@meta.data[["sample_id"]]))

  #if the object contains ground.truth, get it
  if (exists("Ground.truth")==TRUE) {
    if (Ground.truth==TRUE){
      if (is.null(Subset.Seurat@misc$muscat_ground_truth)==FALSE) { #muscat ground truth
        Ground.truth$type="Muscat"
        Ground.truth$genes=Subset.Seurat@misc$muscat_ground_truth
        write("Ground Truth:",file_path)
        write(names(table(Ground.truth$genes$category)),file_path, append = TRUE, sep = "\t", ncolumns = 6)
        write(table(Ground.truth$genes$category),file_path, append = TRUE, sep = "\t", ncolumns = 6)
        write("\n",file_path, append = TRUE)
      } else if (is.null(Subset.Seurat@misc$DE.genes)==FALSE) { #other DE.genes as ground truth
        Ground.truth$type="Other"
        Ground.truth$genes=Subset.Seurat@misc$DE.genes
        write("Ground Truth:",file_path)
        write(length(Subset.Seurat@misc$DE.genes),file_path, append = TRUE)
        write("\n",file_path, append = TRUE)
      } else {
        stop("Ground truth not found")
      }
    }
  }

  return(Subset.Seurat)

}


#' @title Pooling cells by k-mean clustering.
#'
#' @rdname cellPooling.kmean.dev.yiyi
#'
#' @description This function pool cells by k-mean clustering.
#'
#'
#' @param dataset A Seurat object
#' @param n_cells number of cells to pool together
#' @param nstart the nstart in kmeans clustering, which represents how many sets to start with.
#' @param assay_name the assay to pool
#' @param readcounts "mean" or "sum" or "10X". How to treat read counts when pool it. "mean" is mean the read counts and round; "sum" is sum; "10X" is the mean of read counts and 10 times it
#' @param cell_cutoff Cell filtering number, After subset according to sample, group and cluster,if remaining cells lower than this number, these cells will be filtered out.Default 25.
#'
#' @import data.table
#' @import Seurat
#' @import SeuratObject
#' @import tidyverse
#'
#' @importFrom stats kmeans model.matrix
#'
#' @returns
#' @export
#'
#' @examples
cellPooling.kmean.dev.yiyi <- function(dataset, n_cells= 10, nstart=100, assay_name="RNA", readcounts = "mean", cell_cutoff = 25){
  print(dataset)
  print(assay_name)
  pseudo_cell_mtx <- matrix(, nrow=length(dataset[[assay_name]]$counts@Dimnames[[1]]), ncol=0)

  #Here: filter cluster(after all splitting) whose cell number < 25
  meta_data = c()
  group_id = c()
  cluster_id = c()
  sample_id = c()
  ktable = data.frame(row.names = rownames(dataset))
  drop_out_counter = 0

  print("Pooling...")
  for (x in levels(as.factor(dataset$group_id))){ # This is important
    group_subset <- subset(dataset, subset= group_id ==x)

    #for each group...
    for (z in levels(as.factor(dataset$cluster_id))) {
      cluster_subset <- subset(group_subset,subset=cluster_id==z)

      if(z %in% cluster_subset$cluster_id ){
        for(y in levels(as.factor(dataset$sample_id))){

          #If the patient belong to the same group z
          if(y %in% as.factor(cluster_subset$sample_id)){
            counter = 0
            sample_subset <- subset(cluster_subset, subset=sample_id==y)

            #only keep going when cell number > 25
            if(as.integer(length(colnames(sample_subset))) > cell_cutoff){

              k = as.integer(length(colnames(sample_subset))/n_cells)+1
              if (k<2){
                stop("Error: k < 2 for at least one sample. Change n_cells")
              }
              sample_subset@meta.data$kmeans <- stats::kmeans(x = sample_subset@reductions[["pca"]]@cell.embeddings,centers = k, nstart = nstart)$cluster

              for(h in levels(as.factor(sample_subset@meta.data$kmeans))){

                k.clusters <- subset(sample_subset,subset=kmeans==h)

                #pool cells
                cells <- rownames(k.clusters@meta.data) #get cell rownames for kcluster


                cell_number <- length(cells) #get cell number that would be pooled

                if(cell_number > 1){

                  pool <- k.clusters[[assay_name]]$counts[,cells] #get the cell information inside kcluster
                  exp_mtx <- as.matrix(pool)#make a matrix of cell information inside kcluster
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
                  counter = counter + 1
                  meta_data <- append(meta_data, paste(y,"_",counter)) # New cell name

                  sample_id <- append(sample_id, paste(y))
                  group_id <- append(group_id, paste(x))
                  cluster_id <- append(cluster_id,paste(z))


                  ktable.cells <- data.frame(row.names = cells, pooled_cells=rep(paste(y,h,sep = "_"), length(cells)))
                  ktable <- rbind(ktable,ktable.cells)

                }else{
                  drop_out_counter = drop_out_counter + 1
                }


              }
            }
          }
        }
      }
    }
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

  # print("sample_id before dimplot:")
  # print(table(dataset$sample_id))
  # print(Seurat::DimPlot(dataset, split.by = "sample_id",ncol = 4) + NoLegend())

  table(pseudo_cell_seurat@meta.data$sample_id)

  print("Drop out cell number during kmean pooling is:")
  print(drop_out_counter)

  #return(ktable)
  #return(dataset)
  return(pseudo_cell_seurat)
}


#' @title Pooling cells randomly.
#'
#' @rdname random.cellPooling.dev.yiyi
#'
#' @description This function pool cells randomly.
#'
#' @param dataset A seurat object
#' @param n_cells The number of cells to pool together.
#' @param assay_name The assay to pool based_on
#' @param readcounts "sum" or "mean". How to generate read counts when pooling cells together.
#' @param cell_cutoff Cell filtering number, After subset according to sample, group and cluster,if remaining cells lower than this number, these cells will be filtered out.Default 25.
#'
#' @import data.table
#' @import Seurat
#' @import SeuratObject
#' @import tidyverse
#'
#' @returns
#' @export
#'
#' @examples
random.cellPooling.dev.yiyi <- function(dataset, n_cells= 10, assay_name="RNA", readcounts = "mean", cell_cutoff = 25){
  pseudo_cell_mtx <- matrix(, nrow=length(dataset[[assay_name]]$counts@Dimnames[[1]]), ncol=0)

  meta_data = c()
  group_id = c()
  cluster_id = c()
  sample_id = c()
  rtable = data.frame(row.names = rownames(dataset))

  print("Pooling...")
  for (x in levels(as.factor(dataset$group_id))){
    ##group, A and B, by ident()
    print(x)
    group_subset <- subset(dataset, subset= group_id ==x)

    #for each group...
    for (z in levels(as.factor(group_subset$cluster_id))) {
      cluster_subset <- subset(group_subset,subset=cluster_id==z)

      #if the group has the cluster (?! Not sure why this step, to handle multiple clusters?)
      if(z %in% cluster_subset$cluster_id ){

        #For each sample/patient/replicate
        for(y in levels(as.factor(dataset$sample_id))){
          counter = 0

          #If the patient belong to the same group z
          if(y %in% as.factor(cluster_subset$sample_id)){
            #pool cells
            sample_subset <- subset(cluster_subset, subset= sample_id ==y)  # subsets according to the sample/replicate

            if(as.integer(length(colnames(sample_subset))) > cell_cutoff){

              real.cells <- rownames(sample_subset@meta.data) #get cell rowname to pool
              cluster_counter = 0
              if (length(real.cells) >=n_cells) {


                while (length(real.cells) >=n_cells){  #when there are more than n cells in the cluster
                  pool<- sample(real.cells, n_cells, replace = FALSE) #randomly pool n cells from the subset

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
              }else{ #when cell number is less than cell number we want to pool
                #create pseudo cell directly
                # object: cluster_subset
                # cell name inside one sample: real.cells
                cell_id_to_pool <- real.cells
                pool <- cluster_subset[[assay_name]]$counts[,real.cells] #get the n cells readcounts(before was only names)
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
  row.names(pseudo_cell_mtx) <- dataset[[assay_name]]$counts@Dimnames[[1]]
  colnames(pseudo_cell_mtx) <- meta_data
  pseudo_cell_seurat <- CreateSeuratObject(counts = pseudo_cell_mtx)
  pseudo_cell_seurat$group_id <- group_id
  pseudo_cell_seurat$sample_id <- sample_id
  pseudo_cell_seurat$cluster_id <- cluster_id

  #split UMAP
  dataset <- AddMetaData(dataset, metadata = rtable, col.name = "Pooled_randomly_cells")
  Idents(dataset) <- "Pooled_randomly_cells"

  # print("sample_id before dimplot:")
  # print(table(dataset$sample_id))
  # print("dataset object:")
  # print(colnames(dataset@meta.data))
  # print(Seurat::DimPlot(dataset, split.by = "sample_id",ncol = 4))

  table(pseudo_cell_seurat@meta.data$sample_id)

  return(pseudo_cell_seurat)
}


#' Pseudobulk Function
#' @description This function pseudobulk the input sce object, and run DE analysis as contrast indicate
#'
#' @param sce The SingleCellExperiment object want to be processed
#' @param contrasts The contrast matrix used to indicate comparison user want to make
#' @param verbose Show the input comparison matrix or not
#' @param ... Places for other further possible parameter
#'
#' @import muscat
#' @import data.table
#' @import tidyverse
#'
#' @importFrom SummarizedExperiment metadata metadata<-
#' @importFrom stats model.matrix
#' @importFrom SingleCellExperiment int_colData
#' @importFrom utils globalVariables
#'
#' @returns
#' @export
#'
#' @examples
pseudobulk_dge_pool <- function(sce, contrasts, verbose= FALSE, ...) {
  if(missing(sce) | missing(contrasts)) {
    stop('sce and contrasts are required arguments')
  }
  if(any(c('cluster_id', 'sample_id', 'group_id') %in% names(sce@colData) == FALSE)) {
    stop('SingleCellExperiment must contain colData with columns: cluster_id, sample_id, group_id')
  }

  # Make valid names for contrast matrix; check that the new names do not
  # generate duplicates
  ngrp <- length(unique(sce@colData$group_id))
  sce@colData$group_id <- as.factor(make.names(sce@colData$group_id))

  if(length(unique(sce@colData$group_id)) != ngrp) {
    stop('Duplicate group_id generated while making valid names.\nRename values in sce@colData$group_id to be valid names')
  }

  # Create one pseudobulk for each cluster_id. Each pseudobulk has counts
  # aggregrated within sample_id
  pb <- aggregateData(sce, assay= "counts", fun= "sum", by= c("cluster_id", "sample_id"))
  # add a matrix we need


  # Create a design matrix - make sure order is consistent with pseudobulk
  design <- unique(as.data.table(sce@colData[, c('group_id', 'sample_id')]))

  samples <- colnames(pb@assays@data[[1]])
  design <- design[match(samples, design$sample_id)]

  if(verbose) {
    cat('Design:\n')
    print(design)
    cat('\n')
  }

  # Prepare for DGE
  mm <- model.matrix(~ 0 + design$group_id)
  colnames(mm) <- sub('design$group_id', '', colnames(mm), fixed= TRUE)
  rownames(mm) <- design$sample_id

  if(verbose) {
    cat('Model matrix:\n')
    print(mm)
    cat('\n')
  }

  contrast_mat <- limma::makeContrasts(contrasts= contrasts, levels = mm)
  if(verbose) {
    cat('Contrast matrix:\n')
    print(contrast_mat)
    cat('\n')
  }

  # muscat requires `experiment_info` to be not NULL but I'm not sure it is
  # needed. Using `pb@metadata[['experiment_info']] <- NA` should also do
  pb@metadata[['experiment_info']] <- design

  ##muscat requires another n_cells parameter in pb
  .n_cells <- function(x) {
    y <- int_colData(x)$n_cells
    if (is.null(y)) return(NULL)
    if (length(metadata(x)$agg_pars$by) == 2)
      y <- as.matrix(data.frame(y, check.names = FALSE))
    return(as.table(y))
  }
  metadata(pb)$n_cells <- .n_cells(pb)

  # Run DGE for each contrast and cluster_id
  print(pb@metadata$n_cells)
  print(pb@metadata)
  dge <- pbDS(pb, design= mm, verbose= verbose, contrast= contrast_mat,min_cells = 1,...)

  # Collate results in a single table
  tab <- list()
  cntr <- names(dge$table)
  cls <- names(dge$table[[1]])
  for(x in cntr) {
    for(y in cls) {
      tab[[length(tab) + 1]] <- dge$table[[x]][[y]]
    }
  }
  tab <- rbindlist(tab)
  return(tab)
}

utils::globalVariables("metadata")

#' Function to run DE analysis: FindMakrer + P-value filtering
#' @description This function perform DE analysis with FindMarker
#'
#' @param dataset Input seurat object
#' @param test.use The DE analysis method
#' @param DE.ident.1 The first ident to compare
#' @param DE.ident.2 The second ident to compare
#' @param DE.group The metadata where comparison based on - eg: group_id
#' @param assay The assay to run DE, default "RNA"
#' @param p_cutoff The adjusted p_value cutoff, 0.05 as default
#'
#' @import Seurat
#' @importFrom utils head
#'
#' @returns
#' @export
#'
#' @examples
return.DE <- function(dataset, test.use = "wilcox", DE.ident.1, DE.ident.2, DE.group, assay="RNA", p_cutoff = 0.05) {
  markers <- Seurat::FindMarkers(dataset, ident.1 = DE.ident.1, ident.2 = DE.ident.2, group.by = DE.group, test.use = test.use)
  DE <- rownames(subset(markers, markers$p_val_adj < p_cutoff, na.rm = T))
  print(length(DE))
  print(head(DE))
  return(DE)
}



#' Funciton to run CellDEEP with DE analysis
#' @description This function run CellDEEP, pooling cells then run DE analysis
#'
#' @param object The input seurat object (cells user want to pool)
#' @param ident.1 The first ident to compare for DE
#' @param ident.2 The second ident to compare for DE
#' @param group.by The group comparison for DE (metadata name where the DE will be applied)
#' @param test.use The DE method to use
#' @param Pool Pool cells or not, default TRUE, when user want to run DE withour any pooling set to FASLE
#' @param readcounts The way to treat readcounts while pooling, options: "sum" and "mean" (For Kmean clustering pooling, one extra option "10X", which is multiple mean readcounts 10 times )
#' @param n_cells How many cells to pool together, this parameter affects method performance, check paper for details. 10 as defatult.
#' @param assay The assay to process. "RNA" as default.
#' @param cell_cutoff Cell filtering number, After subset according to sample, group and cluster,if remaining cells lower than this number, these cells will be filtered out.Default 25.
#' @param pool_way The pooling atrategy, Two option: "kmean" and "random".
#' @param pvalue_cutoff The adjusted p-value cutoff for DE analysis.
#' @param ... Places for other further possible parameter
#'
#' @import Seurat
#' @import SeuratObject
#' @import data.table
#'
#' @returns
#' @export
#'
#' @examples
FindMarker.CellDEEP <- function(object,
                                ident.1 = NULL,
                                ident.2 = NULL,
                                group.by = NULL,
                                test.use = "wilcox",
                                Pool = TRUE,
                                readcounts = "sum",
                                n_cells = 10,
                                assay = "RNA",
                                cell_cutoff = 25,
                                pool_way = "kmean",
                                pvalue_cutoff = 0.05,
                                ...
){

  #if Pool = FALSE, run FindMarker directly
  if (Pool == FALSE) {
    de.markers <- return.DE(object,
                            DE.ident.1 = ident.1,
                            DE.ident.2 = ident.2,
                            DE.group = group.by,
                            assay = assay,
                            test.use = test.use,
                            ...)
  }else #When pool = TRUE, run CellDEEP progress
  {

    print("Start Pooling.....")

    if (pool_way == "kmean") {
      #Pooling the cell
      pooled.object <- cellPooling.kmean.dev.yiyi(object, readcounts = readcounts, n_cells= n_cells, assay_name = assay, cell_cutoff = cell_cutoff )

      print("FindMarker running.....")

      print("1st ident is:")
      print(ident.1)
      print("2nd ident is:")
      print(ident.2)
      print("group by:")
      print(group.by)

      pooled.object <- Seurat::NormalizeData(pooled.object)
      pooled.object <- Seurat::FindVariableFeatures(pooled.object, selection.method = "vst", nfeatures = 2000)
      all.genes <- rownames(pooled.object)
      pooled.object <- Seurat::ScaleData(pooled.object, features = all.genes)

      #Run FindMarker and return DE gene list
      de.markers <- return.DE(pooled.object,
                              DE.ident.1 = ident.1,
                              DE.ident.2 = ident.2,
                              DE.group = group.by,
                              test.use = test.use,
                              p_cutoff = pvalue_cutoff,
                              ...)

    }else if(pool_way == "random"){
      #Pooling the cell
      pooled.object <- random.cellPooling.dev.yiyi(object, readcounts = readcounts, n_cells= n_cells, assay_name = assay, cell_cutoff = cell_cutoff )

      print("FindMarker running.....")

      print("1st ident is:")
      print(ident.1)
      print("2nd ident is:")
      print(ident.2)
      print("group by:")
      print(group.by)

      pooled.object <- Seurat::NormalizeData(pooled.object)
      pooled.object <- Seurat::FindVariableFeatures(pooled.object, selection.method = "vst", nfeatures = 2000)
      all.genes <- rownames(pooled.object)
      pooled.object <- Seurat::ScaleData(pooled.object, features = all.genes)

      #Run FindMarker and return DE gene list
      de.markers <- return.DE(pooled.object,
                              DE.ident.1 = ident.1,
                              DE.ident.2 = ident.2,
                              DE.group = group.by,
                              test.use = test.use,
                              p_cutoff = pvalue_cutoff,
                              ...)
    }else{
      print("Wrong pooling way, kmean or random")
    }

  }

  return(de.markers)
}


