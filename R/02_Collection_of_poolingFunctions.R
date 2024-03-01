library(unixtools)
set.tempdir("/data/Dom/temporary")
ulimit::memory_limit(300000)

##------ Testing ---------

# # Uncomment to Read your dataset and test it

# Real.Seurat.subset <- readRDS("/data/Dom/datathon/Macrophage_cluster0_onlyHealhtyResistant.rds")
# Real.Seurat.subset
# DimPlot(Real.Seurat.subset)
# DimPlot(Real.Seurat.subset, split.by = "group_id")
#
# #test <- cellPooling.kmean.dev.dom2(Real.Seurat.subset, n_cells = 10, readcounts = "sum")
# test <- nn_pseudobulk(Real.Seurat.subset)
# test
# table(test[["original"]]@meta.data$sample_id)
# table(test[["pseudobulk"]]@meta.data$sample_id)
# table(test[["original"]]@meta.data$group_id)
# table(test[["pseudobulk"]]@meta.data$group_id)
# table(test[["original"]]@meta.data$cluster_id)
# table(test[["pseudobulk"]]@meta.data$cluster_id)
#
# Idents(test[["original"]]) <- "sample_id"
# only1sample <- subset(test[["original"]], idents = "SA144")
# DimPlot(only1sample, split.by = "bulk")

##------cell pooling functions ---------

#levels: replicate_id = patient id etc, condition_id = condition, cluster_id = cluster
nn_pseudobulk <- function(dataset, assay_name="RNA", min_virtual_cells=10, round_mtx=FALSE, verbose=FALSE){
  exp_mtx <- matrix(,nrow = length(row.names(dataset@assays[[assay_name]])), ncol = 0)
  row.names(exp_mtx) <- dataset[[assay_name]]@counts@Dimnames[[1]]
  meta_data <- data.frame(group_label = "i", sample_label = "a", cluster_label = "c", cells_in_bulk = 0)
  cells_meta_data <- data.frame(bulk = rep("c", length(colnames(dataset))), row.names = colnames(dataset))
  counter <- 1
  for (group in unique(dataset@meta.data$group_id)){
    group_subset <- subset(dataset, group_id == as.character(group))
    #print(group)

    for (sample in unique(group_subset@meta.data$sample_id)){
      # print(sample)
      sample_subset <- subset(group_subset, sample_id == as.character(sample))

      for (cluster in unique(sample_subset@meta.data$cluster_id)){
        #print(paste0("cluster: ",cluster))
        cluster_subset <- subset(sample_subset, cluster_id == as.character(cluster))
        res <- 0.7

        cluster_subset <- FindNeighbors(cluster_subset ,verbose = F)
        cluster_subset <- FindClusters(cluster_subset,resolution = res,verbose = F)

        while (length(unique(cluster_subset@active.ident)) < min_virtual_cells){

          if (verbose) {
            print(DefaultAssay(cluster_subset))
            print(cluster_subset)
          }

          cluster_subset <- FindNeighbors(cluster_subset, verbose = F)

          if (verbose) {
            print(DefaultAssay(cluster_subset))
            print(cluster_subset)
          }

          cluster_subset <- FindClusters(cluster_subset,resolution = res, verbose = F)
          res <- res + 0.1
        }

        #print("end")
        #print(levels(cluster_subset@active.ident))

        for (sub_cluster in unique(cluster_subset@active.ident)){
          #print(paste0("sub_cluster: ",sub_cluster))
          cells <- colnames(cluster_subset)[which(cluster_subset@active.ident == sub_cluster)]

          cells_meta_data[cells,1] <- paste0("c", counter)

          subset_mtx <- cluster_subset@assays[[assay_name]]@counts[, cells]

          pseudobulk <- rowMeans(subset_mtx)

          pseudobulk <- data.frame(pseudobulk,row.names = rownames(exp_mtx))

          exp_mtx <- cbind(exp_mtx,pseudobulk$pseudobulk)

          meta_data <- rbind(meta_data, data.frame(group_label = group, sample_label = sample,cluster_label = cluster, cells_in_bulk = length(cells)))

          counter <- counter + 1
        }

      }

    }

  }

  if (round_mtx == TRUE){
    exp_mtx <- round(exp_mtx)
    #print("TRUE")
  } else if (round_mtx == FALSE){
    #don't do anything, just use the round_mtx as it is
  }
  else if (round_mtx == "10X"){
    exp_mtx <- round(10*exp_mtx)
    #print("10X")
  } else {
    stop("Error: readcounts parameter not known")
  }

  #row.names(exp_mtx) <- row.names(dataset)
  meta_data <- meta_data[-1,]
  #print(meta_data)
  colnames(exp_mtx) <- paste0("c",seq(1, length(meta_data[,1]),1))
  row.names(meta_data) <- paste0("c",seq(1, length(meta_data[,1]),1))
  output_obj <- CreateSeuratObject(exp_mtx, project ="pseudobulk", meta.data= as.data.frame(meta_data))
  dataset <- AddMetaData(dataset, cells_meta_data)
  output_obj$sample_id <- output_obj$sample_label
  output_obj$sample_label <- NULL
  output_obj$cluster_id <- output_obj$cluster_label
  output_obj$cluster_label <- NULL
  output_obj$group_id <- output_obj$group_label
  output_obj$group_label <- NULL
  output_list <- list(output_obj, dataset)
  names(output_list) <- c("pseudobulk", "original")

  return(output_list)

}



cellPooling.kmean.dev.yiyi <- function(dataset, n_cells= 10, nstart=100, assay_name="RNA", readcounts = "mean"){

  pseudo_cell_mtx <- matrix(, nrow=length(dataset[[assay_name]]@counts@Dimnames[[1]]), ncol=0)

  #Here: filter cluster(after all splitting) whose cell number < 25
  meta_data = c()
  group_id = c()
  cluster_id = c()
  sample_id = c()
  ktable = data.frame(row.names = rownames(dataset))

  print("Pooling...")
  for (x in levels(as.factor(dataset$group_id))){ # This is important
    ##group, A and B, by ident()
    #print("here is group_id:")
    #print(x)
    group_subset <- subset(dataset, subset= group_id ==x)

    #for each group...
    for (z in levels(as.factor(dataset$cluster_id))) {
      #print("here is z, aka cluster_id(after group_id)")
      #print(z)
      #cluster_counter = 0
      cluster_subset <- subset(group_subset,subset=cluster_id==z)

      #if the group has the cluster (?! Not sure why this step, to handle multiple clusters?)
      if(z %in% cluster_subset$cluster_id ){
        #print("IF Z...")

        #For each sample/patient/replicate
        for(y in levels(as.factor(dataset$sample_id))){
          # a <- strsplit(y,".",-1)
          # a <- a[[1]][2]
          # if(a==x){
          #print("here is y, aka sample_id:")
          #print(y)
          #counter = 0

          #If the patient belong to the same group z
          if(y %in% as.factor(cluster_subset$sample_id)){
            #print("IF Y...")

            #if(length(group.cluster$sample_id) > 10){ # DOM: should be n_cells...but I don't think you can control that anymore. Now you can choose k
            #filter too small c2
            #print("then is group cluster")
            #print(group.cluster)

            #for(z in levels(as.factor(group.cluster$sample_id))){
            #if(z %in% as.factor(group.cluster$sample_id)){
            counter = 0
            sample_subset <- subset(cluster_subset, subset=sample_id==y)
            #print("then is sample cluster")
            #print(head(sample_subset@meta.data))
            #print(sample_subset)

            #only keep going when cell number > 25
            if(as.integer(length(colnames(sample_subset))) > 25){

            #Dom: choose k accordingly to the n of cells to pool in the parameter
            k = as.integer(length(colnames(sample_subset))/n_cells)+1
            #print("here is sample_subset:")
            #print(colnames(sample_subset))
            #print("Then is k:")
            #print(k)
            if (k<2){
              stop("Error: k < 2 for at least one sample. Change n_cells")
            }
            sample_subset@meta.data$kmeans <- kmeans(x = sample_subset@reductions[["pca"]]@cell.embeddings,centers = k, nstart = nstart)$cluster

            #in previous part, we choose k according to n_cells
            #and run kmeans() with k, which means get k cluster centers
            # we would get k subsubcluster(cell list), each should contain ~ 10 cells
            # call them k_clusters
            # and for each k_cluster(h is the k_cluster name), we generate its cell information


            #print("here is sample_subset:")
            #print(head(sample_subset@meta.data))
            #use k-mean clustering to cluster c3
            for(h in levels(as.factor(sample_subset@meta.data$kmeans))){

              # print("here is h:")
              # print(h)
              # print("here is subsample:")
              # print(sample_subset)
              k.clusters <- subset(sample_subset,subset=kmeans==h)

              #pool cells
              cells <- rownames(k.clusters@meta.data) #get cell rownames for kcluster

              # print("here are cell list!")
              # print(cells)
              # print("cell list end!!")
              # print(length(cells))

              cell_number <- length(cells) #get cell number that would be pooled
              pool <- k.clusters[[assay_name]]@counts[,cells] #get the cell information inside kcluster
              exp_mtx <- as.matrix(pool) #make a matrix of cell information inside kcluster
              sum_total <- rowSums(exp_mtx)

              if (readcounts == "round") {
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
              #cluster_counter= cluster_counter+1
              meta_data <- append(meta_data, paste(y,"_",counter)) # New cell name
              sample_id <- append(sample_id, paste(y))
              group_id <- append(group_id, paste(x))
              cluster_id <- append(cluster_id,paste(z))

              #print("cells:")
              #print(cells)

              ktable.cells <- data.frame(row.names = cells, pooled_cells=rep(paste(y,h,sep = "_"), length(cells)))
              #print("here is ktable.cells")
              #print(ktable.cells)
              ktable <- rbind(ktable,ktable.cells)
              #print("DONE a kmeans")
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

  # print("here is dataset:")
  # print(dataset)
  # print("here is ketable：")
  # print(ktable)

  dataset <- AddMetaData(dataset, metadata = ktable, col.name = "Pooled_kmeans_cells")
  Idents(dataset) <- "Pooled_kmeans_cells"
  print(DimPlot(dataset, split.by = "sample_id",ncol = 4) + NoLegend())

  table(pseudo_cell_seurat@meta.data$sample_id)

  #return(ktable)
  #return(dataset)
  return(pseudo_cell_seurat)
}


#random pooling + mean/sum option + 25 cells filter
random.cellPooling.dev.yiyi <- function(dataset, n_cells= 10, assay_name="RNA", readcounts = "mean"){
  pseudo_cell_mtx <- matrix(, nrow=length(dataset[[assay_name]]@counts@Dimnames[[1]]), ncol=0)

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
      #print(z)
      #cluster_counter = 0
      cluster_subset <- subset(group_subset,subset=cluster_id==z)

      #if the group has the cluster (?! Not sure why this step, to handle multiple clusters?)
      if(z %in% cluster_subset$cluster_id ){
        #print("IF Z...")

        #For each sample/patient/replicate
        for(y in levels(as.factor(dataset$sample_id))){
          # a <- strsplit(y,".",-1)
          # a <- a[[1]][2]
          # if(a==x){
          #print(y)
          counter = 0

          #If the patient belong to the same group z
          if(y %in% as.factor(cluster_subset$sample_id)){
            #print("IF Y...")
            #pool cells
            sample_subset <- subset(cluster_subset, subset= sample_id ==y)  # subsets according to the sample/replicate

            if(as.integer(length(colnames(sample_subset))) > 25){

            real.cells <- rownames(sample_subset@meta.data) #get cell rowname to pool
            cluster_counter = 0

            while (length(real.cells) >=n_cells){  #when there are more than n cells in the cluster
              pool<- sample(real.cells, n_cells, replace = FALSE) #randomly pool n cells from the subset

              cell_id_to_pool <- pool
              # print("Hi Hi here are cell list to pool!")
              # print(cell_id_to_pool)
              # print("cell list end!!")

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

              # print("here is cluster counter:")
              # print(cluster_counter)

              meta_data <- append(meta_data, paste(y,"_",counter)) # New cell name
              sample_id <- append(sample_id, paste(y))
              group_id <- append(group_id, paste(x))
              cluster_id <- append(cluster_id,paste(z))

              rtable.cells <- data.frame(row.names = cell_id_to_pool, pooled_cells=rep(paste(y,cluster_counter,sep = "_"), length(cell_id_to_pool)))
              # print(head(rtable.cells))
              rtable <- rbind(rtable,rtable.cells)
              #print("DONE a kmeans")
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
  print(DimPlot(dataset, split.by = "sample_id",ncol = 4) + NoLegend())

  table(pseudo_cell_seurat@meta.data$sample_id)

  return(pseudo_cell_seurat)
}

# split and give it new replicate for pesudobulk problem
replicate.splitting <- function(dataset, sample_to_split = "ZC3H20_KO", subid1 = "ZC3H20_KO_1",
                                 subid2 = "ZC3H20_KO_2"){
  #Subset.Seurat <- Seurat.tbrucil
  dataset$sample_id <- as.character(dataset$sample_id)

  #when sample_id == sample_to_split, replace sample_id to subid1 and subid2
  seurat_to_split = subset(x = dataset, subset = sample_id == sample_to_split)
  len = length(seurat_to_split$sample_id)/2
  cells <- rownames(seurat_to_split@meta.data) #get all cell rowname to split
  #randomly get len cells
  split1 <- sample(cells, len, replace = FALSE) # split1 have cell name that need to replace by subid1
  split2 <- c() # split2 have cell name that need to replace by subid2
  for(i in 1:length(cells)){
    if (! (cells[i] %in% split1)){
      split2 <- c(split2, cells[i])
    }
  }

  #replace their sample_id

  for (i in 1:length(rownames(dataset@meta.data))) {
    print("cell ID:")
    print(rownames(dataset@meta.data)[i])
    if(rownames(dataset@meta.data)[i] %in% split1) {
      print("for split1")
      print(dataset$sample_id[[i]])
      dataset$sample_id[[i]] <- subid1
      print("after replace")
      print(dataset$sample_id[[i]])
    } else if(rownames(dataset@meta.data)[i] %in% split2){
      dataset$sample_id[[i]] = subid2
    }
  }

  dataset$sample_id <- as.factor(dataset$sample_id)
  return(dataset)

}

#random pooling + mean/sum/10X option + 25 cells filter
random.cellPooling.dev.yiyi.2 <- function(dataset, n_cells= 10, assay_name="RNA", readcounts = "mean"){
  pseudo_cell_mtx <- matrix(, nrow=length(dataset[[assay_name]]@counts@Dimnames[[1]]), ncol=0)

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
      #print(z)
      #cluster_counter = 0
      cluster_subset <- subset(group_subset,subset=cluster_id==z)

      #if the group has the cluster (?! Not sure why this step, to handle multiple clusters?)
      if(z %in% cluster_subset$cluster_id ){
        #print("IF Z...")

        #For each sample/patient/replicate
        for(y in levels(as.factor(dataset$sample_id))){
          # a <- strsplit(y,".",-1)
          # a <- a[[1]][2]
          # if(a==x){
          #print(y)
          counter = 0

          #If the patient belong to the same group z
          if(y %in% as.factor(cluster_subset$sample_id)){
            #print("IF Y...")
            #pool cells
            sample_subset <- subset(cluster_subset, subset= sample_id ==y)  # subsets according to the sample/replicate

            if(as.integer(length(colnames(sample_subset))) > 25){

              real.cells <- rownames(sample_subset@meta.data) #get cell rowname to pool
              cluster_counter = 0

              while (length(real.cells) >=n_cells){  #when there are more than n cells in the cluster
                pool<- sample(real.cells, n_cells, replace = FALSE) #randomly pool n cells from the subset

                cell_id_to_pool <- pool
                # print("Hi Hi here are cell list to pool!")
                # print(cell_id_to_pool)
                # print("cell list end!!")

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
                } else if ( readcounts == "10X"){
                  mean_total <- round(10*(sum_total/n_cells))
                  mean_total <- data.frame(mean_total) #make a dataframe
                  pseudo_cell_mtx <- cbind(pseudo_cell_mtx, mean_total$mean_total)
                }else {
                  stop("Error: readcounts parameter not known")
                }

                #increase counters
                counter = counter + 1
                cluster_counter= cluster_counter+1

                # print("here is cluster counter:")
                # print(cluster_counter)

                meta_data <- append(meta_data, paste(y,"_",counter)) # New cell name
                sample_id <- append(sample_id, paste(y))
                group_id <- append(group_id, paste(x))
                cluster_id <- append(cluster_id,paste(z))

                rtable.cells <- data.frame(row.names = cell_id_to_pool, pooled_cells=rep(paste(y,cluster_counter,sep = "_"), length(cell_id_to_pool)))
                # print(head(rtable.cells))
                rtable <- rbind(rtable,rtable.cells)
                #print("DONE a kmeans")
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
  print(DimPlot(dataset, split.by = "sample_id",ncol = 4) + NoLegend())

  table(pseudo_cell_seurat@meta.data$sample_id)

  return(pseudo_cell_seurat)
}


##------ old pooling functions ---------
#Can be deleted??

#The original random pooling function
cellPooling2 <- function(dataset, n_cells= 10, assay_name="RNA"){
  pseudo_cell_mtx <- matrix(, nrow=length(dataset[[assay_name]]@counts@Dimnames[[1]]), ncol=0)
  ##pseudo_cell_mtx--------total gene num,4000
  meta_data = c()
  cluster_id = c()

  for (x in levels(as.factor(dataset@active.ident))){ # This is important
    ##group, A and B,by ident()
    print(x)
    counter = 0
    group_subset <- subset(dataset, subset= group_id ==x)
    for (z in levels(as.factor(dataset$cluster_id))) {
      cluster_counter = 0
      cluster_subset <- subset(group_subset,subset=cluster_id==z)
      if(z %in% cluster_subset$cluster_id ){
        for(y in levels(as.factor(dataset$sample_id))){
          # a <- strsplit(y,".",-1)
          # a <- a[[1]][2]
          # if(a==x){
          if(y %in% as.factor(group_subset$sample_id)){
            real.cells_condition <- subset(cluster_subset, subset= sample_id ==y)  #-----get subsets according to their replicates
            real.cells <- rownames(real.cells_condition@meta.data) #-----get their rowname
            while (length(real.cells) >=n_cells){  #---when there are more than 10 cells in the cluster
              pool<- sample(real.cells, n_cells, replace = FALSE) #---randomly get 10 cells from the subset
              real.cells <- subset(real.cells, !(real.cells %in% pool)) #---delete 10 cells in the subset???
              pool <- cluster_subset[[assay_name]]@counts[,pool] #---get the 10 cells(before only name, this step get other information)
              exp_mtx <- as.matrix(pool) #---the 10 cells became a matrix
              sum_total <- rowSums(exp_mtx) #10?
              mean_total <- sum_total/n_cells #10/10=1?
              mean_total <- data.frame(mean_total) #emmmmmm……
              pseudo_cell_mtx <- cbind(pseudo_cell_mtx, mean_total$mean_total)
              counter = counter + 1
              cluster_counter= cluster_counter+1
              meta_data <- append(meta_data, paste(y,"_",counter)) # This is important
              cluster_id <- append(cluster_id,paste(z,"_",cluster_counter))
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
  pseudo_cell_seurat$cluster_id <- cluster_id
  return(pseudo_cell_seurat)
}

#First version of k means working
#It works, but create the same number of virtual cells for each cluster (in this case 25), while I want something different for each cluster
cellPooling.kmean.dev.dom <- function(dataset, k=25, nstart=100, assay_name="RNA"){
  pseudo_cell_mtx <- matrix(, nrow=length(dataset[[assay_name]]@counts@Dimnames[[1]]), ncol=0)

  meta_data = c()
  group_id = c()
  cluster_id = c()
  sample_id = c()
  ktable = data.frame(row.names = rownames(dataset))

  for (x in levels(as.factor(dataset@active.ident))){ # This is important
    ##group, A and B, by ident()
    print(x)
    group_subset <- subset(dataset, subset= group_id ==x)

    #for each group...
    for (z in levels(as.factor(dataset$cluster_id))) {
      print(z)
      #cluster_counter = 0
      cluster_subset <- subset(group_subset,subset=cluster_id==z)

      #if the group has the cluster (?! Not sure why this step, to handle multiple clusters?)
      if(z %in% cluster_subset$cluster_id ){
        print("IF Z...")

        #For each sample/patient/replicate
        for(y in levels(as.factor(dataset$sample_id))){
          # a <- strsplit(y,".",-1)
          # a <- a[[1]][2]
          # if(a==x){
          print(y)
          #counter = 0

          #If the patient belong to the same group z
          if(y %in% as.factor(group_subset$sample_id)){
            print("IF Y...")

            #if(length(group.cluster$sample_id) > 10){ # DOM: should be n_cells...but I don't think you can control that anymore. Now you can choose k
              #filter too small c2
              #print("then is group cluster")
              #print(group.cluster)

            #for(z in levels(as.factor(group.cluster$sample_id))){
            #if(z %in% as.factor(group.cluster$sample_id)){
            counter = 0
            sample_subset <- subset(cluster_subset, subset=sample_id==y)
            #print("then is sample cluster")
            #print(head(sample_subset@meta.data))
            #print(sample_subset)

            sample_subset@meta.data$kmeans <- kmeans(x = sample_subset@reductions[["pca"]]@cell.embeddings,centers = k, nstart = nstart)$cluster
            #print(head(sample_subset@meta.data))
            #use k-mean clustering to cluster c3
            for(h in levels(as.factor(sample_subset@meta.data$kmeans))){
              print(h)
              #print(sample_subset)
              k.clusters <- subset(sample_subset,subset=kmeans==h)

            #pool cells
            #real.cells_condition <- subset(cluster_subset, subset= sample_id ==y)  # subsets according to the sample/replicate
              cells <- rownames(k.clusters@meta.data) #get cell rownames
            #while (length(real.cells) >=n_cells){  #when there are more than n cells in the cluster
              #pool<- sample(k.clusters, n_cells, replace = FALSE) #randomly pool n cells from the subset
              #real.cells <- subset(real.cells, !(real.cells %in% pool)) #delete those n cells from the subset
              pool <- k.clusters[[assay_name]]@counts[,cells] #get the n cells readcounts(before was only names)
              exp_mtx <- as.matrix(pool) #make a matrix and the mean
              sum_total <- rowSums(exp_mtx)
              mean_total <- sum_total/n_cells
              mean_total <- data.frame(mean_total) #make a dataframe
              pseudo_cell_mtx <- cbind(pseudo_cell_mtx, mean_total$mean_total)

              #increase counters
              counter = counter + 1
              #cluster_counter= cluster_counter+1
              meta_data <- append(meta_data, paste(y,"_",counter)) # New cell name
              sample_id <- append(sample_id, paste(y))
              group_id <- append(group_id, paste(x))
              cluster_id <- append(cluster_id,paste(z))

              ktable.cells <- data.frame(row.names = cells, rep(h, length(cells)))
              print(ktable.cells)
              ktable <- merge(ktable,ktable.cells)
              #print("DONE a kmeans")
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
  print(head(ktable))
  print(tail(ktable))
  return(pseudo_cell_seurat)
}

#You can decide what to do with readcounts of pooled cells: "mean" or "sum" them
cellPooling.kmean.dev.dom2 <- function(dataset, n_cells= 10, nstart=100, assay_name="RNA", readcounts = "mean"){
  pseudo_cell_mtx <- matrix(, nrow=length(dataset[[assay_name]]@counts@Dimnames[[1]]), ncol=0)

  meta_data = c()
  group_id = c()
  cluster_id = c()
  sample_id = c()
  ktable = data.frame(row.names = rownames(dataset))

  print("Pooling...")
  for (x in levels(as.factor(dataset$group_id))){ # This is important
    ##group, A and B, by ident()
    print("here is group_id:")
    print(x)
    group_subset <- subset(dataset, subset= group_id ==x)

    #for each group...
    for (z in levels(as.factor(dataset$cluster_id))) {
      #print("here is z, aka cluster_id(after group_id)")
      #print(z)
      #cluster_counter = 0
      cluster_subset <- subset(group_subset,subset=cluster_id==z)

      #if the group has the cluster (?! Not sure why this step, to handle multiple clusters?)
      if(z %in% cluster_subset$cluster_id ){
        #print("IF Z...")

        #For each sample/patient/replicate
        for(y in levels(as.factor(dataset$sample_id))){
          # a <- strsplit(y,".",-1)
          # a <- a[[1]][2]
          # if(a==x){
          print("here is y, aka sample_id:")
          print(y)
          #counter = 0

          #If the patient belong to the same group z
          #if(y %in% as.factor(group_subset$sample_id)){
          if(y %in% as.factor(cluster_subset$sample_id)){
            #print("IF Y...")

            #if(length(group.cluster$sample_id) > 10){ # DOM: should be n_cells...but I don't think you can control that anymore. Now you can choose k
            #filter too small c2
            #print("then is group cluster")
            #print(group.cluster)

            #for(z in levels(as.factor(group.cluster$sample_id))){
            #if(z %in% as.factor(group.cluster$sample_id)){
            counter = 0
            sample_subset <- subset(cluster_subset, subset=sample_id==y)
            print("then is sample cluster")
            print(head(sample_subset@meta.data))
            #print(sample_subset)

            #Dom: choose k accordingly to the n of cells to pool in the parameter
            k = as.integer(length(colnames(sample_subset))/n_cells)+1
            print("here is sample_subset:")
            print(colnames(sample_subset))
            print("Then is k:")
            print(k)
            if (k<2){
              stop("Error: k < 2 for at least one sample. Change n_cells")
            }
            sample_subset@meta.data$kmeans <- kmeans(x = sample_subset@reductions[["pca"]]@cell.embeddings,centers = k, nstart = nstart)$cluster
            #print(head(sample_subset@meta.data))
            #use k-mean clustering to cluster c3
            for(h in levels(as.factor(sample_subset@meta.data$kmeans))){
              #print(h)
              #print(sample_subset)
              k.clusters <- subset(sample_subset,subset=kmeans==h)

              #pool cells
              #real.cells_condition <- subset(cluster_subset, subset= sample_id ==y)  # subsets according to the sample/replicate
              cells <- rownames(k.clusters@meta.data) #get cell rownames
              #while (length(real.cells) >=n_cells){  #when there are more than n cells in the cluster
              #pool<- sample(k.clusters, n_cells, replace = FALSE) #randomly pool n cells from the subset
              #real.cells <- subset(real.cells, !(real.cells %in% pool)) #delete those n cells from the subset
              pool <- k.clusters[[assay_name]]@counts[,cells] #get the n cells readcounts(before was only names)
              exp_mtx <- as.matrix(pool) #make a matrix and the mean
              sum_total <- rowSums(exp_mtx)

              if (readcounts == "round") {
                mean_total <- round(sum_total/n_cells)
                mean_total <- data.frame(mean_total) #make a dataframe
                pseudo_cell_mtx <- cbind(pseudo_cell_mtx, mean_total$mean_total)
              } else if (readcounts == "sum") {
                sum_total <- data.frame(sum_total)
                pseudo_cell_mtx <- cbind(pseudo_cell_mtx, sum_total$sum_total)
              } else if (readcounts == "10X") {
                mean_total <- round(10*(sum_total/n_cells))
                mean_total <- data.frame(mean_total) #make a dataframe
                pseudo_cell_mtx <- cbind(pseudo_cell_mtx, mean_total$mean_total)
              }
              else {
                stop("Error: readcounts parameter not known")
              }

              #increase counters
              counter = counter + 1
              #cluster_counter= cluster_counter+1
              meta_data <- append(meta_data, paste(y,"_",counter)) # New cell name
              sample_id <- append(sample_id, paste(y))
              group_id <- append(group_id, paste(x))
              cluster_id <- append(cluster_id,paste(z))

              ktable.cells <- data.frame(row.names = cells, pooled_cells=rep(paste(y,h,sep = "_"), length(cells)))
              #print(head(ktable.cells))
              ktable <- rbind(ktable,ktable.cells)
              #print("DONE a kmeans")
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

  dataset <- AddMetaData(dataset, metadata = ktable, col.name = "Pooled_kmeans_cells")
  Idents(dataset) <- "Pooled_kmeans_cells"
  print(DimPlot(dataset, split.by = "sample_id") + NoLegend())

  table(pseudo_cell_seurat@meta.data$sample_id)

  #return(ktable)
  #return(dataset)
  return(pseudo_cell_seurat)
}

#First version of random pooling + mean/sum option
random.cellPooling.dev <- function(dataset, n_cells= 10, assay_name="RNA", readcounts = "mean"){
  pseudo_cell_mtx <- matrix(, nrow=length(dataset[[assay_name]]@counts@Dimnames[[1]]), ncol=0)

  meta_data = c()
  group_id = c()
  cluster_id = c()
  sample_id = c()

  print("Pooling...")
  for (x in levels(as.factor(dataset$group_id))){
    ##group, A and B, by ident()
    print(x)
    group_subset <- subset(dataset, subset= group_id ==x)

    #for each group...
    for (z in levels(as.factor(dataset$cluster_id))) {
      #print(z)
      #cluster_counter = 0
      cluster_subset <- subset(group_subset,subset=cluster_id==z)

      #if the group has the cluster (?! Not sure why this step, to handle multiple clusters?)
      if(z %in% cluster_subset$cluster_id ){
        #print("IF Z...")

        #For each sample/patient/replicate
        for(y in levels(as.factor(dataset$sample_id))){
          # a <- strsplit(y,".",-1)
          # a <- a[[1]][2]
          # if(a==x){
          #print(y)
          counter = 0

          #If the patient belong to the same group z
          if(y %in% as.factor(group_subset$sample_id)){
            #print("IF Y...")
            #pool cells
            real.cells_condition <- subset(cluster_subset, subset= sample_id ==y)  # subsets according to the sample/replicate

            real.cells <- rownames(real.cells_condition@meta.data) #get cell rowname
            while (length(real.cells) >=n_cells){  #when there are more than n cells in the cluster
              pool<- sample(real.cells, n_cells, replace = FALSE) #randomly pool n cells from the subset
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
              #cluster_counter= cluster_counter+1
              meta_data <- append(meta_data, paste(y,"_",counter)) # New cell name
              sample_id <- append(sample_id, paste(y))
              group_id <- append(group_id, paste(x))
              cluster_id <- append(cluster_id,paste(z))
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
  return(pseudo_cell_seurat)
}

