# Here is the function to help user prepare the data format CellDEEP required.
# Input shall be an seurat object, output shall be an object with *sample_id*, *group_id* and *cluster_id*.
# User input not only the ids it requires, but also does it need normalization/scale/UMAP generation
# Output object will be a seurat object with sample_id, group_id and cluster_id.



################ This function allow you to easily subset every object to test them ############
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
#' @return Subset Seurat object
#' @export
#'
#' @examples
#'
#'
Subset.seurat.object <- function(Seurat.object, Ident.to.subset, Specific.Idents.to.subset) {
  Moment.Seurat.object = Seurat.object

  print(Moment.Seurat.object)
  for (i in 1:length(Ident.to.subset)) {
    #print(i)
    #print(Ident.to.subset[[i]])
    Seurat::Idents(Moment.Seurat.object) <- Ident.to.subset[[i]]
    #print(Specific.Idents.to.subset[[i]])
    Moment.Seurat.object <- subset(Moment.Seurat.object, idents= Specific.Idents.to.subset[[i]])
    #Moment.Seurat.object <- subset(Moment.Seurat.object, subset = Specific.Idents.to.subset[[i]])
    Moment.Seurat.object <- ScaleData(Moment.Seurat.object)
  }
  #print("#######Done")
  print(Moment.Seurat.object)
  return(Moment.Seurat.object)
}


############### This function check a ground truth (muscat) table vs a list of genes DE ###############
#' @title check a ground truth (muscat) table vs a list of genes DE
#'
#' @rdname check.ground.truth
#'
#' @description This function take a seurat object as input, and subset it as user input.
#'
#' @param ground.truth.table table with ground truth genes
#' @param genes.DE The list of DE genes that user want to check
#' @param title The title to put in output, eg: "MAST"
#' @param verbose Show the info of ground truth or not
#'
#' @return ground truth table VS the DE gene list
#' @export
#'
#' @examples
#'
#'
check.ground.truth <- function(ground.truth.table, genes.DE, title, verbose=FALSE) {
  #ground.truth.table = table with ground truth genes
  #genes.DE = list of genes DE that you want to check

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

##### This function adjust input data ######

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
#' @return
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
    Original.Seurat.test <- CreateSeuratObject(matrix)
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
    #Subset.Seurat <- subset(x = Original.Seurat, subset = seurat_clusters == "2")
    Subset.Seurat@meta.data[[group_id]] <- droplevels(as.factor(Subset.Seurat@meta.data[[group_id]]))
    Subset.Seurat@meta.data[[cluster_id]] <- droplevels(as.factor(Subset.Seurat@meta.data[[cluster_id]]))
    Subset.Seurat@meta.data[[sample_id]] <- droplevels(as.factor(Subset.Seurat@meta.data[[sample_id]]))
  }
  #Subset.Seurat <- subset(x = Original.Seurat, subset = seurat_clusters == "2")

  table(Subset.Seurat@meta.data[[sample_id]])
  table(Subset.Seurat@meta.data[[group_id]])
  table(Subset.Seurat@meta.data[[cluster_id]])

  if (Need.to.Norm_Transf == TRUE) {
    Subset.Seurat <- NormalizeData(Subset.Seurat)
    Subset.Seurat <- FindVariableFeatures(Subset.Seurat, selection.method = "vst", nfeatures = 2000)
    Subset.Seurat <- ScaleData(Subset.Seurat, features = rownames(Subset.Seurat))
  }

  if (Need.to.Scale == TRUE) {
    Subset.Seurat <- ScaleData(Subset.Seurat, features = rownames(Subset.Seurat))
  }

  if (Need.UMAP.generation == TRUE) {
    Subset.Seurat <- RunPCA(object = Subset.Seurat, npcs = 30)
    ElbowPlot(object = Subset.Seurat, ndims  = 30)
    DimPlot(Subset.Seurat, reduction = "pca")
    Subset.Seurat <- FindNeighbors(Subset.Seurat, dims = 1:20)
    Subset.Seurat <- FindClusters(Subset.Seurat,resolution = 0.1)
    Subset.Seurat <- RunUMAP(Subset.Seurat, dims = 1:20)
  }

  #Subset.Seurat <- subset(x = emma, subset = seurat_clusters == "2")

  #SaveRDS
  if (is.null(Ident.to.subset) == FALSE) {
    saveRDS(Subset.Seurat, paste(folder_middle_path, Subset.filename.tosave, sep = "/"))
  }

  rm(Original.Seurat) #You don't need it and we save memory

  Subset.Seurat@meta.data[["group_id"]] <- Subset.Seurat@meta.data[[group_id]]
  Subset.Seurat@meta.data[["sample_id"]] <- Subset.Seurat@meta.data[[sample_id]]
  Subset.Seurat@meta.data[["cluster_id"]] <- Subset.Seurat@meta.data[[cluster_id]]

  # drop extra levels in Subset.Seurat object
  Subset.Seurat@meta.data[["group_id"]] <- droplevels(as.factor(Subset.Seurat@meta.data[["group_id"]]))
  Subset.Seurat@meta.data[["cluster_id"]] <- droplevels(as.factor(Subset.Seurat@meta.data[["cluster_id"]]))
  Subset.Seurat@meta.data[["sample_id"]] <- droplevels(as.factor(Subset.Seurat@meta.data[["sample_id"]]))

  # #Check if folders exist otherwise create them
  # if (dir.exists(folder_Results_path)==FALSE) {
  #   dir.create(folder_Results_path)
  # }
  #
  # folder_middle_path = paste(folder_Results_path,folder_name,sep = "/")
  # if (dir.exists(folder_middle_path)==FALSE) {
  #   dir.create(folder_middle_path)
  # }
  #
  # file_path = paste(folder_middle_path,output_filename,sep = "/")
  #
  # folder_Results_plot=paste(folder_middle_path,"Plots",sep = "/")
  # if (dir.exists(folder_Results_plot)==FALSE) {
  #   dir.create(folder_Results_plot)
  # }

  # if (DE.group == group_id){
  #   DE.group = "group_id"
  # } else if (DE.group == sample_id){
  #   DE.group = "sample_id"
  # } else if (DE.group == cluster_id){
  #   DE.group = "cluster_id"
  # }

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
    # else {
    #   write("No Ground Truth.\n",file_path)
    #   rm(Ground.truth)
    # }
  }
  # else {
  #   write("No Ground Truth.\n",file_path)
  # }

  return(Subset.Seurat)

}


