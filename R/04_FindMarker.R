


############### This function return genes DE ###############

#' @title The function to return DE genes
#'
#' @param dataset A seurat object
#' @param test.use method to run DE analysis.eg: "MAST", "DESeq2"....(check in seurat FindMarker command)
#' @param ident.1 The 1st ident to compare, can be missing
#' @param ident.2 The 2nd ident to comapre, can be missing
#' @param group.by The group to run DE analysis
#' @param assay The assay to run analysis, "RNA" as default
#' @param p_cutoff P-value cutoff for DE genes, 0.05 as default
#'
#' @return
#' @export
#'
#' @examples
return.DE <- function(dataset, test.use = "wilcox", DE.ident.1, DE.ident.2, DE.group, assay="RNA", p_cutoff = 0.05, name.only = TRUE) {
  if(name.only == TRUE){
    markers <- Seurat::FindMarkers(dataset, ident.1 = DE.ident.1, ident.2 = DE.ident.2, group.by = DE.group, test.use = test.use)
    #DE <- rownames(markers[markers$p_val_adj<0.05,])
    DE <- rownames(subset(markers, p_val_adj < p_cutoff, na.rm = T))
    print(length(DE))
    print(head(DE))
    return(DE)
  }else if(name.only == FALSE){
    markers <- Seurat::FindMarkers(dataset, ident.1 = DE.ident.1, ident.2 = DE.ident.2, group.by = DE.group, test.use = test.use)
    markers <- subset(markers, p_val_adj < p_cutoff, na.rm = T)
    return(markers)
  }

}

##########This function running FindMarker with Pooling#####
#' @title FindMarker with CellDEEP
#'
#' @param object A Seurat Object
#' @param ident.1 1st ident to compare
#' @param ident.2 2nd ident to compare
#' @param group.by  the group to compare
#' @param slot
#' @param reduction
#' @param features
#' @param logfc.threshold
#' @param pseudocount.use
#' @param test.use
#' @param Pool
#' @param readcounts
#' @param n_cells
#' @param assay
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
FindMarker.CellDEEP <- function(object,
                       ident.1 = NULL,
                       ident.2 = NULL,
                       group.by = NULL,
                       test.use = "wilcox",
                       Pool = FALSE,
                       readcounts = "sum",
                       n_cells = 10,
                       assay = "RNA",
                       cell_cutoff = 25,
                       pool_way = "kmean",
                       name.only = TRUE,
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
                            name.only = name.only,
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
                              name.only = name.only,
                              ...)
    }else if(pool_way == "random"){
      #Pooling the cell
      pooled.object <- random.cellPooling.dev.yiyi(object, readcounts = readcounts, n_cells= n_cells, assay_name = assay)

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
                              name.only = name.only,
                              ...)
    }else{
      print("Wrong pooling way, kmean or random")
    }

    }

  #return(pooled.object)
  return(de.markers)
}
