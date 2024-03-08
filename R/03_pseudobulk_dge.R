# library(muscat)
# library(Seurat)
# library(data.table)
# library(edgeR)
#
# VERSION <- '0.1.0'

# ------Used version-------------------

# Run differential gene expression analysis on each cluster of the
# singleCellExperiment.
#
# sce: SingleCellExperiment object with cluster_id, sample_id, group_id in slot
# colData. group_id
#
# contrasts: Vector of contrasts comparing `group_id` levels. E.g.
# `c('Heathy-Resistant')`
#
# ... Further arguments for differential expression analysis passed to
# muscat::pbDS.  E.g. `method='limma-voom'`

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
    design <- design[match(samples, sample_id)]

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

    contrast_mat <- makeContrasts(contrasts= contrasts, levels = mm)
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

    ##adjust pb's ident & orig.ident
    # pb$ident <- pb$virtual_sample_id
    # pb$orig.ident <- pb$virtual_sample_id

    # Run DGE for each contrast and cluster_id
    print(pb@metadata$n_cells)
    print(pb@metadata)
    dge <- pbDS(pb, design= mm, verbose= verbose, contrast= contrast_mat,min_cells = 1,...)
    #dge <- pbDS(pb, design= mm, verbose= verbose, contrast= contrast_mat,min_cells = 1)
    #print("Let's see what's the problem")
    #print(dge$table)
    #print(names(dge$table)) #A-B
    #print(names(dge$table[[1]])) #Cluster0

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

pre_pesudobulk <- function(seurat_object){
  #NEW function - add virtual cell
  #make a loop
  new_ids = c()
  for (i in 1:length(seurat_object$sample_id)) {
    #print(i)
    #print(paste(seurat_object$sample_id[i][[1]], names(seurat_object$sample_id[i]), sep = "__"))
    #new_ids=c(new_ids, paste(seurat_object$sample_id[i][[1]], names(seurat_object$sample_id[i]), sep = "__"))
    name <- names(seurat_object$sample_id[i])
    new_ids=c(new_ids, gsub(" ","",name))
    #print(new_ids)
  }

  #new_ids
  seurat_object$virtual_sample_id <- new_ids

  ##Transfer the new ids

  #id transfer for running pesudobulk part
  seurat_object$real_sample_id <- seurat_object$sample_id
  seurat_object$sample_id <- seurat_object$virtual_sample_id

  return(seurat_object)
}

##------Previous version-------

pseudobulk_dge <- function(sce, contrasts, verbose= FALSE, ...) {
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
  design <- design[match(samples, sample_id)]

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

  contrast_mat <- makeContrasts(contrasts= contrasts, levels = mm)
  if(verbose) {
    cat('Contrast matrix:\n')
    print(contrast_mat)
    cat('\n')
  }

  # muscat requires `experiment_info` to be not NULL but I'm not sure it is
  # needed. Using `pb@metadata[['experiment_info']] <- NA` should also do
  pb@metadata[['experiment_info']] <- design

  # Run DGE for each contrast and cluster_id
  dge <- pbDS(pb, design= mm, verbose= verbose, contrast= contrast_mat, ...)

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

pseudobulk_dge_test <- function(sce, contrasts, verbose= FALSE, ...) {
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

  # Create a design matrix - make sure order is consistent with pseudobulk
  design <- unique(as.data.table(sce@colData[, c('group_id', 'sample_id')]))

  samples <- colnames(pb@assays@data[[1]])
  design <- design[match(samples, sample_id)]

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

  contrast_mat <- makeContrasts(contrasts= contrasts, levels = mm)
  if(verbose) {
    cat('Contrast matrix:\n')
    print(contrast_mat)
    cat('\n')
  }

  # muscat requires `experiment_info` to be not NULL but I'm not sure it is
  # needed. Using `pb@metadata[['experiment_info']] <- NA` should also do
  pb@metadata[['experiment_info']] <- design

  # Run DGE for each contrast and cluster_id
  dge <- pbDS(pb, design= mm, verbose= verbose, contrast= contrast_mat, ...)

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

pseudobulk_dge_pool_withoutfilter <- function(sce, contrasts, verbose= FALSE, ...) {
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
  design <- design[match(samples, sample_id)]

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

  contrast_mat <- makeContrasts(contrasts= contrasts, levels = mm)
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

  ##adjust pb's ident & orig.ident
  # pb$ident <- pb$virtual_sample_id
  # pb$orig.ident <- pb$virtual_sample_id

  # Run DGE for each contrast and cluster_id
  print(pb@metadata$n_cells)
  print(pb@metadata)
  dge <- pbDS(pb, design= mm, verbose= verbose, contrast= contrast_mat,min_cells = 1,filter = "none", ...)
  #dge <- pbDS(pb, design= mm, verbose= verbose, contrast= contrast_mat,min_cells = 1)
  #print("Let's see what's the problem")
  #print(dge$table)
  #print(names(dge$table)) #A-B
  #print(names(dge$table[[1]])) #Cluster0

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
