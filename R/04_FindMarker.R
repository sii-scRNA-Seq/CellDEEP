


############### This function return genes DE ###############
return.DE <- function(dataset, test.use, ident.1 = DE.ident.1, ident.2 = DE.ident.2, group.by = DE.group, assay="RNA") {
  markers <- FindMarkers(dataset, ident.1 = DE.ident.1, ident.2 = DE.ident.2, group.by = DE.group, test.use = test.use)
  #DE <- rownames(markers[markers$p_val_adj<0.05,])
  DE <- rownames(subset(markers, p_val_adj < 0.05, na.rm = T))
  print(length(DE))
  print(head(DE))
  return(DE)
}

##########This function running FindMarker with Pooling#####
FindMarker <- function(object,
                       ident.1 = NULL,
                       ident.2 = NULL,
                       group.by = NULL,
                       subset.ident = NULL,
                       assay = NULL,
                       slot = "data",
                       reduction = NULL,
                       features = NULL,
                       logfc.threshold = 0.1,
                       pseudocount.use = 1,
                       test.use = "wilcox",
                       min.pct = 0.01,
                       min.diff.pct = -Inf,
                       verbose = TRUE,
                       only.pos = FALSE,
                       max.cells.per.ident = Inf,
                       random.seed = 1,
                       latent.vars = NULL,
                       min.cells.feature = 3,
                       min.cells.group = 3,
                       mean.fxn = NULL,
                       fc.name = NULL,
                       base = 2,
                       densify = FALSE,
                       Pool = FALSE,
                       readcounts = "sum",
                       n_cells = 10,
                       assay = "RNA",
                       ...
                       ){

  #if Pool = FALSE, run FindMarker directly
  if (Pool == FALSE) {
    de.markers <- return.DE(object,
                             ident.1 = ident.1,
                             ident.2 = ident.2,
                             group.by = group.by,
                             subset.ident = subset.ident,
                             assay = assay,
                             slot = slot,
                             reduction = reduction,
                             features = features,
                             logfc.threshold = logfc.threshold,
                             pseudocount.use = pseudocount.use,
                             test.use = test.use,
                             min.pct = min.pct,
                             min.diff.pct = min.diff.pct,
                             verbose = verbose,
                             only.pos = only.pos,
                             max.cells.per.ident = max.cells.per.ident,
                             random.seed = random.seed,
                             latent.vars = latent.vars,
                             min.cells.feature = min.cells.feature,
                             min.cells.group = min.cells.group,
                             mean.fxn = mean.fxn,
                             fc.name = fc.name,
                             base = base,
                             densify = densify,
                             ...)
  }else #When pool = TRUE, run CellDEEP progress
    {

    print("Start Pooling.....")

    #Pooling the cell
    pooled.object <- cellPooling.kmean.dev.yiyi(object, readcounts = readcounts, n_cells= n_cells, assay_name = assay)

    Print("FindMarker running.....")

    #Run FindMarker and return DE gene list
    de.markers <- return.DE(pooled.object,
                            ident.1 = ident.1,
                            ident.2 = ident.2,
                            group.by = group.by,
                            subset.ident = subset.ident,
                            assay = assay,
                            slot = slot,
                            reduction = reduction,
                            features = features,
                            logfc.threshold = logfc.threshold,
                            pseudocount.use = pseudocount.use,
                            test.use = test.use,
                            min.pct = min.pct,
                            min.diff.pct = min.diff.pct,
                            verbose = verbose,
                            only.pos = only.pos,
                            max.cells.per.ident = max.cells.per.ident,
                            random.seed = random.seed,
                            latent.vars = latent.vars,
                            min.cells.feature = min.cells.feature,
                            min.cells.group = min.cells.group,
                            mean.fxn = mean.fxn,
                            fc.name = fc.name,
                            base = base,
                            densify = densify,
                            ...)
    }

  return(de.markers)
}
