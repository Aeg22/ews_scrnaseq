


FindNeighbors_FindClusters_PrintUMAP <- function(
  
  # Parameters
  sobj = NULL, # Seurat object
    # Wtihin sobj, this function expects umap calculated and metadata to have: "Phase"
  variable = NULL,
  neighbors_reduction = "pca",
  neighbors_dim = 1:20, # Use ElbowPlot(sobj, ndims = 40)
  neighbors_k.param = 20,
  resolutionsCalculated = seq(0.1, 1, by = 0.1)

){
  
  print("It may be especially desirable to specify neighbors_dim (default = 1:20) and ElbowPlot(sobj, ndims = 40) (default = seq(0.1, 1, by = 0.1)). ")
  
  if(is.null(sobj)){
    stop("sobj is required. Use a Seurat object.")
  }
  if(is.null(sobj)){
    stop("variable is required. Use a character/factor column present in sobj (ex. 'Sample'.")
  }
  
  library(gridExtra)
  library(Seurat)
  
  
  sobj_function <- FindNeighbors(sobj, 
                                 reduction = neighbors_reduction,
                                 dims = neighbors_dim, # 
                                 k.param = neighbors_k.param) # Constructs k and shared nearest neighbor graphs
  
  sobj_function <- FindClusters(sobj_function, 
                                resolution = resolutionsCalculated) # This is running Louvain
  
  assay <- DefaultAssay(sobj_function)
  
  # look at what clusters would be labelled as at each resolution
  for (i in resolutionsCalculated){
    # set resolution
    Idents(object = sobj_function) <- paste(assay, "_snn_res." , i, sep = "") # name is picked from column name of resolutions
    
    # general umap
    allClusters <- DimPlot(sobj_function,
                           reduction = "umap",
                           label = TRUE,
                           label.size = 6) + ggtitle(paste(assay, "_snn_res." , i, sep = "")) + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
    
    # UMAP of cells in each cluster by variable - split
    UMAP_splitByVariable <- DimPlot(sobj_function, 
                             label = F, 
                             split.by = variable)  + NoLegend() + ggtitle(paste0("Clusters split by ", variable)) + theme(plot.title = element_text(hjust = 0.5))
    
    # UMAP of variable
    UMAPvariable <- DimPlot(sobj_function, 
                             label = F,
                            group.by = variable) + ggtitle(paste0("Clusters by ", variable)) + theme(plot.title = element_text(hjust = 0.5))
    
    # Explore whether clusters segregate by cell cycle phase
    splitByPhase <- DimPlot(sobj_function,
                            label = T, 
                            split.by = "Phase")  + NoLegend() + ggtitle("Clusters split by cell cycle phase") + theme(plot.title = element_text(hjust = 0.5))
    
    # Print plot
    grid.arrange(allClusters, UMAPvariable, splitByPhase, UMAP_splitByVariable + remove("x.text"), 
                 ncol = 2, nrow = 2)
    
  }
    
    return(sobj_function) # Return the new sobj
  
    print("It may be especially desirable to specify neighbors_dim (default = 1:20) and ElbowPlot(sobj, ndims = 40) (default = seq(0.1, 1, by = 0.1)). ")

}

