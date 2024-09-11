# Making changes to save classDf as RDS
# Also adding print steps so I can see where the pipeline might break

pipelineCNA_AG <- function (count_mtx, sample = "", par_cores = 20, norm_cell = NULL, 
          SUBCLONES = TRUE, beta_vega = 0.5, ClonalCN = TRUE, plotTree = TRUE, 
          AdditionalGeneSets = NULL, SCEVANsignatures = TRUE, organism = "human") 
{
  dir.create(file.path("./output"), showWarnings = FALSE)
  start_time <- Sys.time()
  normalNotKnown <- length(norm_cell) == 0
  res_proc <- SCEVAN:::preprocessingMtx(count_mtx, sample, par_cores = par_cores, 
                               findConfident = normalNotKnown, AdditionalGeneSets = AdditionalGeneSets, 
                               SCEVANsignatures = SCEVANsignatures, organism = organism)
  if (normalNotKnown) 
    norm_cell <- names(res_proc$norm_cell)
  res_class <- SCEVAN::classifyTumorCells(res_proc$count_mtx_norm, 
                                  res_proc$count_mtx_annot, sample, par_cores = par_cores, 
                                  ground_truth = NULL, norm_cell_names = norm_cell, SEGMENTATION_CLASS = TRUE, 
                                  SMOOTH = TRUE, beta_vega = beta_vega)
  print(paste("found", length(res_class$tum_cells), "tumor cells"))
  classDf <- data.frame(class = rep("filtered", length(colnames(count_mtx))), 
                        row.names = colnames(count_mtx))
  classDf[colnames(res_class$CNAmat)[-(1:3)], "class"] <- "normal"
  classDf[res_class$tum_cells, "class"] <- "tumor"
  classDf[res_class$confidentNormal, "confidentNormal"] <- "yes"
  end_time <- Sys.time()
  print(paste("time classify tumor cells: ", end_time - start_time))
  if (ClonalCN) 
    SCEVAN:::getClonalCNProfile(res_class, res_proc, sample, par_cores, 
                       organism = organism)
  mtx_vega <- SCEVAN:::segmTumorMatrix(res_proc, res_class, sample, 
                              par_cores, beta_vega)
  print("Ran ClonalCN")
  if (SUBCLONES) {
    res_subclones <- SCEVAN:::subcloneAnalysisPipeline(count_mtx, 
                                              res_class, res_proc, mtx_vega, sample, par_cores, 
                                              classDf, beta_vega, plotTree, organism)
    print("Ran subcloneAnalysisPipeline")
    FOUND_SUBCLONES <- res_subclones$FOUND_SUBCLONES
    print("FOUND_SUBCLONES")
    classDf <- res_subclones$classDf
    print("classDf created")
  }
  else {
    FOUND_SUBCLONES <- FALSE
  }
  print("ran subclones")
  if (!FOUND_SUBCLONES) 
    SCEVAN:::plotCNclonal(sample, ClonalCN, organism)
  print("ran plotCNclonal")
  count_mtx_annot <- res_proc$count_mtx_annot
  save(count_mtx_annot, file = paste0("./output/", sample, 
                                      "_count_mtx_annot.RData"))
  mtx_vega_files <- list.files(path = "./output/", pattern = "_mtx_vega")
  
  # AG added here
  save(classDf, file = paste0("./output/", sample, 
                                      "_classDf.RData"))
  print("saved files")
  
  sapply(mtx_vega_files, function(x) file.remove(paste0("./output/", 
                                                        x)))
  return(classDf)
}