#' multiSampleComparisonClonalCN Compare the clonal Copy Number of multiple samples.
#'
#' @param listCountMtx Named list of raw count matrix of samples
#' @param analysisName Name of the analysis (default "all")
#' @param organism Organism to be analysed (optional - "mouse" or "human" - default "human")
#' @param par_cores number of cores (default 20)
#'
#' @return
#' @export
#'
#' @examples 
#' 
multiSampleComparisonClonalCN_AGreorder <- function(listCountMtx, analysisName = "all", organism = "human" , par_cores = 20){
  
  resList <- lapply(names(listCountMtx), function(x) {
    pipelineCNA(listCountMtx[[x]], sample = x, SUBCLONES = FALSE, ClonalCN = TRUE, par_cores = par_cores, organism=organism)
  })
  names(resList) <- names(listCountMtx)
  
  print("resList generated")
  
  sampleAlterList <- lapply(names(listCountMtx), function(x) {
    SCEVAN:::analyzeSegm(x, nSub = 0)
  })
  names(sampleAlterList) <- paste0(names(listCountMtx),"_subclone", 1:length(names(listCountMtx)))
  
  names(sampleAlterList) <- paste0(analysisName,"_subclone", 1:length(names(listCountMtx)))
  
  diffList <- SCEVAN:::diffSubclones(sampleAlterList, analysisName, nSub = length(names(listCountMtx)))
  
  diffList <- SCEVAN:::testSpecificAlteration(diffList, length(names(listCountMtx)), analysisName)
  
  print("diffList generated")
  
  genesMtx <- lapply(listCountMtx, function(x) rownames(x))
  
  genesMtx <- sort(unique(unlist(genesMtx)))
  
  genesMtx <- data.frame(x = rep(0, length(genesMtx)), row.names = genesMtx)
  
  annot_mtx <- SCEVAN::annotateGenes(genesMtx)
  
  print("genes annotated")
  
  for(i in 1:length(names(listCountMtx))){
    names(diffList) <- gsub(paste0("subclone",i), names(listCountMtx)[i], names(diffList))
  }
  
  names(diffList) <- gsub("clone", "shared", names(diffList))
  
  outputAnalysis <- list(resList, diffList)
  
  save(outputAnalysis, file = paste0("./output/",analysisName,"_outputAnalysis.RData"))
  
  print("outputAnalysis saved")
  
  oncoHeat <- SCEVAN:::annoteBandOncoHeat(annot_mtx, diffList, length(names(listCountMtx)), organism)
  
  rownames(oncoHeat) <- names(listCountMtx)
  
  print("oncoHeat complete")
  
  SCEVAN:::plotAllClonalCN(names(listCountMtx), name = analysisName)
  
  print("plotAllClonalCN complete")
  
  outputAnalysis
  
  #SCEVAN:::plotOncoHeatSubclones(oncoHeat, length(names(listCountMtx)), analysisName, NULL, organism)
  
  #print("plotOncoHeatSubclones complete")
  
}
