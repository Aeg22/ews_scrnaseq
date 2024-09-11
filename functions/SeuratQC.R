
SeuratQC <- function(so = NULL,
                     category = NULL,
                     outdir = NULL,
                     violin_pt_size = 0,
                     QC_umi_violin = T,
                     QC_genes_violin = T,
                     QC_MT_violin = T,
                     QC_Sscore_violin = T){
  
  
  ####################################################################
  ##### ----- Arguments
  ####################################################################
  
  # so                    seurat object
  # category              name of column in meta data to perform QC across
  # outdir                Output directory
  # violin_pt_size        Size of data points for violin plots. 0 will show no points.
  # QC_...                Should this QC be run?
  
  
  ####################################################################
  ##### ----- Input warnings
  ####################################################################
  
  # Null warnings
  if(is.null(so) == TRUE){stop("so must be specified")}
  if(is.null(category) == TRUE){stop("category must be specified")}
  if(is.null(outdir) == TRUE){stop("outdir must be specified")}
  
  # Check for outdir existence 
  if(file.exists(outdir) == FALSE){stop("outdir does not exist")}
  
  
  ####################################################################
  ##### ----- Violin plots
  ####################################################################
  
  # Violin plot vectors
  violin_QC <- c(QC_umi_violin,
                 QC_genes_violin,
                 QC_MT_violin,
                 QC_Sscore_violin)
  features <- c("nCount_RNA",
                "nFeature_RNA",
                "percent_mito",
                "S.Score")
  # Subset features to only those requested
  features <- features[violin_QC]
  
  # Make violin plots
  for(feature in features){
    violin_plot <- VlnPlot(
      so,
      features = feature,
      ncol     = 1,
      pt.size  = violin_pt_size,
      group.by = category
    ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    ggsave(paste0(outdir, "/", feature, ".pdf"),
           violin_plot)
  }
  
}