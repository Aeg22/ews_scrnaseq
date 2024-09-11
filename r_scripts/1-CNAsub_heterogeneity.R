################################################################################
####### -------                    Goals                         ------- #######

# CNA subclusters were identified in some primary samples 

# Here, I will CNA heterogeneity by comparing these subclusters against
# several metrics

# This script will loop through each sample independently 

################################################################################
####### -------                 Load R packages                  ------- #######

library(Seurat)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(msigdbr)
library(dplyr)
library(plyr)
library(clusterProfiler)
library(openxlsx)
library(pheatmap)
library(ggpubr)

source("functions/SeuratQC.R")


################################################################################
####### -------                  Define params                   ------- #######

setwd("/Volumes/Partition 2/Core/Hayashi_Ewing_scRNAseq_Dec2022/")

output <- "Analysis/Annotation/Primary/Normalized/CNAsub_heterogeneity/"

so_all <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_primary_celltype_CNAsubs.RDS")

sample_col <- "Sample"
subcluster_col <- "CNAsubclusters"
nontumor <- c("none") 

samples <- c("005_Primary", "061_Primary")

subcluster_FAM_FDR_threshold <- 0.25

################################################################################
####### -------     Loop through each sample with subclusters    ------- #######

for(sample in samples){
  
  print(sample)
  output_sample <- paste0(output, "/", sample, "/")
  dir.create(output_sample, recursive = T)
  
  ##############################################################################
  ####### -------         Subset so to sample-level              ------- #######
  
  # Subset to sample-level and to only tumor
  so_all_meta <- so_all@meta.data
  so <- so_all[,so_all_meta[,sample_col] %in% sample & 
                 !(so_all_meta[,subcluster_col] %in% nontumor)]
  
  subcluster_numbers <- data.frame(unclass(table(so@meta.data[,subcluster_col])))
  subcluster_numbers <- cbind.data.frame(Cluster = rownames(subcluster_numbers),
                                         Cells = subcluster_numbers[,1])
  write.table(subcluster_numbers,
              file = paste0(output_sample, "CNA_subcluster_cell_numbers.txt"),
              sep = "\t", row.names = F, quote = F)
  
  ##############################################################################
  ####### -------     scale, dimension reduction, cluster        ------- #######
  
  # Sometimes this won't be necessary if performed before but this is the first
  # time each sample has been subset by itself
  
  # Using same parameters as the whole dataset (pre-integration) except n_dim
  # because there is less heterogeneity now
  
  assay <- "RNA"
  DefaultAssay(so) <- assay
  so <- FindVariableFeatures(so,
                             selection.method = "vst")
  so <- ScaleData(object = so, 
                  vars.to.regress = c("nCount_RNA",
                                      "percent_mito"))
  so <- RunPCA(so,
               features = VariableFeatures(so),
               verbose = F,
               ndims.print = 0)
  ElbowPlot(so,
            ndims = 40)
  n_dim <- 10 # Previously, I used 30 but there is less heterogeneity now
  so <- RunUMAP(so, 
                reduction = "pca", 
                dims = 1:n_dim)
  so <- FindNeighbors(so, 
                      reduction = "pca",
                      dims = 1:n_dim, 
                      k.param = 20) 
  resolutionsCalculated <- c( 0.5)
  so <- FindClusters(so, 
                     resolution = resolutionsCalculated)
  resolutionChosen <- 0.5
  so$CNArnaCluster <- so@meta.data[,paste(assay, "_snn_res.", resolutionChosen, sep = "")] 
  
  
  ##############################################################################
  ####### -------         Plot UMAP by CNA and RNA cluster       ------- #######
  
  umap_cluster <- DimPlot(so, reduction = "umap",
                          label = T,
                          group.by = "CNArnaCluster") # labels by cluster
  umap_cluster
  ggsave(paste0(output_sample, "/UMAP_CNArnaCluster.PDF"), 
         plot = umap_cluster)
  
  umap_CNA <- DimPlot(so, reduction = "umap",
                          label = T,
                          group.by = subcluster_col) # labels by cluster
  umap_CNA
  ggsave(paste0(output_sample, "/UMAP_CNAsubcluster.PDF"), 
         plot = umap_CNA)
  
  ##############################################################################
  ####### -------              QC plots by cluster              ------- #######
  
  output_res_qc <- paste0(output_sample, "/QC/")
  dir.create(output_res_qc)
  
  SeuratQC(so = so,
           category = subcluster_col,
           outdir = output_res_qc)
  
  
  ##############################################################################
  ####### -------             Save RDS by sample                 ------- #######
  
  saveRDS(so,
          paste0(output_sample, sample, "_so.RDS"))
  
}






