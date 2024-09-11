################################################################################
####### -------                    Goals                         ------- #######

# Processing the primary samples only with inferCNV

# Running inferCNV in this manner:
#    As a whole dataset where normals are pooled across patient and each 
#    tumor is separated by patient (ex. 001_EWS, 005_EWS, etc.)
#    This will use analysis_mode='samples'


################################################################################
####### -------                 Load R packages                  ------- #######

library(Seurat)
library(ggplot2)
library(infercnv)
library(httr)
library(plyr)
options(scipen = 100)


################################################################################
####### -------                  Define params                   ------- #######

setwd("/Volumes/Partition 2/Core/Hayashi_Ewing_scRNAseq_Dec2022/")

so <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_primary_celltype.RDS")
Idents(so)

output <- "/Volumes/Partition 2/Core/Hayashi_Ewing_scRNAseq_Dec2022/Analysis/Annotation/Primary/Normalized/inferCNV/"
dir.create(output, recursive = T)

tumor_cells <- c("Ewing sarcoma", "Ewing sarcoma proliferating") # Used only for clone step
cell_type_col <- "Cell.Type"

################################################################################
####### -------             Make metadata columns                ------- #######

so$normal_tumor <- so@meta.data[,cell_type_col] %in% tumor_cells
so$normal_tumor <- gsub("TRUE", "Tumor", so$normal_tumor)
so$normal_tumor <- gsub("FALSE", "Normal", so$normal_tumor)
so$normal_tumor_patient <- paste0(so$Sample, "_", so$normal_tumor)
so$normal_tumor_patient[grep("Normal", so$normal_tumor_patient)] <- "Normal"
table(so$normal_tumor_patient)

# In both cases, "Normal" denotes the normal cells - but each run has a different column
normal_cells <- "Normal"


################################################################################
####### -------              inferCNV                            ------- #######

# Whole dataset
# Grouped by samples' tumor cells
# Normals defined

column_identity <- "normal_tumor_patient"

output_run <- paste0(output, "Primary_whole_dataset_grouped_by_sample_w_normal_defined/")
dir.create(output_run,
           recursive = T)
output_run_select <- paste0(output_run, "/select_results/")
dir.create(output_run_select,
           recursive = T)

#################   Run inferCNV
  
  #################    Create a cell annotations file
  cell_ann <- cbind.data.frame(Cell = rownames(so@meta.data),
                               Cluster = so@meta.data[,column_identity])
  # Need to remove and clusters with <2 cells since inferCNV won't do it automatically
  cell_ann_table <- data.frame(table(cell_ann$Cluster))
  clusters_keep <- as.character(cell_ann_table$Var1[cell_ann_table$Freq>1])
  cell_ann <- cell_ann[cell_ann$Cluster %in% clusters_keep ,]
  cell_ann_out <- paste0(output, "cell_ann/")
  dir.create(cell_ann_out,
             recursive = T)
  write.table(cell_ann,
              file = paste0(cell_ann_out, "whole_dataset.txt"),
              row.names = F, col.names = F, sep = "\t", quote = F)
  
  # Need to remove any references that aren't in this sample
  normal_cells_sub <- normal_cells[normal_cells %in% cell_ann$Cluster]
  
  #################   Get raw counts - need to re-subset based on above annotation file
  
  so_sub <- subset(so,
                   cells = colnames(so)[colnames(so) %in% cell_ann$Cell])
  exp.rawdata <- as.matrix(so@assays$RNA@counts)
  
  # Create infercnv object
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix= exp.rawdata,
                                      annotations_file=paste0(cell_ann_out, "whole_dataset.txt"),
                                      delim="\t",
                                      gene_order_file="helper_data/infercnv_gene_position.txt",
                                      ref_group_names= normal_cells_sub)
    
  # perform infercnv operations to reveal cnv signal
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir= output_run,  # dir is auto-created for storing outputs
                               cluster_by_groups=T,   # cluster
                               HMM_type='i6', 
                               denoise=T,
                               HMM=T,
                               analysis_mode = "samples",
                               num_threads = 1,
                               BayesMaxPNormal = 0.25,
                               output_format = "pdf" 
  )


  
  