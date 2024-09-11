################################################################################
####### -------                    Goals                         ------- #######

# Processing the primary samples only with inferCNV

# This script is optimized to run an inferCNV subclusters protocol that
# I developed. 

# More details below but I essentially run inferCNV by sample using random_trees
# and then I prune those branches to a more reasonable number based on CNA similarity

# More detailed steps:
# Step 1 - inferCNV by sample using the subclusters mode and random_trees
# Step 2 - inferCNV in sample mode using the random_tree (RT, RF in script) classification
#          This is because inferCNV does not plot those subclusters in Step 1
# Step 2.5 (not scripted) - Manually simplification of the tree structures from Step 1 (viewing in Step 2)
# Step 3 - final inferCNV run using the simplified RT subclusters
# Step 4 - uphyloplot2 using simplified RT subclusters
# Step 5 (not scripted) - Manual CNV annotation of the tree diagrams

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

# Only using raw counts so I could have done a previous version if it were annotated with clusters
so <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_primary_celltype.RDS")

output <- "/Volumes/Partition 2/Core/Hayashi_Ewing_scRNAseq_Dec2022/Analysis/Annotation/Primary/Normalized/inferCNV/"
dir.create(output, recursive = T)

tumor_cells <- c("Ewing sarcoma", "Ewing sarcoma proliferating") # Used only for clone step
cell_type_col <- "Cell.Type"

path_uphyloplot2 <- "~/Documents/Uphyloplot2/"

################################################################################
####### -------             Make metadata columns                ------- #######

# For inferCNV run #1
so$normal_tumor <- so@meta.data[,cell_type_col] %in% tumor_cells
so$normal_tumor <- gsub("TRUE", "Tumor", so$normal_tumor)
so$normal_tumor <- gsub("FALSE", "Normal", so$normal_tumor)
so$normal_tumor_patient <- paste0(so$Sample, "_", so$normal_tumor)
so$normal_tumor_patient[grep("Normal", so$normal_tumor_patient)] <- "Normal"
table(so$normal_tumor_patient)

# In both cases, "Normal" denotes the normal cells - but each run has a different column
normal_cells <- "Normal"


################################################################################
####### -------  inferCNV steps 1 and 2 - random trees w/ pval          ------- #######

samples <- c("005_Primary", "061_Primary")

# By sample
# Normals defined
column_identity <- "normal_tumor_patient"

output_run <- paste0(output, "Step1_randomTrees_subclusPval0000000001/")
dir.create(output_run, recursive = T)
output_run_select <- paste0(output_run, "/select_results/")
dir.create(output_run_select, recursive = T)
cell_ann_out <- paste0(output, "cell_ann/")
dir.create(cell_ann_out, recursive = T)
output_run_step2 <- paste0(output, "Step2_subclusters_all/")
dir.create(output_run_step2, recursive = T)
output_run_select_step2 <- paste0(output_run_step2, "/select_results/")
dir.create(output_run_select_step2, recursive = T)
output_run_step3 <- paste0(output, "Step3_final_subclusters/")
dir.create(output_run_step3, recursive = T)
output_run_select_step3 <- paste0(output_run_step3, "/select_results/")
dir.create(output_run_select_step3, recursive = T)
output_run_step4 <- paste0(output, "Step4_uphyloplot2/")
dir.create(output_run_step4, recursive = T)

#################   Loop through each sample for steps 1 and 2 
# The next step requires manually curation after this loop

for (sample in samples){
  
  #################   Get raw counts
  
  #sample <- "051_Primary" # If I just want to test, skip the loop and subset a sample like this
  
  print(paste0("###########---------------------------########## ", sample, " of ", length(samples), " total"))
  
  so_sub <- subset(so,
                   cells = colnames(so)[so$Sample %in% sample])
  
  #################    Create a cell annotations file based on cell type
  cell_ann <- cbind.data.frame(Cell = rownames(so_sub@meta.data),
                               Cluster = so_sub@meta.data[,column_identity])
  # Need to remove and clusters with <2 cells since inferCNV won't do it automatically
  cell_ann_table <- data.frame(table(cell_ann$Cluster))
  clusters_keep <- as.character(cell_ann_table$Var1[cell_ann_table$Freq>1])
  cell_ann <- cell_ann[cell_ann$Cluster %in% clusters_keep ,]
  
  write.table(cell_ann,
              file = paste0(cell_ann_out, sample, ".txt"),
              row.names = F, col.names = F, sep = "\t", quote = F)
  
  # Need to remove any references that aren't in this sample
  normal_cells_sub <- normal_cells[normal_cells %in% cell_ann$Cluster]
  
  #################   Get raw counts - need to re-subset based on above annotation file
  
  so_sub <- subset(so_sub,
                   cells = colnames(so_sub)[colnames(so_sub) %in% cell_ann$Cell])
  exp.rawdata <- as.matrix(so_sub@assays$RNA@counts)
  
  # Create infercnv object
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix= exp.rawdata,
                                      annotations_file=paste0(cell_ann_out, sample, ".txt"),
                                      delim="\t",
                                      gene_order_file="helper_data/infercnv_gene_position.txt",
                                      ref_group_names= normal_cells_sub)
  
  # perform infercnv operations to reveal cnv signal
  sample_output <- paste0(output_run, sample)
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir= sample_output,  # dir is auto-created for storing outputs
                               cluster_by_groups=F,   # cluster (Changed for uphyloplot2)
                               plot_steps=T, # Changed from false for uphyloplot2
                               HMM_type='i6', # Changed for uphyloplot2
                               denoise=T,
                               HMM=T,
                               analysis_mode = "subclusters",
                               num_threads = 6,
                               tumor_subcluster_partition_method = "random_trees",
                               tumor_subcluster_pval = 0.0000000001,
                               output_format = "pdf" 
  )
  
  # There are so many large files created by this process
  # I am copying one to another folder for simplicity 
  file.copy(from = paste0(sample_output, "/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings"),
            to = paste0(output_run_select, "cell_groupings_",  sample, ".txt"),
            overwrite = T)
  
  # Could use this to delete the mass of files for each run
  #unlink(sample_output,
  #       recursive = T)
  
  
  #################    Create a cell annotations file based on subcluster RF tree
  
  rf_annotation <- read.table(file = paste0(sample_output, "/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings"),
                              sep = "\t", header = T, stringsAsFactors = F)
  rf_annotation <- cbind.data.frame(Cell = rf_annotation$cell,
                                    Subcluster = rf_annotation$cell_group_name)
  
  sum(!(cell_ann$Cell %in% rf_annotation$Cell)) # None should be missing
  
  # Create a file containing the unique subcluster names
  subclusters_unique <- cbind.data.frame(Subcluster = unique(rf_annotation$Subcluster)[-grep(normal_cells, unique(rf_annotation$Subcluster))],
                                         Simplified = "")
  
  # Merge previous annotation with subcluster
  cell_ann_step2 <- join(cell_ann,
                         rf_annotation)
  # Rename all normal subclusters as normal
  cell_ann_step2$Subcluster[grep(normal_cells, cell_ann_step2$Subcluster)] <- normal_cells
  # Simplify tumor subcluster name
  cell_ann_step2 <- cbind.data.frame(Cell = cell_ann_step2$Cell,
                                     Cluster = cell_ann_step2$Subcluster)
  
  # Some subclusters will only have 1 cell (need to remove prior to running)
  cell_ann_table <- data.frame(table(cell_ann_step2$Cluster))
  clusters_keep <- as.character(cell_ann_table$Var1[cell_ann_table$Freq>1])
  cell_ann_step2 <- cell_ann_step2[cell_ann_step2$Cluster %in% clusters_keep ,]
  
  write.table(cell_ann_step2,
              file = paste0(cell_ann_out, sample, "_step2.txt"),
              row.names = F, col.names = F, sep = "\t", quote = F)
  
  so_sub <- subset(so_sub,
                   cells = colnames(so_sub)[colnames(so_sub) %in% cell_ann_step2$Cell])
  exp.rawdata <- as.matrix(so_sub@assays$RNA@counts)
  
  # Create infercnv object - Step 2
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = exp.rawdata,
                                      annotations_file = paste0(cell_ann_out, sample, "_step2.txt"),
                                      delim = "\t",
                                      gene_order_file = "/Volumes/Partition 2/Core/Hayashi_Ewing_scRNAseq_Aug2021/Analysis/Annotation/infercnv/gene_position/your_gen_pos.txt",
                                      ref_group_names = normal_cells_sub)
  
  # perform infercnv operations to reveal cnv signal - step 2
  sample_output_step2 <- paste0(output_run_step2, sample)
  dir.create(sample_output_step2)
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir= sample_output_step2,  # dir is auto-created for storing outputs
                               cluster_by_groups=T,   # cluster (Changed for uphyloplot2)
                               #plot_steps=T, # Changed from false for uphyloplot2
                               HMM_type='i6', # Changed for uphyloplot2
                               denoise=T,
                               HMM=T,
                               analysis_mode = "samples",
                               num_threads = 6,
                               output_format = "pdf" # default: "png"
  )
  file.copy(from = paste0(sample_output_step2, "/infercnv.pdf"),
            to = paste0(output_run_select_step2,  sample, "_infercnv.pdf"),
            overwrite = T)
  file.copy(from = paste0(sample_output_step2, "/infercnv.20_HMM_predHMMi6.hmm_mode-samples.Pnorm_0.5.repr_intensities.pdf"),
            to = paste0(output_run_select_step2,  sample, "_infercnv_bayes.pdf"),
            overwrite = T)
  write.table( subclusters_unique,
               file = paste0(output_run_select_step2,  sample, "_unique_subclusters.txt"),
               sep = "\t", row.names = F, quote = F)
  
  # Deleting excess files
  unlink(sample_output,
         recursive = T)
  unlink(sample_output_step2,
         recursive = T)
  
} 

# For each sample, open step 2's infercnv.pdf and infercnv.17_HMM_predHMMi6.hmm_mode-samples.pdf
# Inspect the subclusters and their RF hierarchy and simplify branches in unique_subclusters.txt
# Once the 'Simplified' column in unique_subclusters.txt has been edited, 
# Step 3, pruned subclusters run as 'samples' for final inferCNV plots
# Step 4, replace the high resolution RF subcluster names with unique_subclusters.txt and run uphyloplot2

################################################################################
####### -------  inferCNV step 3 (and files for step 4)                ------- #######


for (sample in samples){
  
  #################   Read in cell annotation and the replacement subcluster names
  
  replacement_subclusters <- read.table( file = paste0(output_run_select_step2,  sample, "_unique_subclusters.txt"),
                                         sep = "\t", header = T, stringsAsFactors = F)
  
  cell_ann <- read.table(file = paste0(output_run_select, "cell_groupings_",  sample, ".txt"),
                         sep = "\t", header = T, stringsAsFactors = F)
  
  # Replace to simplified clusters
  for(i in 1:nrow(replacement_subclusters)){
    cell_ann$cell_group_name <- gsub(replacement_subclusters$Subcluster[i],
                                     replacement_subclusters$Simplified[i],
                                     cell_ann$cell_group_name)
  }
  
  # Shorten tumor subcluster name
  cell_ann$cell_group_name <- gsub("all_observations.all_observations.", "subcluster.", cell_ann$cell_group_name)
  cell_ann$cell_group_name[grep(normal_cells, cell_ann$cell_group_name)] <- normal_cells
  
  # Write file for uphyloplot2
  uphylo <- cell_ann[-grep(normal_cells, cell_ann$cell_group_name),] # Only want tumor cells
  write.table(uphylo,
              paste0(output_run_step4, sample, "_uphyloplot2.txt"),
              sep = "\t", row.names = F, quote = F)
  
  # Write cell_ann for step 3 run
  cell_ann <- cell_ann[,c(2,1)]
  write.table(cell_ann,
              file = paste0(cell_ann_out, sample, "_step3.txt"),
              row.names = F, col.names = F, sep = "\t", quote = F)
  
  so_sub <- subset(so,
                   cells = colnames(so)[colnames(so) %in% cell_ann$cell])
  exp.rawdata <- as.matrix(so_sub@assays$RNA@counts)
  
  # Create infercnv object - Step 3
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = exp.rawdata,
                                      annotations_file = paste0(cell_ann_out, sample, "_step3.txt"),
                                      delim = "\t",
                                      gene_order_file = "/Volumes/Partition 2/Core/Hayashi_Ewing_scRNAseq_Aug2021/Analysis/Annotation/infercnv/gene_position/your_gen_pos.txt",
                                      ref_group_names = normal_cells)
  
  # perform infercnv operations to reveal cnv signal - step 3
  sample_output_step3 <- paste0(output_run_step3, sample)
  dir.create(sample_output_step3)
  bayes <- 0.25
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir= sample_output_step3,  # dir is auto-created for storing outputs
                               cluster_by_groups=T,   # cluster (Changed for uphyloplot2)
                               #plot_steps=T, # Changed from false for uphyloplot2
                               HMM_type='i6', # Changed for uphyloplot2
                               # uphyloplot2 recommends noise_filter=0.12 
                               denoise=T,
                               HMM=T,
                               analysis_mode = "samples",
                               num_threads = 6,
                               BayesMaxPNormal = bayes,
                               output_format = "pdf" # default: "png"
  )
  file.copy(from = paste0(sample_output_step3, "/infercnv.pdf"),
            to = paste0(output_run_select_step3,  sample, "_infercnv.pdf"),
            overwrite = T)
  file.copy(from = paste0(sample_output_step3, "/infercnv.20_HMM_predHMMi6.hmm_mode-samples.Pnorm_", 
                          bayes, ".repr_intensities.pdf"),
            to = paste0(output_run_select_step3,  sample, "_infercnv_bayes.pdf"),
            overwrite = T)
  
}



################################################################################
####### -------  Step 4 - uphyloplot2                ------- #######

for(sample in samples){
  
  # First, clear all results already in uphyloplot2
  system(paste0("rm ", path_uphyloplot2, "CNV_files/*"))
  system(paste0("rm ", path_uphyloplot2, "Inputs/*"))
  system(paste0("rm ", path_uphyloplot2, "output.svg"))
  
  file.copy(from = paste0(output_run_step4, sample, "_uphyloplot2.txt"),
            to = paste0(path_uphyloplot2, "Inputs/", sample, ".cell_groupings"),
            overwrite = T)
  
  # Run uphyloplot2
  setwd(path_uphyloplot2)
  system(paste0("python uphyloplot2.py -c 1" ))
  
  # Copy tree
  file.copy(to = paste0(output_run_step4, sample, "_uphyloplot2.svg"),
            from = paste0(path_uphyloplot2, "output.svg"),
            overwrite = T)
  
  # Copy key
  file.copy(to = paste0(output_run_step4, sample, "_uphyloplot2_key.csv"),
            from = paste0(path_uphyloplot2, "CNV_files/", sample, ".cell_groupings.csv"),
            overwrite = T)
  
}

################################################################################
####### -------          Annotate so with subclusters            ------- #######

# Add subclones for the samples that contain them

so_subclusters_meta <- so@meta.data
so_subclusters_meta$CNVsubclusters <- "none"

for(sample in samples){
  
  subclusters <- read.table(file = paste0(output_run_step4, sample, "_uphyloplot2.txt"),
              sep = "\t", header = T, stringsAsFactors = F)
  
  for(i in 1:nrow(subclusters)){
    so_subclusters_meta$CNVsubclusters[so_subclusters_meta$Sample == sample &
                                         rownames(so_subclusters_meta) == subclusters$cell[i]] <-
      paste0(sample, "_", subclusters$cell_group_name[i])
  }
  
}

# Add to so (should use CNA, not CNV)
so$CNAsubclusters <- so_subclusters_meta$CNVsubclusters

so_test <- subset(so, Sample == "005_Primary")
table(so_test$CNAsubclusters,
      so_test$normal_tumor)

so_test <- subset(so, Sample == "061_Primary")
table(so_test$CNAsubclusters,
      so_test$normal_tumor)

################################################################################
####### -------                   write so                      ------- #######

setwd("/Volumes/Partition 2/Core/Hayashi_Ewing_scRNAseq_Dec2022/")

saveRDS(so,
        "Data/Seurat_objects/AggrNoNorm_filtered_primary_celltype_CNAsubs.RDS")



