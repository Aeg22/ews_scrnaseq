
################################################################################
####### -------                   Goals                          ------- #######

# Read in cellranger aggr outs/
# Make a seurat object
# Add sample names to metadata

################################################################################
####### -------        Notes on this analysis                    ------- #######

# Normalization was not performed in cellranger aggr

################################################################################
####### -------        Read in R packages                        ------- #######

library(dplyr)
library(plyr)
library(Seurat)

################################################################################
####### -------        Set paths and parameters                  ------- #######

setwd("/Volumes/Partition 2/Core/Hayashi_Ewing_scRNAseq_Dec2022/") # Main project directory

cellranger_outs_path <- "/Volumes/Partition 2/Core/Hayashi_Ewing_scRNAseq_Aug2021/Data/Aggr_noNorm_12-7-21/" # Location of the cellranger outs directory from the server
seurat_outputs <- "Data/Seurat_objects/"
dir.create(seurat_outputs,
           recursive = T)

min.cells <- 10
min.features <- 100

################################################################################
####### -------                Read in data                      ------- #######

# Load the dataset
data.10x <- Read10X(data.dir = paste0(cellranger_outs_path,
                                       "/count/filtered_feature_bc_matrix/")
)
# The contents of `data.dir` are the filtered_feature_bc_* supplemental files supplied to GEO

# Initialize the Seurat object with the raw (non-normalized data)
so <- CreateSeuratObject(counts = data.10x, 
                           min.cells = min.cells, 
                           min.features = min.features)

################################################################################
####### -------        Add sample names to metadata              ------- #######

# Samples are labels -1, -2, etc., which reflects the order of the samples in the CSV used in cellranger aggr
# Read in this CSV to add the sample info
aggr <- read.table(file = paste0(cellranger_outs_path,
                                 "/aggregation.csv"),
                   sep = ",", header = T, stringsAsFactors = F)
aggr$Sample_num <- 1:nrow(aggr)

barcode_order <- cbind.data.frame(Full = rownames(so@meta.data),
                                  Sample_num = sapply(strsplit(rownames(so@meta.data),
                                                               split = "-"), function(x){
                                                                 x[[2]]
                                                               }))

barcode_order <- join(barcode_order,
                      aggr)

# Add the sample to the metadata
so <- AddMetaData(so, 
                  barcode_order$sample_id, 
                  col.name = "Sample")

# Clean up sample names if necessary
so$Sample <- gsub("_sample", "", so$Sample)
so$Sample <- gsub("_ReSeqDec2021", "", so$Sample)


# Write RDS object
saveRDS(so, paste0(seurat_outputs,
                     "/AggrNoNorm_minFiltered.RDS"))

