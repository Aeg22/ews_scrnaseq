

################################################################################
####### -------                     Set paths                    ------- #######

setwd("/Volumes/Partition 2/Core/Hayashi_Ewing_scRNAseq_Dec2022/")
output <- "GEO/"
dir.create(output)

library(Seurat)
library(ggplot2)
library(plyr)

################################################################################
####### -------        Read in data (merge if necessary)         ------- #######

# More complicated because there are two main seurat objects (primary and CTC)
# and have metadata across multiple seurat objects that I want to combine

####### -------        CTC
ctc <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_celltype_CTC.RDS")
DefaultAssay(ctc) <- "RNA" # Only RNA has the counts
colnames(ctc@meta.data) # Will remove some columns 

####### -------        Primary
# This object should contain all primary cells
so <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_primary_celltype.RDS")
DefaultAssay(so) <- "RNA" # Only RNA has the counts
colnames(so@meta.data) # Will remove some columns later
so_meta <- so@meta.data # Will add subsequent metadata columns to this
so_meta <- cbind.data.frame(cellbarcode = rownames(so_meta), so_meta)
orig_n <- ncol(so_meta) # Any columns in addition to this will be added back to seurat object

# Integrated tumor w/ CNAsubclusters and NMF
integrated_tumor <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_primary_celltype_CNAsubs_integrateP_clustered_pathway_NMF.RDS")
colnames(integrated_tumor@meta.data)
integrated_tumor_meta <- integrated_tumor@meta.data
cols_to_add <- c("CNAsubclusters", "Merge.GEP.EWS",                                                  
                 "Merge.GEP.OxPhos",                                                  
                 "Merge.GEP.UPR",                                                     
                 "Merge.GEP.Proliferation",                                           
                 "Merge.GEP.MAPK",                                                    
                 "Merge.GEP.Ribosome",
                 "merged_program")
integrated_tumor_meta <- cbind.data.frame(cellbarcode = rownames(integrated_tumor_meta),
                                          integrated_tumor_meta[,cols_to_add])
umapCoord <- as.data.frame(Embeddings(object = integrated_tumor[["umap"]]))
integrated_tumor_meta <- cbind(integrated_tumor_meta,
                               umapCoord)
colnames(integrated_tumor_meta)
colnames(integrated_tumor_meta)[10:11] <- paste0("Tumor_integrated_", 
                                                 colnames(integrated_tumor_meta)[10:11])
so_meta <- join(so_meta, integrated_tumor_meta)

# All immune subtype labels (not going to include their UMAP)
immune_subtypes <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_primary_celltype_specific_immune.RDS")
colnames(immune_subtypes@meta.data)
table(immune_subtypes$Sub.Cell.Type)
immune_subtypes_meta <- immune_subtypes@meta.data
cols_to_add <- c("Sub.Cell.Type")
immune_subtypes_meta <- cbind.data.frame(cellbarcode = rownames(immune_subtypes_meta),
                                          immune_subtypes_meta[,cols_to_add])
colnames(immune_subtypes_meta) <- c("cellbarcode", cols_to_add)
so_meta <- join(so_meta, immune_subtypes_meta)

# SCEVAN metadata
scevan <- read.table(file = "Analysis/Annotation/Primary/Normalized/SCEVAN/multiSample/SCEVAN_metadata.pdf",
                     sep = "\t", header = T, stringsAsFactors = F)
scevan <- cbind.data.frame(cellbarcode = scevan$barcodes,
                           scevan[,"SCEVAN_class"])
colnames(scevan) <- c("cellbarcode", "SCEVAN_class")
so_meta <- join(so_meta, scevan)
colnames(so_meta)

# Add new metadata back to seurat object (will still remove some columns later)
for(i in (orig_n+1):ncol(so_meta)){
  print(i)
  so@meta.data[,colnames(so_meta)[i]] <- so_meta[,i]
}
colnames(so@meta.data)

# Integrated primary:
# AggrNoNorm_filtered_primary_celltype_CNAsubs_integrateP_clustered_pathway_NMF.RDS

# All immune subtype labels??
# AggrNoNorm_filtered_primary_celltype_specific_immune.RDS

# Myeloid subset
# AggrNoNorm_filtered_primary_celltype_CNAsubs_sublusteredMyeloid_pathway.RDS

# T subset
# AggrNoNorm_filtered_primary_celltype_CNAsubs_sublusteredTcells_pathway.RDS


################################################################################
####### -------       Write counts matrix        ------- #######

# This one is complicated because there are two main seurat objects (primary and CTC)

####### -------        CTC
counts_matrix <- as.data.frame(as.matrix(GetAssayData(object = ctc, slot = "counts")))
# Above, as.data.frame() is necessary over data.frame()
counts_matrix <- cbind.data.frame(Gene = rownames(counts_matrix),
                                  counts_matrix)
dim(counts_matrix)
counts_matrix_sub <- counts_matrix[1:10, 1:10]
View(counts_matrix_sub)
write.table(counts_matrix, 
            paste0(output, "counts_matrix_blood.txt"), 
            sep = '\t', 
            row.names = F, 
            col.names = T, 
            quote = F)
# Confirm file works and is complete (compare to counts_matrix)
counts_matrix_ReadIn <- read.table(
            file = paste0(output, "counts_matrix_blood.txt"), 
            sep = '\t', 
            header = T,
            stringsAsFactors = F,
            quote = "")
all.equal(counts_matrix,counts_matrix_ReadIn) # Attributes only are fine
which(counts_matrix != counts_matrix_ReadIn, arr.ind=TRUE) # Shows if there are any differences, row col means no differences

####### -------        Primary
counts_matrix <- as.data.frame(as.matrix(GetAssayData(object = so, slot = "counts")))
counts_matrix <- cbind.data.frame(Gene = rownames(counts_matrix),
                                  counts_matrix)
dim(counts_matrix)
counts_matrix_sub <- counts_matrix[1:10, 1:10]
View(counts_matrix_sub)
write.table(counts_matrix, 
            paste0(output, "counts_matrix_primary.txt"), 
            sep = '\t', 
            row.names = F, 
            col.names = T, 
            quote = F)
# Confirm file works and is complete (compare to counts_matrix)
counts_matrix_ReadIn <- read.table(
  file = paste0(output, "counts_matrix_primary.txt"), 
  sep = '\t', 
  header = T,
  stringsAsFactors = F)
all.equal(counts_matrix,counts_matrix_ReadIn) # Attributes only are fine
which(counts_matrix != counts_matrix_ReadIn, arr.ind=TRUE) # Shows if there are any differences, row col means no differences


################################################################################
####### -------                    metadata                      ------- #######

####### -------        CTC
# Add umap coordinates to metadata
umapCoord <- as.data.frame(Embeddings(object = ctc[["umap"]]))
metadata <- ctc@meta.data
metadata <- cbind(metadata,
                  umapCoord)
colnames(metadata)
metadata$orig.ident
remove <- c(colnames(metadata)[grep("RNA_snn_res", colnames(metadata))],
            "orig.ident", "cluster_full")
metadata_sub <- metadata[,!(colnames(metadata) %in% remove)]
# Add Cell
metadata_sub <- cbind.data.frame(Cell = rownames(metadata_sub),
                                 metadata_sub)
write.table(metadata_sub, 
            paste0(output, "cell_metadata_blood.txt"), 
            sep = '\t', 
            row.names = F, 
            col.names = T, 
            quote = F)
plot(metadata_sub$UMAP_1,
     metadata_sub$UMAP_2) # Just confirm I have the correct UMAP coordinates
metadata_sub_ReadIn <- read.table(paste0(output, "cell_metadata_blood.txt"), 
                                  sep = '\t', header = T)
all.equal(metadata_sub,metadata_sub_ReadIn) # Attributes only are fine
which(metadata_sub != metadata_sub_ReadIn, arr.ind=TRUE) # Shows if there are any differences, row col means no differences


####### -------        Primary
# Add umap coordinates to metadata
umapCoord <- as.data.frame(Embeddings(object = so[["umap"]]))
metadata <- so@meta.data
metadata <- cbind(metadata,
                  umapCoord)
colnames(metadata)
remove <- c(colnames(metadata)[grep("RNA_snn_res", colnames(metadata))],
            "orig.ident", "cluster_full")
metadata_sub <- metadata[,!(colnames(metadata) %in% remove)]
# Add Cell
metadata_sub <- cbind.data.frame(Cell = rownames(metadata_sub),
                                 metadata_sub)
write.table(metadata_sub, 
            paste0(output, "cell_metadata_primary.txt"), 
            sep = '\t', 
            row.names = F, 
            col.names = T, 
            quote = F)
plot(metadata_sub$UMAP_1,
     metadata_sub$UMAP_2) # Just confirm I have the correct UMAP coordinates
metadata_sub_ReadIn <- read.table(paste0(output, "cell_metadata_primary.txt"), 
                                  sep = '\t', header = T, stringsAsFactors = F)
all.equal(metadata_sub,metadata_sub_ReadIn) # Attributes only are fine
which(metadata_sub != metadata_sub_ReadIn, arr.ind=TRUE) # Shows if there are any differences, row col means no differences

