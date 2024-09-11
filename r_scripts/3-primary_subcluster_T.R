
################################################################################
####### -------                   Goals                          ------- #######

# Subtype T cells

################################################################################
####### -------        Read in data                              ------- #######

library(dplyr)
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
library(reshape2)
library(singleseqgset)
library(pheatmap)
library(openxlsx)

# Source custom functions
source("functions/FindNeighbors_FindClusters_PrintUMAP.R")
source("functions/SeuratQC.R")

setwd("/Volumes/Partition 2/Core/Hayashi_Ewing_scRNAseq_Dec2022/")
sobj <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_primary_celltype.RDS")

output <- "Analysis/Annotation/Primary/subcluster_Tcells/"
dir.create(output, recursive = T)
output_cell_type <- paste0(output, "cell_type/")
dir.create(output_cell_type, recursive = T)
output_cell_type_exploration <- paste0(output, "exploration_cell_type/")
dir.create(output_cell_type_exploration, recursive = T)

sobj_meta <- sobj@meta.data

annotation <- read.xlsx("helper_data/primary_T_cell_types.xlsx")
cluster_type <- as.list(annotation$Cell.Type)
names(cluster_type) <- annotation$Cluster
cluster_broad <- as.list(annotation$Broad.Cell.Type)
names(cluster_broad) <- annotation$Cluster

################################################################################
####### -------      Subset data        ------- #######

unique(sobj$Cell.Type)
cell_types <- c("T")
sobj_sub <- subset(sobj,
             Cell.Type %in% cell_types)
so <- sobj_sub

################################################################################
####### -------                 Normalize                        ------- #######

so <- NormalizeData(so, normalization.method = "LogNormalize")

################################################################################
####### -------               Add metrics                        ------- #######

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
so <- CellCycleScoring(so, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

################################################################################
####### -------               UMAP                               ------- #######

so <- FindVariableFeatures(so, selection.method = "vst")

so <- ScaleData(object = so, 
                vars.to.regress = c("nCount_RNA",
                                    "percent_mito"))
so <- RunPCA(so,
             features = VariableFeatures(so),
             verbose = F,
             ndims.print = 0)
ElbowPlot(so,
          ndims = 40)
n_dim <- 10
so <- RunUMAP(so, 
              reduction = "pca", 
              dims = 1:n_dim)

gg_sample <- DimPlot(so, 
                     group.by = "Sample",
                     reduction = "umap")
ggsave(paste0(output, "UMAP_sample.pdf"), gg_sample)

gg_sample_split <- DimPlot(so, 
                           group.by = "Sample",
                           reduction = "umap",
                           split.by = "Sample")
ggsave(paste0(output, "UMAP_sample_split.pdf"), gg_sample_split, height = 3, width = 15)

gg_celltype <- DimPlot(so, 
                     group.by = "Cell.Type",
                     reduction = "umap")
ggsave(paste0(output, "UMAP_CellType.pdf"), gg_celltype)


################################################################################
####### -------                   Clustering                     ------- #######

immune.combined <- so 

# Cluster, choose resolution
immune.combined <- FindNeighbors_FindClusters_PrintUMAP(sobj = immune.combined,
                                                 variable = "Sample",
                                                 neighbors_dim = 1:15,
                                                 resolutionsCalculated = c(0.1,
                                                0.2, 0.3, 0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1))
resolutionChosen <- 0.9
assay <- DefaultAssay(immune.combined)
immune.combined$seurat_clusters <- immune.combined@meta.data[,paste(assay, "_snn_res.", resolutionChosen, sep = "")] 

Idents(object = immune.combined) <- "seurat_clusters"
so <- immune.combined
# Compare distribution of new clusters and samples
prop_df_sample <- data.frame(prop.table(table(so@meta.data[,"Sample"], 
                                              so@meta.data[,"seurat_clusters"]), margin = 2))
gg_stackbar_sample <- ggplot(prop_df_sample, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="stack", stat="identity") +ylab("Proportion") + xlab("New cluster") +
  theme_classic()
gg_stackbar_sample
ggsave(paste0(output, "/", "Sample_stacked_barplot.pdf"),
       gg_stackbar_sample)

################################################################################
####### -------              Heterogeneity analysis              ------- #######

clustering_resolutions <- c("seurat_clusters")
n_samples_positive <- 2 # Number of samples that must have a gene positive FC

all = msigdbr(species = "Homo sapiens")

other_gene_sets_literature <- c(
  "GOBP_CD8_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION" ,
  "GOBP_CD8_POSITIVE_ALPHA_BETA_T_CELL_DIFFERENTIATION",
  "GOBP_CD8_POSITIVE_ALPHA_BETA_T_CELL_PROLIFERATION",                                             
  "GOBP_NEGATIVE_REGULATION_OF_CD8_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION",                         
  "GOBP_POSITIVE_REGULATION_OF_CD8_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION"     ,                    
  "GOBP_POSITIVE_REGULATION_OF_CD8_POSITIVE_ALPHA_BETA_T_CELL_DIFFERENTIATION"    ,                
  "GOBP_REGULATION_OF_CD8_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION"             ,                     
  "GOBP_REGULATION_OF_CD8_POSITIVE_ALPHA_BETA_T_CELL_DIFFERENTIATION",
  "GOBP_CD4_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION"            ,                                
  "GOBP_CD4_POSITIVE_ALPHA_BETA_T_CELL_CYTOKINE_PRODUCTION"     ,                              
 "GOBP_CD4_POSITIVE_ALPHA_BETA_T_CELL_DIFFERENTIATION"        ,                               
 "GOBP_CD4_POSITIVE_ALPHA_BETA_T_CELL_PROLIFERATION"       ,                                  
 "GOBP_CD4_POSITIVE_CD25_POSITIVE_ALPHA_BETA_REGULATORY_T_CELL_DIFFERENTIATION"   ,           
 "GOBP_CD4_POSITIVE_OR_CD8_POSITIVE_ALPHA_BETA_T_CELL_LINEAGE_COMMITMENT"      ,              
 "GOBP_NEGATIVE_REGULATION_OF_CD4_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION"     ,                
 "GOBP_NEGATIVE_REGULATION_OF_CD4_POSITIVE_ALPHA_BETA_T_CELL_DIFFERENTIATION"  ,              
 "GOBP_NEGATIVE_REGULATION_OF_CD4_POSITIVE_ALPHA_BETA_T_CELL_PROLIFERATION"  ,                
 "GOBP_POSITIVE_REGULATION_OF_CD4_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION"    ,                 
 "GOBP_POSITIVE_REGULATION_OF_CD4_POSITIVE_ALPHA_BETA_T_CELL_DIFFERENTIATION" ,               
 "GOBP_POSITIVE_REGULATION_OF_CD4_POSITIVE_ALPHA_BETA_T_CELL_PROLIFERATION"  ,                
 "GOBP_REGULATION_OF_CD4_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION"      ,                        
 "GOBP_REGULATION_OF_CD4_POSITIVE_ALPHA_BETA_T_CELL_DIFFERENTIATION"    ,                     
 "GOBP_REGULATION_OF_CD4_POSITIVE_CD25_POSITIVE_ALPHA_BETA_REGULATORY_T_CELL_DIFFERENTIATION"
)

other_hallmarks <- unique(all$gs_name[grep("HALLMARK_", all$gs_name, ignore.case = T)])

################################################################################
####### -------      Add EMT and EWS scores           ------- #######

DefaultAssay(so) <- "RNA"
all = msigdbr(species = "Homo sapiens")
gene_sigs <- list()
for(geneset in unique(c(other_gene_sets_literature, other_hallmarks))){
  print(geneset)
  gene_sigs[[paste0(geneset, "_activity")]] <- all$gene_symbol[all$gs_name == geneset]
}

length(gene_sigs)

so <- AddModuleScore(
  so,
  features = gene_sigs,
  name = names(gene_sigs)
)

# Replace the '1' (and other numbers) added to the signature names
colnames(so@meta.data) <- gsub("activity\\d+$", "activity", colnames(so@meta.data))

################################################################################
####### -------           Add cell type labels                   ------- #######

# Update the top section before running below

# Discrete cell type
so$Cell.Type <- "cell"
for(cluster in unique(names(cluster_type))){
  so$Cell.Type[so$seurat_clusters == cluster] <- cluster_type[[cluster]]
}
# Broad cell type
so$Broad.Cell.Type <- "cell"
for(cluster in unique(names(cluster_broad))){
  so$Broad.Cell.Type[so$seurat_clusters == cluster] <- cluster_broad[[cluster]]
}

################################################################################
####### -------                 Cell type UMAP                   ------- #######

umap_cell <- DimPlot(so, reduction = "umap",
                     group.by = "Cell.Type",
                     label = T,
                     label.size = 5) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(paste0(output_cell_type, "UMAP_CellType.PDF"), 
       plot = umap_cell,
       width = 10,
       height = 7.5)

################################################################################
####### -------               Cell type dotplot                  ------- #######

ggdot_celltype_final <- DotPlot(so,
                                features = c("CD8A", 
                                             "CD8B", 
                                             "KLRK1", # effector
                                             "KLRB1", # Memory
                                             "CD4",
                                             # Helper
                                             "MAF",
                                             # Regulatory:
                                             "FOXP3", "CTLA4"
                                             
                                ),
                                cluster.idents=T,
                                group.by = "Cell.Type") + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
ggdot_celltype_final

ggsave(paste0(output_cell_type_exploration, "Dotplot_CellType.PDF"), 
       plot = ggdot_celltype_final)


so$cluster.Cell.Type <- paste0(so$seurat_clusters, ".", so$Cell.Type)
ggdot_celltype_extra <- DotPlot(so,
                                features = c("CD8A", 
                                             "CD8B", 
                                             "KLRK1", # effector
                                             "KLRB1", # Memory
                                             "CD4",
                                             # Helper
                                             "MAF",
                                             # Regulatory:
                                             "FOXP3", "CTLA4",
                                             # Added
                                             "CD3D", # Pan-T
                                             "BATF", "TIGIT", # Treg
                                             "CCL4", "CCL5", "NKG7", "KLRD1", "PRF1", # CD8 T
                                             "PDCD1", "LAG3", #Negative
                                             "KLRC1",
                                             "IL7R"
                                ),
                                cluster.idents=T,
                                group.by = "cluster.Cell.Type") + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
ggdot_celltype_extra
ggsave(paste0(output_cell_type_exploration, "Dotplot_CellType_extra.PDF"), 
       plot = ggdot_celltype_extra)


################################################################################
####### -------     Barplot of cell type by sample               ------- #######

so$Sample_n <- so$Sample
for(i in 1:nrow(so@meta.data)){
  print(i)
  
  so$Sample_n[i] <- paste0(so$Sample_n[i], " (n=", 
                           sum(so$Sample == so$Sample_n[i]), ")")
}

# Compare distribution of new clusters and samples
prop_sample_celltype <- data.frame(prop.table(table(so@meta.data[,"Cell.Type"], 
                                                    so@meta.data[,"Sample_n"]), margin = 2))
gg_stackbar_sample <- ggplot(prop_sample_celltype, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="stack", stat="identity") +ylab("Cell type proportion") + 
  xlab("") + labs(fill="Cell type") +
  theme_classic() + theme(axis.text.x=element_text(angle=90,
                                                   hjust=1, vjust = 0.5, size = 12))
gg_stackbar_sample
ggsave(paste0(output_cell_type, "/", "Celltype_stacked_barplot.pdf"),
       gg_stackbar_sample)

################################################################################
####### -------            Final cell type dotplot               ------- #######

ggdot_Jul_w_checkpoint <- DotPlot(so,
                     features = c(c("CD3D", "CD3E"),
                                  c("CD8A", "CD8B", "CD4"),
                                  c("KLRK1", "KLRD1", "CD7"), # Effector
                                  c("CD69", "KLF6", "KLRB1"), # CD161 - more memory-related
                                  c("IL7R", "CD40LG", "MAF"), # CD4+IL7R+
                                  c("FOXP3", "CTLA4", "TIGIT", "BATF"), # Treg
                                  c( "LAG3", "PDCD1", "TOX"), #Exhaustion markers
                                  c("CD28", "CLEC2D"), # Related to TME
                                  c("MIF", "CD74", "CXCR4") # cellchat finding
                     ),
                     cluster.idents=T,
                     group.by = "Cell.Type") + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
ggdot_Jul_w_checkpoint
ggsave(paste0(output_cell_type, "Dotplot_cell_type_w_checkpoint.pdf"), 
       plot = ggdot_Jul_w_checkpoint)


################################################################################
####### -------               Cell type featureplot              ------- #######

# https://github.com/satijalab/seurat/issues/1080
# List of featureplots
featureplot_celltype_final <- FeaturePlot(so,
                                          features = c(c("CD3D", "CD3E"),
                                                       c("CD8A", "CD8B", "CD4"),
                                                       c("KLRK1", "KLRD1", "CD7"),
                                                       c("CD69", "KLF6", "KLRB1" # CD161
                                                       ),
                                                       c("IL7R", "CD40LG", "MAF"),
                                                       c("FOXP3", "CTLA4", "TIGIT", "BATF"), #,
                                                       c("LAG3", "PDCD1"),
                                                       c("CD28", "CLEC2D"), # Related to TME
                                                       c("MIF", "CD74", "CXCR4") # cellchat finding
                                          ),
                                          combine = F)
for(i in 1:length(featureplot_celltype_final)) {
  featureplot_celltype_final[[i]] <- featureplot_celltype_final[[i]] + 
    NoLegend() + NoAxes() + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
}
featureplot_celltype_final_grid <- cowplot::plot_grid(plotlist = featureplot_celltype_final)

ggsave(paste0(output_cell_type, "UMAP_celltype_marker_grid.pdf"), 
       plot = print(featureplot_celltype_final_grid),
       height = 8,
       width = 10)


################################################################################
####### -------       DGE and targeted ORA by cell type          ------- #######

####### -------       Output
output_ORA <- paste0(output, "/ORA_cell_type/")
dir.create(output_ORA)

####### -------       subset
so_sub <- so
Idents(so_sub) <- so_sub$Cell.Type
table(Idents(so_sub))

####### -------       ORA-based method
so.markers_orig <- FindAllMarkers(so_sub, 
                                  only.pos = TRUE, 
                                  min.pct = 0.05, 
                                  logfc.threshold = 0.1)
so.markers <- so.markers_orig[so.markers_orig$p_val_adj < 0.05,] # Only keep sig genes
table(so.markers$cluster) 

# Writing the sig tables
write.table(so.markers,
            file=paste0(output_ORA,
                        "Cell_type_FindAllMarkers_sigOnly.txt"),
            sep = "\t",
            row.names = F)

# Subsetting to top genes
so.markers_ORA <- so.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC) 

# Gene sets
all = msigdbr(species = "Homo sapiens") 
specific_output1 <- "ORA_C8"
gene_sets1 <- all[all$gs_cat %in% "C8",] %>% 
  dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
length(unique(gene_sets1$gs_name))

# Testing out PangloaDB (plus with technical gene sets)
panglao_gene_set_plus <- read.table(file = paste0("helper_data/PanglaoDB_plus_human_GeneSets.csv"),
                                    sep = ",", header = T, stringsAsFactors = F)
panglao_gene_set_plus$Collection_Shortname <- "PanglaoDB_plus"
panglao_gene_set_plus$Gene_Set_Shortname <- panglao_gene_set_plus$gs_name
gene_sets2 = panglao_gene_set_plus[,c("Gene_Set_Shortname", "gene_symbol")]
specific_output2 <- "PanglaoDB_plus"

# Cycle through each cluster and perform ORA
specific_outputs <- c(specific_output1,
                      specific_output2)
gene_sets_list <- list(gene_sets1,
                       gene_sets2) 
for(i in 1:length(specific_outputs)){
  gene_sets <- gene_sets_list[[i]]
  specific_output <- specific_outputs[i]
  
  clusters <- unique(so.markers_ORA$cluster)
  background <- rownames(so_sub) # All gene names in dataset
  output_ORA_collection <- paste0(output_ORA, specific_output, "/")
  dir.create(output_ORA_collection,
             recursive = T)
  
  for(cluster in clusters){
    print(cluster)
    
    genes <- so.markers_ORA$gene[so.markers_ORA$cluster %in% cluster]
    
    # Gene set 1
    ora_enricher <- enricher(genes, TERM2GENE=gene_sets, pAdjustMethod = "fdr", 
                             minGSSize = 1, maxGSSize = 1000, 
                             qvalueCutoff = 1, pvalueCutoff =1,
                             universe = background) 
    ora_result <- data.frame(ora_enricher)
    write.table(ora_result,
                file = paste0(output_ORA_collection,
                              cluster, ".txt"),
                row.names = F, sep = "\t", quote = F)
  }
}

################################################################################
####### -------              Save RDS object                     ------- #######

saveRDS(so,
        "Data/Seurat_objects/AggrNoNorm_filtered_primary_celltype_CNAsubs_sublusteredTcells_pathway.RDS")

