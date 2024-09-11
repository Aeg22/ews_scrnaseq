

################################################################################
####### -------                    Goals                         ------- #######

# Processing the primary samples only

################################################################################
####### -------                 Load R packages                  ------- #######

library(Seurat)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(msigdbr)
library(plyr)
library(dplyr)
library(clusterProfiler)

source("functions/SeuratQC.R")

################################################################################
####### -------                Update variables                  ------- #######

setwd("/Volumes/Partition 2/Core/Hayashi_Ewing_scRNAseq_Dec2022/")

input_rds = "Data/Seurat_objects/AggrNoNorm_filtered.RDS"
samples <- c("051_Primary", "061_Primary", "066_Primary",
             "005_Primary", "010_Primary", "038_Primary", "001_Primary")
output <- paste0("Analysis/Annotation/Primary/")
dir.create(output, recursive = T)

output_full <- paste0(output, "Normalized/")
dir.create(output_full)

# Cluster identities - Only known after the cell type identity section (skip in the beginning)

library(openxlsx)
annotation <- read.xlsx("helper_data/primary_cell_types.xlsx")
cluster_type <- as.list(annotation$Cell.Type)
names(cluster_type) <- annotation$Cluster
cluster_broad <- as.list(annotation$Broad.Cell.Type)
names(cluster_broad) <- annotation$Cluster


################################################################################
####### -------           Read in and subset data                ------- #######

seurat_obj <- readRDS(input_rds)
so <- subset(seurat_obj, Sample %in% samples)
DefaultAssay(so) <- "RNA"

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
so <- RunPCA(so, features = VariableFeatures(so),
             verbose = F, ndims.print = 0)
ElbowPlot(so, ndims = 40)
n_dim <- 30 
so <- RunUMAP(so, reduction = "pca", dims = 1:n_dim)

gg_sample <- DimPlot(so, group.by = "Sample", reduction = "umap")
ggsave(paste0(output_full, "UMAP_sample.pdf"),
       gg_sample)

gg_sample_split <- DimPlot(so, 
                     group.by = "Sample",
                     reduction = "umap",
                     split.by = "Sample")
ggsave(paste0(output_full, "UMAP_sample_split.pdf"),
       gg_sample_split, height = 3, width = 15)

################################################################################
####### -------               UMAP with QC metrics               ------- #######

output_umap_qc <- paste0(output_full, "/qc_umap/")
dir.create(output_umap_qc, recursive = T)

gg_MT <- FeaturePlot(so, 
                     features = "percent_mito",
                     reduction = "umap")
ggsave(paste0(output_umap_qc, "UMAP_MT.pdf"), gg_MT)

gg_CCsscore <- FeaturePlot(so, 
                           features = "S.Score",
                           reduction = "umap")
ggsave(paste0(output_umap_qc, "UMAP_CCsscore.pdf"), gg_CCsscore)

gg_CCg2m <- FeaturePlot(so, 
                           features = "G2M.Score",
                           reduction = "umap")
ggsave(paste0(output_umap_qc, "UMAP_CCg2Mscore.pdf"), gg_CCg2m)

gg_phase <- DimPlot(so, 
                    group.by = "Phase",
                    reduction = "umap")
ggsave(paste0(output_umap_qc, "UMAP_CCPhase.pdf"), gg_phase)

gg_umi <- FeaturePlot(so, 
                      features = "nCount_RNA",
                      reduction = "umap")
ggsave(paste0(output_umap_qc, "UMAP_UMIs.pdf"), gg_umi)


################################################################################
####### -------              Clustering                         ------- #######

assay <- "RNA"
DefaultAssay(so) <- assay

so <- FindNeighbors(so, 
                    reduction = "pca",
                    dims = 1:n_dim, 
                    k.param = 20) 

resolutionsCalculated <- c( 0.1,  0.3,  0.5, 1, 1.5, 1.75, 2)
so <- FindClusters(so, 
                   resolution = resolutionsCalculated) 

for (i in resolutionsCalculated){
  # set resolution
  Idents(object = so) <- paste(assay, "_snn_res." , i, sep = "") # name is picked from column name of resolutions
  
  # general umap
  allClusters <- DimPlot(so,
                         reduction = "umap",
                         label = TRUE,
                         label.size = 6) + ggtitle(paste(assay, "_snn_res." , i, sep = "")) + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
  ## UMAP of cells in each cluster by sample
  splitBySample <- DimPlot(so, 
                           label = TRUE, 
                           split.by = "Sample")  + NoLegend() + ggtitle("Clusters split by sample") + theme(plot.title = element_text(hjust = 0.5))
  
  # Explore whether clusters segregate by cell cycle phase
  splitByPhase <- DimPlot(so,
                          label = TRUE, 
                          split.by = "Phase")  + NoLegend() + ggtitle("Clusters split by cell cycle phase") + theme(plot.title = element_text(hjust = 0.5))
  grid.arrange(allClusters, splitBySample
               #splitByPhase 
               + remove("x.text"), 
               ncol = 2, nrow = 1)
  
}

# Choose a resolution
resolutionChosen <- 1.75 

# tells Seurat which resolution to base the remainder of the analysis
Idents(object = so) <- paste(assay, "_snn_res.", resolutionChosen, sep = "") # name is picked from column name of resolutions
so$seurat_clusters <- so@meta.data[,paste(assay, "_snn_res.", resolutionChosen, sep = "")] # Only necessary since some of my script references seurat_clusters
so$cluster_full <- so$seurat_clusters

umap_cluster <- DimPlot(so, reduction = "umap",
                        label = T) # labels by cluster
ggsave(paste0(output_full, "UMAP_Cluster.PDF"), plot = umap_cluster)

df_cluster_distribution <- data.frame(unclass(table(so$seurat_clusters)))
df_cluster_distribution <- cbind.data.frame(Cluster = rownames(df_cluster_distribution),
                                            Cells = df_cluster_distribution[,1])
write.table(df_cluster_distribution,
            file = paste0(output_full, "Cluster_distribution.txt"),
            sep = "\t", row.names = F, quote = F)

# Compare distribution of new clusters and samples
prop_df_sample <- data.frame(prop.table(table(so@meta.data[,"Sample"], 
                                              so@meta.data[,"seurat_clusters"]), margin = 2))
gg_stackbar_sample <- ggplot(prop_df_sample, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="stack", stat="identity") +ylab("Proportion") + xlab("New cluster") +
  theme_classic()
gg_stackbar_sample
ggsave(paste0(output_full, "/", "Sample_stacked_barplot.pdf"),
       gg_stackbar_sample)


################################################################################
####### -------         Cluster-level QC (and by sample)         ------- #######

output_qc_cluster <- paste0(output_full, "/qc_by_cluster/")
dir.create(output_qc_cluster)

# Whole dataset
SeuratQC(so = so,
         category = "seurat_clusters",
         outdir = output_qc_cluster,
         violin_pt_size = 0.01)

# By sample
for(sample in unique(so$Sample)){
  print(sample)
  temp_dir <- paste0(output_qc_cluster, sample, "/")
  dir.create(temp_dir)
  so_temp <- subset(so, Sample == sample)
  SeuratQC(so = so_temp,
           category = "seurat_clusters",
           outdir = temp_dir,
           violin_pt_size = 0.01)
}

################################################################################
####### -------       Help with determining cell type identity   ------- #######

####### -------       ORA-based method

Idents(so) <- so$seurat_clusters
so.markers_orig <- FindAllMarkers(so, 
                                  only.pos = TRUE, 
                                  min.pct = 0.1, 
                                  logfc.threshold = 0.25)
so.markers <- so.markers_orig[so.markers_orig$p_val_adj < 0.05,] # Only keep sig genes
table(so.markers$cluster) 
output_ORA <- paste0(output_full, "/ORA/")
dir.create(output_ORA)
write.table(so.markers,
            file=paste0(output_ORA,
                        "Clusters_FindAllMarkers_sigOnly.txt"),
            sep = "\t",
            row.names = F)

# Subsetting to top genes
so.markers_ORA <- so.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC) # 

# Gene sets
all = msigdbr(species = "Homo sapiens") 
specific_output <- "ORA_C8"
gene_sets <- all[all$gs_cat %in% "C8",] # C8 only
length(unique(gene_sets$gs_name))
gene_sets <- gene_sets %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

# Cycle through each cluster and perform ORA
clusters <- unique(so.markers_ORA$cluster)
background <- rownames(so) # All gene names in dataset
output_ORA_collection <- paste0(output_ORA, specific_output, "/")
dir.create(output_ORA_collection,
           recursive = T)
for(cluster in clusters){
  print(cluster)
  
  genes <- so.markers_ORA$gene[so.markers_ORA$cluster %in% cluster]
  
  ora_enricher <- enricher(genes, TERM2GENE=gene_sets, pAdjustMethod = "fdr", 
                           minGSSize = 5, maxGSSize = 1000, 
                           qvalueCutoff = 1, pvalueCutoff =1,
                           universe = background) 
  ora_result <- data.frame(ora_enricher)
  
  write.table(ora_result,
              file = paste0(output_ORA_collection,
                            cluster, ".txt"),
              row.names = F, sep = "\t", quote = F)
}


# Testing out PangloaDB (plus with technical gene sets)
panglao_gene_set_plus <- read.table(file = paste0("helper_data/PanglaoDB_plus_human_GeneSets.csv"),
                                    sep = ",", header = T, stringsAsFactors = F)
panglao_gene_set_plus$Collection_Shortname <- "PanglaoDB_plus"
panglao_gene_set_plus$Gene_Set_Shortname <- panglao_gene_set_plus$gs_name
gene_sets_to_analyze = panglao_gene_set_plus[,c("Gene_Set_Shortname", "gene_symbol")]

specific_output <- "PanglaoDB_plus"
output_ORA_collection <- paste0(output_ORA, specific_output, "/")
dir.create(output_ORA_collection,
           recursive = T)
for(cluster in clusters){
  print(cluster)
  
  genes <- so.markers_ORA$gene[so.markers_ORA$cluster %in% cluster]
  
  ora_enricher <- enricher(genes, TERM2GENE=gene_sets_to_analyze, pAdjustMethod = "fdr", 
                           minGSSize = 5, maxGSSize = 1000, 
                           qvalueCutoff = 1, pvalueCutoff =1,
                           universe = background) 
  ora_result <- data.frame(ora_enricher)
  
  write.table(ora_result,
              file = paste0(output_ORA_collection,
                            cluster, ".txt"),
              row.names = F, sep = "\t", quote = F)
}


################################################################################
####### -------       Select markers dotplot by cluster          ------- #######

output_dotplot_celltype <- paste0(output_full, "/dotplot_celltype/")
dir.create(output_dotplot_celltype)

output_dotplot_celltype_exploration <- paste0(output_full, "/dotplot_celltype_exploration/")
dir.create(output_dotplot_celltype_exploration)

ggdot_celltype_clustered_update <- DotPlot(so,
                                    features = c("MAPT", "KDSR", "NKX2-2", "NR0B1", "STEAP1", "PRKCB", # EWS
                                                 "MKI67",  "PCNA", # Proliferating
                                                 "MYL9", "ACTA2", # Smooth muscle
                                                 "JCHAIN", "MS4A1", # B cell
                                                 "CD96", "CD247", # General T cell
                                                 "CD4", "RORA", # CD4+ T cells
                                                 "CD8A", "CD8B", # CD8+ T cells
                                                 "CD34", "PECAM1", # Endothelial
                                                 "PDGFRA", "MFAP5", # Fib
                                                 "CA2" , "SPP1", # Osteoclasts
                                                 "CD14", "CSF1R" # Myeloid
                                                
                                    ),
                                    group.by = "seurat_clusters",
                                    cluster.idents = T) + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
ggdot_celltype_clustered_update

ggsave(paste0(output_dotplot_celltype_exploration, "Dotplot_celltype_clustered_UPDATE.PDF"), 
       plot = ggdot_celltype_clustered_update)

################################################################################
####### -------                Name cells                        ------- #######

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

tumor_cells <- c("Ewing sarcoma", "Ewing sarcoma proliferating")

################################################################################
####### -------          Choose order of cells                   ------- #######

unique(so$Cell.Type)
so$Cell.Type <- factor(so$Cell.Type,
                       levels = c("Ewing sarcoma",
                                  "Ewing sarcoma proliferating",
                                  "Inflammatory CAF",
                                  "Matrix CAF",
                                  "Myeloid",
                                  "Osteoclast",
                                  "Endothelial",
                                  "T",
                                  "B"))

unique(so$Broad.Cell.Type)
so$Broad.Cell.Type <- factor(so$Broad.Cell.Type,
                       levels = c("Ewing sarcoma",
                                  "",
                                  "Fibroblast",
                                  "Myeloid",
                                  "Endothelial",
                                  "T",
                                  "B",
                                  "Osteoclast"))

################################################################################
####### -------                 Cell type UMAP                   ------- #######

umap_cell <- DimPlot(so, reduction = "umap",
                     group.by = "Cell.Type",
                     label = T,
                     label.size = 5) + xlab("UMAP 1") + ylab("UMAP 2")
umap_cell
ggsave(paste0(output_full, "UMAP_CellType.PDF"), 
       plot = umap_cell,
       width = 10,
       height = 7.5)

umap_cell_broad <- DimPlot(so, reduction = "umap",
                     group.by = "Broad.Cell.Type",
                     label = T,
                     label.size = 5) + xlab("UMAP 1") + ylab("UMAP 2")
umap_cell_broad
ggsave(paste0(output_full, "UMAP_BroadCellType.PDF"), 
       plot = umap_cell_broad,
       width = 10,
       height = 7.5)

################################################################################
####### -------               Cell type dotplot                  ------- #######

ggdot_temp <- DotPlot(so,
                                features = c("MAPT", "KDSR", "NKX2-2", "NR0B1", "STEAP1", "PRKCB", # EWS
                                             "MKI67",  "PCNA", # Proliferating
                                             "PDGFRA", "MFAP5", # Fib
                                             "CD14", "CSF1R", # Myeloid
                                             "CD34", "PECAM1", # Endothelial
                                             "CD96", "CD247", # General T cell
                                             "CD4", "RORA", # CD4+ T cells
                                             "CD8A", "CD8B", # CD8+ T cells
                                             "MYL9", "ACTA2", # Smooth muscle
                                             "JCHAIN", "MS4A1", # B cell
                                             "CA2" , "SPP1" # Osteoclasts
                                            
                                ),
                                group.by = "Cell.Type") + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
ggdot_temp
ggsave(paste0(output_dotplot_celltype_exploration, "Dotplot_CellType.PDF"), 
       plot = ggdot_temp)

ggdot_celltype_additional <- DotPlot(so,
                                features = c("MAPT", "KDSR", "NKX2-2", "NR0B1", "STEAP1", "PRKCB", # EWS
                                             "MKI67",  "PCNA", # Proliferating
                                             "PDGFRA", "MFAP5", # Fib
                                             "CD14", "CSF1R", # Myeloid
                                             "CD34", "PECAM1", # Endothelial
                                             "CD96", "CD247", # General T cell
                                             "CD4", "RORA", # CD4+ T cells
                                             "CD8A", "CD8B", # CD8+ T cells
                                             "MYL9", "ACTA2", # Smooth muscle
                                             "JCHAIN", "MS4A1", # B cell
                                             "CA2" , "SPP1", # Osteoclasts
                                             
                                             # Myofibroblasts
                                             "MCAM", "RGS5", "MMP11", "CTHRC1", "COL1A1", "COL3A1",
                                             # Myofibroblasts I added:
                                             "DES", "PDPN", "CD248", "CDH11", "TPM1", "CFL1",
                                             
                                             # Additional genes to look at
                                            # "PSTN", "ACTA2"
                                            "TAGLN", "CALD1",
                                            "VIM", "SMTN"  #, "DES"
                                             
                                ),
                                group.by = "Cell.Type") + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
ggdot_celltype_additional
ggsave(paste0(output_dotplot_celltype_exploration, "Dotplot_CellType_additional_genes.PDF"), 
       plot = ggdot_celltype_additional)

dge <- FindMarkers(so,
                   ident.1 = "Smooth Muscle",
                   ident.2 = "Fibroblast",
                   group.by = "Cell.Type",
                   logfc.threshold = 0.1)
dge_sig <- dge[dge$p_val_adj < 0.25,]


################################################################################
####### -------           Final dot plot                         ------- #######

# Cell type level
ggdot_dotplot_final <- DotPlot(so,
                      features = c("MAPT", "KDSR", "NKX2-2", "NR0B1", "STEAP1", "PRKCB", # EWS
                                   "MKI67",  "PCNA", # Proliferating
                                   "FAP", "COL1A1", # General Fib
                                   "CFD", "CXCL12", # iCAF
                                   "COMP", "CTHRC1", # mCAF
                                   "CD14", "CSF1R", # Myeloid
                                   "CA2" , "SPP1", # Osteoclasts
                                   "CD34", "PECAM1", # Endothelial
                                   "CD3D", "CD3E", # General T cell
                                   "JCHAIN", "MS4A1" # B cell
                      ),
                      group.by = "Cell.Type",
                      cluster.idents = F) + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
ggdot_dotplot_final
ggsave(paste0(output_full, "Dotplot_CellType.PDF"), 
       plot = ggdot_dotplot_final)
ggsave(paste0(output_dotplot_celltype, "Dotplot_CellType_FINAL.PDF"), 
       plot = ggdot_dotplot_final)


# Cluster level
ggdot_dotplot_final_cluster <- DotPlot(so,
                               features = c("MAPT", "KDSR", "NKX2-2", "NR0B1", "STEAP1", "PRKCB", # EWS
                                            "MKI67",  "PCNA", # Proliferating
                                            "FAP", "COL1A1", # General Fib
                                            "CFD", "CXCL12", # iCAF
                                            "COMP", "CTHRC1", # mCAF
                                            "CD14", "CSF1R", # Myeloid
                                            "CA2" , "SPP1", # Osteoclasts
                                            "CD34", "PECAM1", # Endothelial
                                            "CD3D", "CD3E", # General T cell
                                            "JCHAIN", "MS4A1" # B cell
                               ),
                               group.by = "seurat_clusters",
                               cluster.idents = T) + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
ggdot_dotplot_final_cluster
ggsave(paste0(output_full, "Dotplot_Cluster.PDF"), 
       plot = ggdot_dotplot_final_cluster)
ggsave(paste0(output_dotplot_celltype, "Dotplot_Cluster_FINAL.PDF"), 
       plot = ggdot_dotplot_final_cluster)

################################################################################
####### -------               Cell type featureplot              ------- #######

# https://github.com/satijalab/seurat/issues/1080
# List of featureplots
featureplot_celltype_final <- FeaturePlot(so,
                                          features = c("MAPT", "KDSR", "NKX2-2", "NR0B1", "STEAP1", "PRKCB", # EWS
                                                       "MKI67",  "PCNA", # Proliferating
                                                       "FAP", "COL1A1", # General Fib
                                                       "CFD", "CXCL12", # iCAF
                                                       "COMP", "CTHRC1", # mCAF
                                                       "CD14", "CSF1R", # Myeloid
                                                       "CA2" , "SPP1", # Osteoclasts
                                                       "CD34", "PECAM1", # Endothelial
                                                       "CD3D", "CD3E", # General T cell
                                                       "JCHAIN", "MS4A1" # B cell
                                          ),
                                          combine = F)
for(i in 1:length(featureplot_celltype_final)) {
  featureplot_celltype_final[[i]] <- featureplot_celltype_final[[i]] + 
    NoLegend() + NoAxes() + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
}
featureplot_celltype_final_grid <- cowplot::plot_grid(plotlist = featureplot_celltype_final)
ggsave(paste0(output_full, "UMAP_celltype_marker_grid_UPDATE.PDF"), 
       plot = print(featureplot_celltype_final_grid),
       height = 8,
       width = 10)


################################################################################
####### -------     Barplot of cell type by sample               ------- #######

# Compare distribution of new clusters and samples
prop_sample_celltype <- data.frame(prop.table(table(so@meta.data[,"Cell.Type"], 
                                              so@meta.data[,"Sample"]), margin = 2))
gg_stackbar_sample <- ggplot(prop_sample_celltype, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="stack", stat="identity") +ylab("Cell type proportion") + 
  xlab("") + labs(fill="Cell type") +
  theme_classic() + theme(axis.text.x=element_text(angle=90,
                      hjust=1, vjust = 0.5, size = 12))
gg_stackbar_sample
ggsave(paste0(output_full, "/", "Celltype_stacked_barplot.pdf"),
       gg_stackbar_sample)

# Same but styled by proportion heatmap
source("/Volumes/Partition 2/Reproducible and test scripts/scRNAseq/proportion_heatmap/prop_hm.R")
dir.create(paste0(output_full, "proportions_celltype_by_sample/"))
prop_hm(so = so,
        group1 = "Sample",
        group2 = "Cell.Type",
        order_by = T,
        outdir = paste0(output_full, "proportions_celltype_by_sample/"),
        breaks = seq(0, 1, by = 0.01), # Sets max color to 1
        height = 4,
        width = 8
)

# Should we also show proportions of non-EWS cells only?
meta_normal <- so@meta.data
meta_normal <- meta_normal[!meta_normal$Cell.Type %in% tumor_cells,]
meta_normal$Cell.Type <- droplevels(meta_normal$Cell.Type)
prop_sample_celltype_normal <- data.frame(prop.table(table(meta_normal[,"Cell.Type"], 
                                                           meta_normal[,"Sample"]), margin = 2))
gg_stackbar_sample_normal <- ggplot(prop_sample_celltype_normal, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="stack", stat="identity") +ylab("Normal cell type proportion") + 
  xlab("") + labs(fill="Cell type") +
  theme_classic() + theme(axis.text.x=element_text(angle=90,
                                                   hjust=1, vjust = 0.5, size = 12))
gg_stackbar_sample_normal
ggsave(paste0(output_full, "/", "Celltype_NormalOnly_stacked_barplot.pdf"),
       gg_stackbar_sample_normal)

# Same but styled by proportion heatmap
dir.create(paste0(output_full, "proportions_NORMAL_celltype_by_sample/"))
so_norm1 <- subset(so, Cell.Type %in% tumor_cells,
                   invert = T)
table(so_norm1$Cell.Type)
so_norm1$Cell.Type <- as.character(so_norm1$Cell.Type)
prop_hm(so = so_norm1,
        group1 = "Sample",
        group2 = "Cell.Type",
        order_by = T,
        outdir = paste0(output_full, "proportions_NORMAL_celltype_by_sample/"),
        breaks = seq(0, 1, by = 0.01), # Sets max color to 1
        height = 4,
        width = 8
)

################################################################################
####### -------               EWS markers by patient             ------- #######

so_ews <- subset(so,
                 Broad.Cell.Type %in% tumor_cells)

ggdot_ews_sample <- DotPlot(so_ews,
                                features = c("MAPT", "KDSR", "NKX2-2", "NR0B1", 
                                             "STEAP1", "PRKCB" # EWS
                                ),
                                group.by = "Sample") + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5)) +
  ggtitle("EWS cells by sample")
ggdot_ews_sample
ggsave(paste0(output_full, "Dotplot_EWSmarkers_by_sample.PDF"), 
       plot = ggdot_ews_sample)

################################################################################
####### -------       DGE and targeted ORA by cell type          ------- #######

####### -------       Output
output_ORA <- paste0(output_full, "/ORA_cell_type/")
dir.create(output_ORA)

####### -------       subset
so_sub <- so
Idents(so_sub) <- so_sub$Cell.Type
table(Idents(so_sub))

####### -------       ORA-based method
so.markers_orig <- FindAllMarkers(so_sub, 
                                  only.pos = TRUE, 
                                  min.pct = 0.1, 
                                  logfc.threshold = 0.25)
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
# Get gene sets in clusterprofiler format
gene_sets2 = panglao_gene_set_plus[,c("Gene_Set_Shortname", "gene_symbol")]
specific_output2 <- "PanglaoDB_plus"
length(unique(gene_sets2$Gene_Set_Shortname))

# EWS gene sets
ews_genesets <- read.table(file = "helper_data/EWS_gene_sets.txt",
                           sep = "\t", header = T, stringsAsFactors = F)
specific_output3 <- "ORA_EWS_genesets"
gene_sets3 <- all[all$gs_name %in% ews_genesets$GeneSetName,] %>% 
  dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
length(unique(gene_sets3$gs_name))

# Cycle through each cluster and perform ORA
specific_outputs <- c(specific_output1,
                      specific_output2,
                      specific_output3)
gene_sets_list <- list(gene_sets1,
                       gene_sets2,
                       gene_sets3) 
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
####### -------       DGE and targeted ORA by cell type          ------- #######

####### -------       Output
output_ORA <- paste0(output_full, "/ORA_broad_cell_type/")
dir.create(output_ORA)

####### -------       subset
so_sub <- so
Idents(so_sub) <- so_sub$Broad.Cell.Type

####### -------       ORA-based method
so.markers_orig <- FindAllMarkers(so_sub, 
                                  only.pos = TRUE, 
                                  min.pct = 0.1, 
                                  logfc.threshold = 0.25)
so.markers <- so.markers_orig[so.markers_orig$p_val_adj < 0.05,] # Only keep sig genes
write.table(so.markers,
            file=paste0(output_ORA,
                        "Broad_cell_type_FindAllMarkers_sigOnly.txt"),
            sep = "\t",
            row.names = F)

so.markers <- read.table(file=paste0(output_ORA,
                                     "Broad_cell_type_FindAllMarkers_sigOnly.txt"),
                         sep = "\t", stringsAsFactors = F, header = T)

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

# EWS gene sets
ews_genesets <- read.table(file = "helper_data/EWS_gene_sets.txt",
                           sep = "\t", header = T, stringsAsFactors = F)
specific_output3 <- "ORA_EWS_genesets"
# Adding one more gene set here
gene_sets3 <- all[all$gs_name %in% c(ews_genesets$GeneSetName,
                                     "ZHANG_TARGETS_OF_EWSR1_FLI1_FUSION"),] %>% 
  dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
length(unique(gene_sets3$gs_name))

# Cycle through each cluster and perform ORA
specific_outputs <- c(specific_output1,
                      specific_output2,
                      specific_output3)
gene_sets_list <- list(gene_sets1,
                       gene_sets2,
                       gene_sets3) 
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
####### -------             Plot additional genes                ------- #######

output_full_add_markers <- paste0(output_full, "Select_genes_by_celltype_v2/")
dir.create(output_full_add_markers)

for(gene in c("TSPAN8")){
  
  for(sample in c(samples, "ALL_SAMPLES")){

    if(!(sample == "ALL_SAMPLES")){
      so_sub <- subset(so,
                       Sample == sample) 
    }else{
      so_sub <- so
    }
  
  gg_gene_umap <- FeaturePlot(so_sub, reduction = "umap",
                              label = F,
                              features = gene,
                              cols = c("blue", "purple", "red"))
  
  gg_gene_vio_celltype <- so_sub %>%
    VlnPlot(
      features = gene,
      ncol     = 1,
      pt.size  = 0.01,
      group.by = "Cell.Type") + theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  # Aggregate all plots
  gg_arr <- grid.arrange(gg_gene_umap,
                         gg_gene_vio_celltype,
                         ncol = 2, nrow = 1)
  
  ggsave(paste0(output_full_add_markers, gene, "_", sample, ".pdf"),
         gg_arr,
         width = 11,
         height = 7)
  
  }
}


################################################################################
####### -------      Sample 005 DGE at cell type level           ------- #######

# This will be similar to the comparison I made in 005 CTC

output_005_dge_celltype <- paste0(output_full, "DGE_005_celltype/")
dir.create(output_005_dge_celltype)


####### -------      Differential expression
DefaultAssay(so) <- "RNA"

so_005 <- subset(so, Sample == "005_Primary")
table(so_005$Cell.Type)
so_005 <- subset(so_005, Cell.Type == "Ewing sarcoma proliferating",
                 invert = T) # Removing this to not get artifacts

dge_005 <- FindMarkers(so_005,
                       ident.1 = "Ewing sarcoma",
                       group.by = "Cell.Type",
                       logfc.threshold = 0,
                       min.pct = 0.1)
dge_005 <- cbind.data.frame(Gene = rownames(dge_005), dge_005)
write.table(dge_005,
            file = paste0(output_005_dge_celltype, "ews_vs_else_005_Primary_DGE_all_tested_genes.txt"),
            sep = "\t", row.names = F)
dge_005 <- read.table(file = paste0(output_005_dge_celltype, "ews_vs_else_005_Primary_DGE_all_tested_genes.txt"),
                      sep = "\t", header = T, stringsAsFactors = F)


####### -------      Volcano plot
library(ggrepel)

# Genes to label
genes_to_label <- c("KDSR", "NKX2-2", "STEAP1", # Canonical Ewing markers
                    "NPW", "NPY", "SCGB3A2", "TIMP1", "FN1", "SELENOP", "SERPINF1", "IGFBP5", "PAPPA", "CNMD", # Secreted proteins
                    "TSPAN8", "FCGRT", "TSPAN13", "ITM2A", # cell surface proteins
                    "MYOM2", "TPM2", # Cell motility related
                    "FOS", "JUN", "HES1", "CITED2", "TCIM" # Other genes to highlight due to contribution to Ewing or other cancer signaling.
)

adj_val <- 0.05 

df <- dge_005

df$neg_log10_fdr <- -log10(df$p_val_adj)

gg <- ggplot(df) +
  geom_point(data = subset(df, p_val_adj >= adj_val),
             aes(x = avg_log2FC, y = neg_log10_fdr),
             size = 1, color = "grey") + 
  geom_point(data = subset(df, p_val_adj < adj_val & avg_log2FC > 0),
             aes(x = avg_log2FC, y = neg_log10_fdr),
             size = 1, color = "#de6e56") + 
  geom_point(data = subset(df, p_val_adj < adj_val & avg_log2FC < 0),
             aes(x = avg_log2FC, y = neg_log10_fdr),
             size = 1, color = "#22a7f0") + 
  geom_point(data = subset(df, Gene %in% genes_to_label),
             aes(x = avg_log2FC, y = neg_log10_fdr),
             size = 2, color = "black") + 
  theme_classic() + 
  geom_text_repel(data=subset(df, Gene %in% genes_to_label), 
                  aes(label=Gene,
                      x = avg_log2FC, y = neg_log10_fdr),
                  nudge_y      = 0,
                  vjust        = 0.5,
                  angle        = 0,
                  force = 2.5,
                  point.padding = 1,
                  max.overlaps = 20) +
  labs(title = "EwS vs other cells in 005 primary tumor",
       x = "Log2 fold-change",
       y = "-log10(fdr)") 
gg <- gg + theme(axis.text=element_text(colour="black"))
gg
ggsave(paste0(output_005_dge_celltype,
              "DGE_005_ews_vs_else_volcano_EwS_genes_labeled.pdf"),
       gg,
       width = 6, height = 6)



####### -------      Similar to above but labeling top genes in both directions

dge_005_very_sig <- dge_005[dge_005$p_val_adj < 10e-25,] # Had to lower from 
table(dge_005_very_sig$avg_log2FC > 0)
dge_005_very_sig <- dge_005_very_sig[order(dge_005_very_sig$avg_log2FC, decreasing = T),]

to_plot <- 10

genes_to_label <- c(dge_005_very_sig$Gene[1:to_plot],
                    dge_005_very_sig$Gene[(nrow(dge_005_very_sig)-to_plot+1):nrow(dge_005_very_sig)])

gg <- ggplot(df) +
  geom_point(data = subset(df, p_val_adj >= adj_val),
             aes(x = avg_log2FC, y = neg_log10_fdr),
             size = 1, color = "grey") + 
  geom_point(data = subset(df, p_val_adj < adj_val & avg_log2FC > 0),
             aes(x = avg_log2FC, y = neg_log10_fdr),
             size = 1, color = "#de6e56") + 
  geom_point(data = subset(df, p_val_adj < adj_val & avg_log2FC < 0),
             aes(x = avg_log2FC, y = neg_log10_fdr),
             size = 1, color = "#22a7f0") + 
  geom_point(data = subset(df, Gene %in% genes_to_label),
             aes(x = avg_log2FC, y = neg_log10_fdr),
             size = 2, color = "black") + 
  theme_classic() + 
  geom_text_repel(data=subset(df, Gene %in% genes_to_label), 
                  aes(label=Gene,
                      x = avg_log2FC, y = neg_log10_fdr),
                  nudge_y      = 0,
                  vjust        = 0.5,
                  angle        = 0,
                  force = 2.5,
                  point.padding = 1,
                  max.overlaps = 20) +
  labs(title = "EwS vs other cells in 005 primary tumor",
       x = "Log2 fold-change",
       y = "-log10(fdr)") 
gg <- gg + theme(axis.text=element_text(colour="black"))
gg
ggsave(paste0(output_005_dge_celltype,
              "DGE_005_ews_vs_else_volcano_TOP_genes_labeled.pdf"),
       gg,
       width = 6, height = 6)




################################################################################
####### -------             Additional genes                     ------- #######

genes <- c("KDSR", "NKX2-2", "STEAP1", # Canonical Ewing markers
           "NPW", "NPY", "SCGB3A2", "TIMP1", "FN1", "SELENOP", "SERPINF1", "IGFBP5", "PAPPA", "CNMD", # Secreted proteins
           "TSPAN8", "FCGRT", "TSPAN13", "ITM2A", # cell surface proteins
           "MYOM2", "TPM2", # Cell motility related
           "FOS", "JUN", "HES1", "CITED2", "TCIM" # Other genes to highlight due to contribution to Ewing or other cancer signaling.
)

so_sub <- ScaleData(object = so, 
                    features = unique(genes))
Idents(so_sub) <- so_sub$Cell.Type
maxcells  <- min(table(Idents(so_sub)))
hm_genes <- DoHeatmap(subset(so_sub, 
                             downsample = maxcells), 
                      features = genes,
                      angle = 90,
                      group.bar.height = 0)
ggsave(paste0(output_full, "Summary_heatmap_CTC_genes.pdf"),
       hm_genes,
       width = 6,
       height = 3)


################################################################################
####### -------          save seurat object                      ------- #######

saveRDS(so,
        "Data/Seurat_objects/AggrNoNorm_filtered_primary_celltype.RDS")


