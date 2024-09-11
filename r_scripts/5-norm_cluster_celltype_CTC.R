

################################################################################
####### -------                    Goals                         ------- #######

# Processing the CTC samples only

################################################################################
####### -------                 Load R packages                  ------- #######

library(Seurat)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(msigdbr)
library(dplyr)
library(clusterProfiler)
library(openxlsx)
library(enrichplot)
library(ggrepel)

source("functions/SeuratQC.R")

################################################################################
####### -------                Update variables                  ------- #######

setwd("/Volumes/Partition 2/Core/Hayashi_Ewing_scRNAseq_Dec2022/")
input_rds = "Data/Seurat_objects/AggrNoNorm_filtered.RDS"

samples <- c("001_CTC",
             "002_CTC",
             "005_CTC",
             "006_CTC",
             "009_CTC",
             "010_CTC",
             "012_CTC",
             "038_CTC"
             )

output <- paste0("Analysis/Annotation/CTC/")
dir.create(output, recursive = T)

output_full <- paste0(output, "Normalized/")
dir.create(output_full)

# Cluster identities - Only known after the cell type identity section (skip in the beginning)
annotation <- read.xlsx("helper_data/CTC_cell_types.xlsx")
cluster_type <- as.list(annotation$Cell.Type)
names(cluster_type) <- annotation$Cluster
cluster_broad <- as.list(annotation$Broad.Cell.Type)
names(cluster_broad) <- annotation$Cluster

################################################################################
####### -------           Read in and subset data                ------- #######

# read in Seurat object and split
seurat_obj <- readRDS(input_rds)

# Subset to sample
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
so <- RunPCA(so,
             features = VariableFeatures(so),
             verbose = F,
             ndims.print = 0)
ElbowPlot(so,
          ndims = 40)
n_dim <- 20 # Used 30 for primary samples
so <- RunUMAP(so, 
              reduction = "pca", 
              dims = 1:n_dim)

gg_sample <- DimPlot(so, 
                     group.by = "Sample",
                     reduction = "umap")
ggsave(paste0(output_full, "UMAP_sample.pdf"), gg_sample)

gg_sample_split <- DimPlot(so, 
                     group.by = "Sample",
                     reduction = "umap",
                     split.by = "Sample")
ggsave(paste0(output_full, "UMAP_sample_split.pdf"), gg_sample_split, height = 3, width = 15)

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

resolutionsCalculated <- c( 0.1,  0.3,  0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1)
so <- FindClusters(so, resolution = resolutionsCalculated) 

for (i in resolutionsCalculated){
  # set resolution
  Idents(object = so) <- paste(assay, "_snn_res." , i, sep = "") # name is picked from column name of resolutions
  
  # general umap
  allClusters <- DimPlot(so,
                         reduction = "umap",
                         label = TRUE,
                         label.size = 6) + ggtitle(paste(assay, "_snn_res." , i, sep = "")) + NoLegend() + theme(plot.title = element_text(hjust = 0.5))
  print(allClusters)
  
}

# Choose a resolution
resolutionChosen <- 0.8 

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
####### -------  Cluster-level QC (whole dataset and by sample)  ------- #######

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
####### -------   inferCNV supervised 005-only (higher bayes p)  ------- #######

so_sub_infercnv <- subset(so,
                          Sample == "005_CTC")

column_identity <- "seurat_clusters"

output_run <- paste0(output_full, "InferCNV/Supervised_005_subclusters_bayes0.5/")
dir.create(output_run,
           recursive = T)
output_run_select <- paste0(output_run, "/select_results/")
dir.create(output_run_select,
           recursive = T)

#################   Run inferCNV

#################    Create a cell annotations file
cell_ann <- cbind.data.frame(Cell = rownames(so_sub_infercnv@meta.data),
                             Cluster = so_sub_infercnv@meta.data[,column_identity])
# Need to remove and clusters with <2 cells since inferCNV won't do it automatically
cell_ann_table <- data.frame(table(cell_ann$Cluster))
clusters_keep <- as.character(cell_ann_table$Var1[cell_ann_table$Freq>1])
cell_ann <- cell_ann[cell_ann$Cluster %in% clusters_keep ,]
cell_ann_out <- paste0(output_run, "cell_ann/")
dir.create(cell_ann_out,
           recursive = T)
write.table(cell_ann,
            file = paste0(cell_ann_out, "_005.txt"),
            row.names = F, col.names = F, sep = "\t", quote = F)

# Need to remove any references that aren't in this sample
normal_cells_sub <- as.vector(unique(cell_ann$Cluster)[!(unique(cell_ann$Cluster) == "11")])

#################   Get raw counts - need to re-subset based on above annotation file

so_sub <- subset(so_sub_infercnv,
                 cells = colnames(so_sub_infercnv)[colnames(so_sub_infercnv) %in% cell_ann$Cell])
exp.rawdata <- as.matrix(so_sub_infercnv@assays$RNA@counts)

# Create infercnv object
infercnv_obj = infercnv::CreateInfercnvObject(raw_counts_matrix= exp.rawdata,
                                              annotations_file=paste0(cell_ann_out, "_005.txt"),
                                              delim="\t",
                                              gene_order_file="helper_data/infercnv_gene_position.txt",
                                              ref_group_names = normal_cells_sub)

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir= output_run,  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             HMM_type='i6', # Changed for uphyloplot2
                             denoise=T,
                             HMM=T,
                             analysis_mode = "subclusters",
                             num_threads = 6,
                             BayesMaxPNormal = 0.5,
                             output_format = "pdf" # default: "png"
)



################################################################################
####### -------       DGE in sample 005-only.                    ------- #######

####### -------       ORA-based method

so_005 <- subset(so, Sample == "005_CTC")
Idents(so_005) <- so_005$seurat_clusters

so_005.markers_orig <- FindAllMarkers(so_005, 
                                  only.pos = TRUE, 
                                  min.pct = 0.1, 
                                  logfc.threshold = 0.25)
so_005.markers <- so_005.markers_orig[so_005.markers_orig$p_val_adj < 0.05,] # Only keep sig genes
output_dge_005 <- paste0(output_full, "/DGE_005_only/")
dir.create(output_dge_005)
write.table(so_005.markers,
            file=paste0(output_dge_005,
                        "Clusters_FindAllMarkers_sigOnly_005only.txt"),
            sep = "\t",
            row.names = F)

################################################################################
####### -------       Help with cell type identity               ------- #######

####### -------       ORA-based method

Idents(so) <- so$seurat_clusters

so.markers_orig <- FindAllMarkers(so, 
                                  only.pos = TRUE, 
                                  min.pct = 0.1, 
                                  logfc.threshold = 0.25)
so.markers <- so.markers_orig[so.markers_orig$p_val_adj < 0.05,] # Only keep sig genes
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
  slice_max(n = 200, order_by = avg_log2FC) # 

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
####### -------               Cell type dotplot                  ------- #######

ggdot_cluster_final <- DotPlot(so,
                                features = c("KDSR", "NKX2-2",  "STEAP1",  # EWS
                                             "NPW", "TSPAN8", # Suspected EWS
                                             "CD1C", "FCER1A", "HLA-DQA1", "HLA-DQB1", # 13 = Dendritic
                                             "PPBP", "PF4", # Platelets
                                             "HBD",  "AHSP", # Erythrocyte
                                             "FCGR3A",  # 7 = CD16+FCN1+ Monocyte # And most of the markers in 0,6,8
                                             "CD14", "FCN1", "VCAN", "S100A8", "S100A9", # 0,6,8 = CD14+FCN1+ Monocyte
                                             "CD79A", "MS4A1", # B cell
                                             "CD3D", "CD3E", # General T cell
                                             "TCF7", "LTB", "CCR7", "LEF1", # Naive T
                                             "CLEC4C", "IL3RA", # 12 = pDC
                                             "NKG7", "PRF1", "GZMB", "GNLY" # Cytotoxic T/NK cells
                                ),
                                group.by = "seurat_clusters",
                               cluster.idents = T) + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
ggdot_cluster_final
ggsave(paste0(output_dotplot_celltype, "Dotplot_cluster.PDF"), 
       plot = ggdot_cluster_final)

################################################################################
####### -------              Plot select genes                   ------- #######

output_full_add_markers <- paste0(output_full, "Select_gene_plots/")
dir.create(output_full_add_markers)

genes <- c("TSPAN8", 
           "NPW", 
           "SCGB3A2",
           "TIMP1", 
           "CCND1", 
           "FOS", 
           "ID1", 
           "VIM", 
           "IGFBP5")

for(gene in genes){
  
  gg_gene_umap <- FeaturePlot(so, reduction = "umap",
                              label = F,
                              features = gene,
                              order = T)
  gg_gene_vio <- so %>%
    VlnPlot(
      features = gene,
      ncol     = 1,
      pt.size  = 0.01,
      group.by = "seurat_clusters" #,
      #y.max = 15,
    ) + theme(legend.position = "none")
  
  # Aggregate all plots
  gg_arr <- grid.arrange(gg_gene_umap,
                         gg_gene_vio,
                         #gg_gene_vio_celltype,
                         ncol = 2, nrow = 1)
  
  ggsave(paste0(output_full_add_markers, gene, ".pdf"),
         gg_arr,
         width = 11,
         height = 7)
  
}

so <- ScaleData(object = so, 
                features = unique(c(VariableFeatures(so) ,
                                    genes)),
                vars.to.regress = c("nCount_RNA",
                                    "percent_mito"))
maxcells  <- min(table(Idents(so)))
hm_genes <- DoHeatmap(subset(so, downsample = maxcells), 
                      features = genes,
                      angle = 0)
hm_genes
ggsave(paste0(output_full_add_markers, "Summary_heatmap.pdf"),
       hm_genes,
       width = 11,
       height = 7)
hm_genes <- DoHeatmap(subset(so, downsample = maxcells), 
                      features = genes,
                      angle = 0) + NoLegend()
hm_genes
ggsave(paste0(output_full_add_markers, "Summary_heatmap_NoLegend.pdf"),
       hm_genes,
       width = 11,
       height = 7)

ggdot_genes <- DotPlot(so,
                       features = genes,
                       group.by = "seurat_clusters") + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
ggdot_genes 
ggsave(paste0(output_full_add_markers, "Summary_dotplot.pdf"),
       ggdot_genes )

featureplot_multi <- FeaturePlot(so,
                                 features = genes,
                                 combine = F,
                                 order=T)
for(i in 1:length(featureplot_multi)) {
  featureplot_multi[[i]] <- featureplot_multi[[i]] + 
    NoLegend() + NoAxes() + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
}
featureplot_multi_grid <- cowplot::plot_grid(plotlist = featureplot_multi)
ggsave(paste0(output_full_add_markers, "Summary_featureplot_grid.pdf"), 
       plot = print(featureplot_multi_grid),
       height = 9,
       width = 9)




################################################################################
####### -------     Multi-gene summary plots - 005-only          ------- #######

output_full_add_markers <- paste0(output_full, "Select_gene_plots_005/")
dir.create(output_full_add_markers)

so_sub <- subset(so, Sample == "005_CTC")

Idents(so_sub) <- so_sub$Cell.Type

genes <- c("KDSR", "NKX2-2", "STEAP1", # Canonical Ewing markers
           "NPW", "NPY", "SCGB3A2", "TIMP1", "FN1", "SELENOP", "SERPINF1", "IGFBP5", "PAPPA", "CNMD", # Secreted proteins
           "TSPAN8", "FCGRT", "TSPAN13", "ITM2A", # cell surface proteins
           "MYOM2", "TPM2", # Cell motility related
           "FOS", "JUN", "HES1", "CITED2", "TCIM" # Other genes to highlight due to contribution to Ewing or other cancer signaling.
)

for(gene in genes){
  
  gg_gene_umap <- FeaturePlot(so_sub, reduction = "umap",
                              label = F,
                              features = gene,
                              order = T)
  gg_gene_vio <- so_sub %>%
    VlnPlot(
      features = gene,
      ncol     = 1,
      pt.size  = 0.01,
      group.by = "Cell.Type" #,
      #y.max = 15,
    ) + theme(legend.position = "none")
  
  # Aggregate all plots
  gg_arr <- grid.arrange(gg_gene_umap,
                         gg_gene_vio,
                         #gg_gene_vio_celltype,
                         ncol = 2, nrow = 1)
  
  ggsave(paste0(output_full_add_markers, gene, ".pdf"),
         gg_arr,
         width = 11,
         height = 7)
  
}

so_sub <- ScaleData(object = so_sub, 
                features = unique(genes))
maxcells <- sum(so_sub$Cell.Type == "Ewing sarcoma")
hm_genes <- DoHeatmap(subset(so_sub, 
                      downsample = maxcells), 
                      features = genes,
                      angle = 90,
                      group.bar.height = 0)
ggsave(paste0(output_full_add_markers, "Summary_heatmap.pdf"),
       hm_genes,
       width = 8,
       height = 16)
hm_genes <- DoHeatmap(subset(so_sub, downsample = maxcells), 
                      features = genes,
                      angle = 90,
                      group.bar.height = 0) + NoLegend()
ggsave(paste0(output_full_add_markers, "Summary_heatmap_NoLegend.pdf"),
       hm_genes,
       width = 8,
       height = 16)

ggdot_genes <- DotPlot(so,
                       features = genes,
                       group.by = "Cell.Type") + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
ggdot_genes 
ggsave(paste0(output_full_add_markers, "Summary_dotplot.pdf"),
       ggdot_genes )

featureplot_multi <- FeaturePlot(so_sub,
                                 features = genes,
                                 combine = F,
                                 order=T)
for(i in 1:length(featureplot_multi)) {
  featureplot_multi[[i]] <- featureplot_multi[[i]] + 
    NoLegend() + NoAxes() + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
}
featureplot_multi_grid <- cowplot::plot_grid(plotlist = featureplot_multi)
ggsave(paste0(output_full_add_markers, "Summary_featureplot_grid.pdf"), 
       plot = print(featureplot_multi_grid),
       height = 9,
       width = 9)

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

so$Cell.Type <- factor(so$Cell.Type,
                       levels = c("Ewing sarcoma",
                                  "CD16+FCN1+ Monocyte",
                                  "CD14+FCN1+ Monocyte",
                                  "Plasmacytoid dendritic",
                                  "Dendritic",
                                  "Naive T",
                                  "Cytotoxic T/NK",
                                  "Platelet",
                                  "B",
                                  "Erythrocyte"))

################################################################################
####### -------            Cluster-level marker dotplot         ------- #######

so$Cluster.Cell.Type <- paste0(so$seurat_clusters, ".", so$Cell.Type)

ggdot_cluster_final <- DotPlot(so,
                               features = c("CD14", "FCN1", "VCAN", "S100A8", "S100A9", # 0,6,8 = CD14+FCN1+ Monocyte
                                            "CD3D", "CD3E", # General T cell
                                            "TCF7", "LTB", "CCR7", "LEF1", # Naive T
                                            "PPBP", "PF4", # Platelets
                                            "NKG7", "PRF1", "GZMB", "GNLY", # Cytotoxic T/NK cells
                                            
                                            "HBD",  "AHSP", # Erythrocyte
                                            "CLEC4C", "IL3RA", # 12 = pDC
                                            "CD79A", "MS4A1", # B cell
                                            "FCGR3A",  # 7 = CD16+FCN1+ Monocyte # And most of the markers in 0,6,8
                                            "KDSR", "NKX2-2",  "STEAP1",  # EWS
                                            "NPW", "TSPAN8", # Suspected EWS
                                            "CD1C", "FCER1A", "HLA-DQA1", "HLA-DQB1" # 13 = Dendritic
                               ),
                               group.by = "Cluster.Cell.Type",
                               cluster.idents = T) + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
ggdot_cluster_final
ggsave(paste0(output_dotplot_celltype, "Dotplot_cluster_cell_type.PDF"), 
       plot = ggdot_cluster_final)

################################################################################
####### -------                 Cell type UMAP                   ------- #######

umap_cell <- DimPlot(so, reduction = "umap",
                     group.by = "Cell.Type",
                     label = T,
                     label.size = 5) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(paste0(output_full, "UMAP_CellType.PDF"), 
       plot = umap_cell,
       width = 10,
       height = 7.5)

umap_cell_broad <- DimPlot(so, reduction = "umap",
                     group.by = "Broad.Cell.Type",
                     label = T,
                     label.size = 5) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave(paste0(output_full, "UMAP_BroadCellType.PDF"), 
       plot = umap_cell_broad,
       width = 10,
       height = 7.5)

################################################################################
####### -------               Cell type dotplot                  ------- #######

ggdot_celltype_final <- DotPlot(so,
                          features = c("KDSR", "NKX2-2",  "STEAP1",  # EWS
                                       "NPW", "TSPAN8", # Suspected EWS
                                       "FCGR3A",  # 7 = CD16+FCN1+ Monocyte # And most of the markers in 0,6,8
                                       "CD14", "FCN1", "VCAN", "S100A8", "S100A9", # 0,6,8 = CD14+FCN1+ Monocyte
                                       "CLEC4C", "IL3RA", # 12 = pDC
                                       "CD1C", "FCER1A", "HLA-DQA1", "HLA-DQB1", # 13 = Dendritic
                                       "CD3D", "CD3E", # General T cell
                                       "TCF7", "LTB", "CCR7", "LEF1", # Naive T
                                       "NKG7", "PRF1", "GZMB", "GNLY", # Cytotoxic T/NK cells
                                       "PPBP", "PF4", # Platelets
                                       "CD79A", "MS4A1", # B cell
                                       "HBD",  "AHSP" # Erythrocyte
                          ),
                          group.by = "Cell.Type") + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
ggdot_celltype_final
ggsave(paste0(output_full, "Dotplot_CellType.PDF"), 
       plot = ggdot_celltype_final)
ggsave(paste0(output_dotplot_celltype, "Dotplot_CellType.PDF"), 
       plot = ggdot_celltype_final)

################################################################################
####### -------               Cell type featureplot              ------- #######

# https://github.com/satijalab/seurat/issues/1080
# List of featureplots
featureplot_celltype_final <- FeaturePlot(so,
                                          features = c("KDSR", "NKX2-2",  "STEAP1",  # EWS
                                                       "NPW", "TSPAN8", # Suspected EWS
                                                       "FCGR3A",  # 7 = CD16+FCN1+ Monocyte # And most of the markers in 0,6,8
                                                       "CD14", "FCN1", "VCAN", "S100A8", "S100A9", # 0,6,8 = CD14+FCN1+ Monocyte
                                                       "CLEC4C", "IL3RA", # 12 = pDC
                                                       "CD1C", "FCER1A", "HLA-DQA1", "HLA-DQB1", # 13 = Dendritic
                                                       "CD3D", "CD3E", # General T cell
                                                       "TCF7", "LTB", "CCR7", "LEF1", # Naive T
                                                       "NKG7", "PRF1", "GZMB", "GNLY", # Cytotoxic T/NK cells
                                                       "PPBP", "PF4", # Platelets
                                                       "CD79A", "MS4A1", # B cell
                                                       "HBD",  "AHSP" # Erythrocyte
                                          ),
                                          order =T,
                                          combine = F)
for(i in 1:length(featureplot_celltype_final)) {
  featureplot_celltype_final[[i]] <- featureplot_celltype_final[[i]] + 
    NoLegend() + NoAxes() + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
}
featureplot_celltype_final_grid <- cowplot::plot_grid(plotlist = featureplot_celltype_final)

ggsave(paste0(output_full, "UMAP_celltype_marker_grid.PDF"), 
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
ggsave(paste0(output_dotplot_celltype, "/", "Celltype_stacked_barplot.pdf"),
       gg_stackbar_sample)

# Same but styled by proportion heatmap
source("functions/proportion_heatmap/prop_hm.R")
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
ggsave(paste0(output_dotplot_celltype, "/", "Celltype_NormalOnly_stacked_barplot.pdf"),
       gg_stackbar_sample_normal)

################################################################################
####### -------               EWS markers by patient             ------- #######

so_ews <- subset(so, Broad.Cell.Type %in% tumor_cells)

ggdot_ews_sample <- DotPlot(so_ews,
                                features = c("MAPT", "KDSR", "NKX2-2", "NR0B1", 
                                             "STEAP1", "PRKCB" # EWS
                                ),
                                group.by = "Sample") + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5)) +
  ggtitle("EWS cells by sample")
ggsave(paste0(output_full, "Dotplot_EWSmarkers_by_sample.PDF"), 
       plot = ggdot_ews_sample)

################################################################################
####### -------              Plot select genes                   ------- #######

output_full_add_markers <- paste0(output_full, "Select_gene_plots_by_cell_type/")
dir.create(output_full_add_markers)

genes <- c("TSPAN8", 
           "NKX2-2",
           "NPW", 
           "SCGB3A2",
           "TIMP1", 
           "CCND1", 
           "FOS", 
           "ID1", 
           "VIM", 
           "IGFBP5")

for(gene in genes){
  
  gg_gene_umap <- FeaturePlot(so, reduction = "umap",
                              label = F,
                              features = gene,
                              order = T,
                              cols = c("blue", "purple", "red"))
  gg_gene_vio <- so %>%
    VlnPlot(
      features = gene,
      ncol     = 1,
      pt.size  = 0.01,
      group.by = "Cell.Type"
    ) + theme(legend.position = "none")
  
  # Aggregate all plots
  gg_arr <- grid.arrange(gg_gene_umap,
                         gg_gene_vio,
                         ncol = 2, nrow = 1)
  
  ggsave(paste0(output_full_add_markers, gene, ".pdf"),
         gg_arr,
         width = 11,
         height = 7)
  
}

so <- ScaleData(object = so, 
                features = unique(c(VariableFeatures(so) ,
                                    genes)),
                vars.to.regress = c("nCount_RNA",
                                    "percent_mito"))
maxcells  <- min(table(so$Cell.Type))
so_ident <- so
Idents(so_ident) <- so_ident$Cell.Type
hm_genes <- DoHeatmap(subset(so_ident, downsample = maxcells), 
                      features = genes,
                      angle = 45,
                      group.by = "Cell.Type")
ggsave(paste0(output_full_add_markers, "Summary_heatmap.pdf"),
       hm_genes,
       width = 11,
       height = 7)
hm_genes <- DoHeatmap(subset(so_ident, downsample = maxcells), 
                      features = genes,
                      angle = 45) + NoLegend()
ggsave(paste0(output_full_add_markers, "Summary_heatmap_NoLegend.pdf"),
       hm_genes,
       width = 11,
       height = 7)

ggdot_genes <- DotPlot(so_ident,
                       features = genes) + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
ggdot_genes 
ggsave(paste0(output_full_add_markers, "Summary_dotplot.pdf"),
       ggdot_genes )

featureplot_multi <- FeaturePlot(so_ident,
                                 features = genes,
                                 combine = F,
                                 order=T)
for(i in 1:length(featureplot_multi)) {
  featureplot_multi[[i]] <- featureplot_multi[[i]] + 
    NoLegend() + NoAxes() + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
}
featureplot_multi_grid <- cowplot::plot_grid(plotlist = featureplot_multi)
ggsave(paste0(output_full_add_markers, "Summary_featureplot_grid.pdf"), 
       plot = print(featureplot_multi_grid),
       height = 9,
       width = 9)

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

# Genee sets
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
gene_sets3 <- all[all$gs_name %in% ews_genesets$GeneSetName,] %>% 
  dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

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
                              gsub("/", "_", cluster), ".txt"),
                row.names = F, sep = "\t", quote = F)
  }
}


################################################################################
####### -------      Sample 005 DGE at cell type level           ------- #######

output_005_dge_celltype <- paste0(output_full, "DGE_005_celltype/")
dir.create(output_005_dge_celltype)

####### -------      Differential expression
DefaultAssay(so) <- "RNA"
so_005 <- subset(so, Sample == "005_CTC")

dge_005 <- FindMarkers(so_005,
                       ident.1 = "Ewing sarcoma",
                       group.by = "Cell.Type",
                       logfc.threshold = 0,
                       min.pct = 0.1)
dge_005 <- cbind.data.frame(Gene = rownames(dge_005), dge_005)
write.table(dge_005,
            file = paste0(output_005_dge_celltype, "ews_vs_else_005_CTC_DGE_all_tested_genes.txt"),
            sep = "\t", row.names = F)
dge_005 <- read.table(file = paste0(output_005_dge_celltype, "ews_vs_else_005_CTC_DGE_all_tested_genes.txt"),
                      sep = "\t", header = T, stringsAsFactors = F)

####### -------      Volcano plot

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
    labs(title = "EwS vs other cells in 005 CTC",
         x = "Log2 fold-change",
         y = "-log10(fdr)") 
gg <- gg + theme(axis.text=element_text(colour="black"))
ggsave(paste0(output_005_dge_celltype,
              "DGE_005_ews_vs_else_volcano_EwS_genes_labeled.pdf"),
         gg,
         width = 6, height = 6)

####### -------      Similar to above but labeling top genes in both directions

dge_005_very_sig <- dge_005[dge_005$p_val_adj < 10e-30,]
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
  labs(title = "EwS vs other cells in 005 CTC",
       x = "Log2 fold-change",
       y = "-log10(fdr)") 
gg <- gg + theme(axis.text=element_text(colour="black"))
gg
ggsave(paste0(output_005_dge_celltype,
              "DGE_005_ews_vs_else_volcano_TOP_genes_labeled.pdf"),
       gg,
       width = 6, height = 6)



####### -------      Similar to above but proposed immunotherapeutic EWS targets

# surfaceome

# Gather genes
# https://aacrjournals.org/clincancerres/article/30/5/1022/734306/Surface-and-Global-Proteome-Analyses-Identify
# The prioritization generated 39 candidates in Group 1 (Z-scores ≥ 1) and 56 candidates in Group 2 (Z-scores > 0; Fig. 3C and D; Supplementary Data S3). 
it_df <- read.xlsx("Data/helper_data/mooney_surfaceome/immune_score_ranking_modified_from_supp_data_s3.xlsx")
group1 <- it_df$symbol[it_df$Group == "Zscore>1"]
group2 <- it_df$symbol[it_df$Group == "Zscore>0"]

# Label DGE table
dge_005_it <- dge_005
dge_005_it$group1 <- dge_005_it$Gene %in% group1
dge_005_it$group2 <- dge_005_it$Gene %in% group2
dge_005_it_sub <- dge_005_it[dge_005_it$group1 == "TRUE" | dge_005_it$group2 == "TRUE",]


####### -------      GSEA of EWS-associated gene sets in 005 CTC vs else

ews_genesets <- read.table(file = "/Volumes/Partition 2/Core/Hayashi_Ewing_scRNAseq_Dec2022/Data/helper_data/EWS_gene_sets.txt",
                           sep = "\t", header = T, stringsAsFactors = F)
gene_sets3 <- all[all$gs_name %in% ews_genesets$GeneSetName,] %>% 
  dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
output_gsea <- paste0(output_005_dge_celltype, "/gsea_ews_genesets/")
dir.create(output_gsea, recursive = T)
ranked <- as.numeric(dge_005$avg_log2FC)
names(ranked) <- dge_005$Gene
ranked <- ranked[order(ranked,
                       decreasing = T)]
edo2 <- GSEA(geneList = ranked,
             TERM2GENE = gene_sets3,
             minGSSize = 1,
             pvalueCutoff = 1,
             maxGSSize = 5000)
edo2_df <- data.frame(edo2)
write.table(edo2_df,
            file = paste0(output_gsea, "GSEA_table.txt"),
            sep = "\t", row.names = F, quote = F)
gg <- gseaplot2(edo2, geneSetID = 1:4, pvalue_table = T)
ggsave(paste0(output_gsea, "enrichment_all_ews_genesets.pdf"), print(gg))


################################################################################
####### ---- Merging primary and ctc to compare 005 ctc vs primary  ---- #######

output_005_comp <- paste0(output_full, "DGE_CTC005_vs_Primary005/")
dir.create(output_005_comp)

so_p <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_primary_celltype.RDS")
so_ctc <- so
so_ctc$cell <- so_ctc$seurat_clusters
so_combined <- merge(so_ctc, so_p)
so_combined$Cell.Type2 <- gsub(" proliferating", "", so_combined$Cell.Type)

# Subset to 005 only and only EWS
table(so_combined$Sample)
so_combined <- subset(so_combined,
                      Sample %in% c("005_Primary", "005_CTC"))
table(so_combined$Cell.Type)
so_combined <- subset(so_combined,
                      Cell.Type %in% c("Ewing sarcoma", "Ewing sarcoma proliferating"))

so_combined <- NormalizeData(so_combined, 
                             normalization.method = "LogNormalize")
Idents(so_combined) <- so_combined$Sample

dge <- FindMarkers(so_combined,
                   ident.1 = "005_CTC",
                   ident.2 = "005_Primary",
                   logfc.threshold = 0)
dge <- cbind.data.frame(Gene = rownames(dge),
                        dge,
                        diff.pct = dge$pct.1 - dge$pct.2)

write.table(dge,
            file = paste0(output_005_comp, "DGE_005_CTC_vs_primary_in_all_EwS_cells_all_genes_tested.txt"),
            sep = "\t", row.names = F, quote = F)
dge <- read.table(file = paste0(output_005_comp, "DGE_005_CTC_vs_primary_in_all_EwS_cells_all_genes_tested.txt"),
                      sep = "\t", header = T, stringsAsFactors = F)


####### -------      Volcano plot

adj_val <- 0.05 
df <- dge
df$neg_log10_fdr <- -log10(df$p_val_adj)

# Labeling top genes in both directions
dge_005_very_sig <- dge[dge$p_val_adj < 10e-20,]
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
  labs(title = "005 CTC vs Primary in EWS cells",
       x = "Log2 fold-change",
       y = "-log10(fdr)") 
gg <- gg + theme(axis.text=element_text(colour="black"))
gg
ggsave(paste0(output_005_comp,
              "DGE_005_CTC_vs_primary_in_EWS_volcano_TOP_genes_labeled.pdf"),
       gg,
       width = 6, height = 6)

####### -------      fGSEA using merged gene programs
output_gsea <- paste0(output_005_comp, "/gsea_merged_GEPs/")
dir.create(output_gsea, recursive = T)
merged_geps <- read.table(file = paste0("Analysis/Annotation/Primary/integrate_P_RPCA_CPM/cNMF_v3/merged_programs/Merged_gene_programs.txt"),
                          sep = "\t", header = T)
gene_sets_geps <- cbind.data.frame(gs_name = merged_geps$NMF,
                                   gene_symbol = merged_geps$Top_genes)
ranked <- as.numeric(dge$avg_log2FC)
names(ranked) <- dge$Gene
ranked <- ranked[order(ranked,
                       decreasing = T)]
edo2 <- GSEA(geneList = ranked,
             TERM2GENE = gene_sets_geps,
             minGSSize = 1,
             pvalueCutoff = 1,
             maxGSSize = 5000)
edo2_df <- data.frame(edo2)
write.table(edo2_df,
            file = paste0(output_gsea, "GSEA_table.txt"),
            sep = "\t", row.names = F, quote = F)
gg <- gseaplot2(edo2, geneSetID = 1:6, pvalue_table = T)
##-   Single enrichment plot
gg2 <- gseaplot2(edo2, geneSetID = 2, 
          title = edo2$Description[2],
          pvalue_table = T)
ggsave(paste0(output_gsea, "enrichment_all_geps.pdf"),
       print(gg))
ggsave(paste0(output_gsea, "enrichment_ews_gep.pdf"),
       print(gg2))

################################################################################
####### -------          save seurat object                      ------- #######

saveRDS(so,
        "Data/Seurat_objects/AggrNoNorm_filtered_celltype_CTC.RDS")
