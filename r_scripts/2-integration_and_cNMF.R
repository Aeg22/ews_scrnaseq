
################################################################################
####### -------                   Goals                          ------- #######

# Integrating EWS cells in primary with RPCA CPM

# cNMF

################################################################################
####### -------        Read in data (to be integrated)           ------- #######

library(dplyr) 
library(Seurat)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(msigdbr)
library(plyr)
library(clusterProfiler)
library(openxlsx)
library(pheatmap)
library(ggpubr)
library(NMF)
library(dorothea)
library(tibble)
library(viper)
library(reshape2)

# Source custom functions
source("/functions/FindNeighbors_FindClusters_PrintUMAP.R")

setwd("/Volumes/Partition 2/Core/Hayashi_Ewing_scRNAseq_Dec2022/")
sobj <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_primary_celltype.RDS")

output <- "Analysis/Annotation/Primary/integrate_P_RPCA_CPM/"
dir.create(output, recursive = T)

sobj_meta <- sobj@meta.data

################################################################################
####### -------      Subset data        ------- #######

samples <- c("051_Primary",
             "061_Primary",
             "066_Primary",
             "005_Primary",
             "010_Primary",
             "038_Primary",
             "001_Primary")
ews_clusters <- c("Ewing sarcoma", "Ewing sarcoma proliferating")
sobj_sub <- subset(sobj, Cell.Type %in% ews_clusters)

################################################################################
####### -------  Integrate USING RPCA while normalizing with CPM ------- #######

k.anchor <- 10 
int_dim <- 20
k.weight <- 25 

split_group <- "Sample"

orig <- Sys.time()
ifnb.list <- SplitObject(sobj_sub, split.by = split_group)

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE,
              npcs = 40) # Had to lower for the sample with few cells
})

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, 
                                         anchor.features = features,
                                         reduction = "rpca",
                                         dims = 1:int_dim, # Default 30
                                         k.anchor = k.anchor # Default is 10
)
immune.combined <- IntegrateData(anchorset = immune.anchors, 
                                         k.weight=k.weight # Normally set to 100. Was 25
)
Sys.time()-orig 

output_int_rpca <- paste0(output,
                          "integrated",
                          "/")
dir.create(output_int_rpca)

DefaultAssay(immune.combined) <- "integrated"

immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined,
                          verbose = FALSE)
ElbowPlot(immune.combined, ndims = 40)
immune.combined <- RunUMAP(immune.combined, 
                           reduction = "pca", 
                           dims = 1:25)

umap_CellType <- DimPlot(immune.combined, reduction = "umap",
                          group.by = "Cell.Type")
umap_sample <- DimPlot(immune.combined, reduction = "umap",
                       group.by = "Sample")
umap_cnasub <- DimPlot(immune.combined, reduction = "umap",
                       group.by = "CNAsubclusters")
umap_cnasub_split <- subset(immune.combined,
                            CNAsubclusters == "none",
                            invert=T) %>% DimPlot( reduction = "umap",
                       group.by = "CNAsubclusters",
                       split.by = "Sample",
                       order = "005_Primary_subcluster.1.2")

ggsave(paste0(output_int_rpca, "UMAP_CellType.pdf"), umap_CellType)
ggsave(paste0(output_int_rpca, "UMAP_sample.pdf"), umap_sample)
ggsave(paste0(output_int_rpca, "UMAP_CNA_subclusters.pdf"), umap_cnasub)
ggsave(paste0(output_int_rpca, "UMAP_CNA_subclusters_split.pdf"), umap_cnasub_split)

################################################################################
####### -------                   Clustering                     ------- #######

# Cluster, choose resolution
immune.combined <- FindNeighbors_FindClusters_PrintUMAP(sobj = immune.combined,
                                                 variable = "Sample",
                                                 neighbors_dim = 1:15,
                                                 resolutionsCalculated = c(0.1,
                                                0.2, 0.3, 0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1))
resolutionChosen <- 0.8
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
ggsave(paste0(output_int_rpca, "/", "Sample_stacked_barplot.pdf"),
       gg_stackbar_sample)

saveRDS(immune.combined,
        "Data/Seurat_objects/AggrNoNorm_filtered_primary_celltype_CNAsubs_integrateP_clustered.RDS")

################################################################################
####### -----                        cNMF                         ------ #######

# https://github.com/dylkot/cNMF

####### -----            Read in seurat object and define primary output directory

so <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_primary_celltype_CNAsubs_integrateP_clustered_pathway.RDS")
DefaultAssay(so) <- "RNA"

output_nmf <- paste0(output, "cNMF_v3/")
dir.create(output_nmf)

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


####### -----            Write counts files for each sample

# Will write raw counts per sample as well as raw counts per sample with only
# like genes across samples. SPECIFY WHICH TO USE WITH input_counts_folder

output_nmf_counts <- paste0(output_nmf, "/counts/")
dir.create(output_nmf_counts)
output_nmf_counts_like_genes <- paste0(output_nmf, "/counts_like_genes/")
dir.create(output_nmf_counts_like_genes)

genes_across_sample <- c()
for(sample in unique(so$Sample)){
  
  # Subset
  print(sample)
  so_sub <- subset(so, Sample == sample)
  so_sub_atkins <- CreateSeuratObject(counts = so_sub@assays$RNA@counts, 
                                      meta.data = so_sub@meta.data)
  
  # Process and write raw counts
  # cNMF is expecting cells as rows and genes as columns
  counts_matrix <- data.frame(t(as.matrix(GetAssayData(object = so_sub, slot = "counts"))),
                              check.names=F)
  counts_matrix <- counts_matrix[,colSums(counts_matrix) > 0] # Removing genes without counts
  counts_matrix <- cbind.data.frame(Cell = rownames(counts_matrix),
                                    counts_matrix)
  colnames(counts_matrix)[1] <- ""
  genes_across_sample <- c(genes_across_sample,
                           colnames(counts_matrix)[2:ncol(counts_matrix)])
  write.table(counts_matrix, 
              paste0(output_nmf_counts, sample, "_raw_counts.txt"), 
              sep = '\t', 
              row.names = F, 
              col.names = T, 
              quote = F)
}

# Also writing counts where only like genes across all samples are included
genes_across_sample <- data.frame(unclass(table(genes_across_sample)))
genes_to_keep <- rownames(genes_across_sample)[genes_across_sample[,1] == length(unique(so$Sample))]
for(sample in unique(so$Sample)){
  
  # Subset
  print(sample)
  so_sub <- subset(so, Sample == sample)
  so_sub_atkins <- CreateSeuratObject(counts = so_sub@assays$RNA@counts, 
                                      meta.data = so_sub@meta.data)
  
  # Process and write raw counts
  # cNMF is expecting cells as rows and genes as columns
  counts_matrix <- data.frame(t(as.matrix(GetAssayData(object = so_sub, slot = "counts"))),
                              check.names=F)
  counts_matrix <- counts_matrix[,colSums(counts_matrix) > 0] # Removing genes without counts
  counts_matrix <- cbind.data.frame(Cell = rownames(counts_matrix),
                                    counts_matrix)
  counts_matrix <- counts_matrix[,colnames(counts_matrix) %in% c("Cell", genes_to_keep)]
  colnames(counts_matrix)[1] <- ""
  write.table(counts_matrix, 
              paste0(output_nmf_counts_like_genes, sample, "_raw_counts.txt"), 
              sep = '\t', 
              row.names = F, 
              col.names = T, 
              quote = F)
}


####### -----            Run cNMF 

# Prepare params
outdir_low_iter <- paste0(output_nmf, "/NMF_by_sample_low_iter/")
k_range <- paste(as.character(seq(2,10)), collapse=" ") 
iterations <- 50 # 200 for real analysis
threads <- 4
seed <- 14
num_genes <- 2000
beta_loss <- "frobenius"

input_counts_folder <- output_nmf_counts_like_genes 

# Running cnmf commands prepare to k_selection_plot
# Need to select appropriate k's per sample after this loop
for(sample in unique(so$Sample)){
  print(sample)
  
  # Additional cnmf parameters
  input_file <- paste0(input_counts_folder, sample, "_raw_counts.txt")
  
  # cnmf prepare to k selection plot
  cmd_prepare <- paste0("cnmf prepare --output-dir ", outdir_low_iter,
                        " --name ", sample,
                        " -c ", input_file,
                        " -k ", k_range,
                        " --n-iter ", iterations,
                        " --total-workers ", threads,
                        " --seed ", seed,
                        " --numgenes ", num_genes,
                        " --beta-loss ", beta_loss) 
  system(cmd_prepare)
  cmd_factorize <- paste0("cnmf factorize --output-dir ", outdir_low_iter,
                          " --name ", sample,
                          " --worker-index ", 0)
  system(cmd_factorize)
  cmd_combine <- paste0("cnmf combine --output-dir ", outdir_low_iter,
                        " --name ", sample)
  system(cmd_combine)
  cmd_k_selection_plot <- paste0("cnmf k_selection_plot --output-dir ", outdir_low_iter,
                                 " --name ", sample)
  system(cmd_k_selection_plot)
}


## Next, running w/ 200 iterations and a single k selected by sample

# SELECT K PER SAMPLE BASED ON ABOVE PLOT
# Ideally want lowest error with high stability
samples <- unique(so$Sample)
k_per_sample <-c(3, 4, 3, 3, 5, 3, 4) 
# Prepare params
outdir <- paste0(output_nmf, "/NMF_by_sample/")
dir.create(outdir)

# Make sure to match all parameters to first run except k and iterations!
iterations <- 200 # 200 for real analysis

# Running cnmf commands prepare to k_selection_plot
for(i in 1:length(samples)){
  print(samples[i])
  sample <- samples[i]
  
  # Additional cnmf parameters
  input_file <- paste0(input_counts_folder, sample, "_raw_counts.txt")
  
  # cnmf prepare to k selection plot
  cmd_prepare <- paste0("cnmf prepare --output-dir ", outdir,
                        " --name ", sample,
                        " -c ", input_file,
                        " -k ", k_per_sample[i],
                        " --n-iter ", iterations,
                        " --total-workers ", threads,
                        " --seed ", seed,
                        " --numgenes ", num_genes,
                        " --beta-loss ", beta_loss) 
  system(cmd_prepare)
  cmd_factorize <- paste0("cnmf factorize --output-dir ", outdir,
                          " --name ", sample,
                          " --worker-index ", 0)
  system(cmd_factorize)
  cmd_combine <- paste0("cnmf combine --output-dir ", outdir,
                        " --name ", sample)
  system(cmd_combine)
  cmd_k_selection_plot <- paste0("cnmf k_selection_plot --output-dir ", outdir,
                                 " --name ", sample)
  system(cmd_k_selection_plot)
}

# Run cnmf consensus with the selected k WITHOUT additional filtering first
for(i in 1:length(samples)){
  print(samples[i])
  
  cmd_consensus <- paste0("cnmf consensus --output-dir ", outdir,
                          " --name ", samples[i],
                          " --components ", k_per_sample[i],
                          " --local-density-threshold ", 2, 
                          " --show-clustering ")
  system(cmd_consensus) # 066_Primary test used k = 6
  
  #     Setting the density threshold to anything >= 2.00 (the maximum possible 
  #     distance between two unit vectors) ensures that nothing will be filtered. "
}

# SELECT THE OUTLIER FILTER PER SAMPLE

# Only some samples had even minimal peaks of outliers
samples
outlier_cutoff <- c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01) 

# Run cnmf consensus with the selected k and outlier filter
for(i in 1:length(samples)){
  print(samples[i])
  
  cmd_consensus <- paste0("cnmf consensus --output-dir ", outdir,
                          " --name ", samples[i],
                          " --components ", k_per_sample[i],
                          " --local-density-threshold ", outlier_cutoff[i], 
                          " --show-clustering ")
  system(cmd_consensus) # 066_Primary test used k = 6
}


####### -----           Summarize cNMF results across sample

n_genes <- 50

top_programs <- cbind.data.frame(Program_holder = 1:n_genes) # Only top genes
genome_scores <- data.frame()
barcode_score <- cbind.data.frame(barcode = NA,
                                  sample = NA,
                                  program_cNMF_call = NA,
                                  program_cNMF_score = NA)
total_var_genes <- cbind.data.frame(sample = NA,
                                    var_genes = NA) # Curious how many variable genes are shared across tumor
for(i in 1:length(samples)){
  
  sample <- samples[i]
  print(sample)
  
  # Read in full genome scores
  scores <- read.table(file = paste0(outdir, sample, "/",
                                sample, ".gene_spectra_score.k_", k_per_sample[i], 
                                ".dt_", gsub("\\.", "_", outlier_cutoff[i]), ".txt"),
                       sep = "\t", header = T, stringsAsFactors = F, row.names = 1,
                       check.names=F)
  scores <- data.frame(t(scores),
                       check.names=F)
  colnames(scores) <- paste0(sample, ".", gsub("X", "program", colnames(scores)))
  print(nrow(scores)) # Some tumors could have 0 reads for genes
  
  # Saving full genome correlation per program
  scores_temp <- cbind.data.frame(Gene = rownames(scores),
                                  scores)
  if(nrow(genome_scores) < 2){
    genome_scores <- scores_temp
  }else{
    genome_scores <- join(genome_scores,
                          scores_temp)
  } # Saving after loop
  
  # Top correlated genes per program
  for(ii in 1:ncol(scores)){
    scores <- scores[order(scores[,ii], decreasing = T),]
    top_programs <- cbind(top_programs, rownames(scores)[1:n_genes])
    colnames(top_programs)[ncol(top_programs)] <- colnames(scores)[ii]
  } # Writing after loop
  
  # Read in scores of barcode by program
  barcode_scores <- read.table(file = paste0(outdir, sample, "/",
                                     sample, ".usages.k_", k_per_sample[i], 
                                     ".dt_", gsub("\\.", "_", outlier_cutoff[i]), ".consensus.txt"),
                       sep = "\t", header = T, stringsAsFactors = F, row.names = 1,
                       check.names=F)
  colnames(barcode_scores) <- paste0(sample, ".", gsub("X", "program", colnames(barcode_scores)))
  barcode_scores$program_cNMF_call <- NA
  barcode_scores$program_cNMF_score <- NA
  for(ii in 1:nrow(barcode_scores)){
    n <- as.numeric(which.max(barcode_scores[ii,]))
    barcode_scores$program_cNMF_call[ii] <- colnames(barcode_scores)[n]
    barcode_scores$program_cNMF_score[ii] <- as.numeric(barcode_scores[ii,n])
  }
  barcode_scores <- cbind.data.frame(barcode = rownames(barcode_scores),
                                     sample = sample,
                                     program_cNMF_call = barcode_scores$program_cNMF_call,
                                     program_cNMF_score = barcode_scores$program_cNMF_score)
  if(nrow(barcode_score) < 2){
    barcode_score <- barcode_scores
  }else{
    barcode_score <- rbind(barcode_score, barcode_scores)
  } # Writing after loop
  
  # Read in variable genes and save for each tumor
  var_genes <- read.table(file = paste0(outdir, sample, "/",
                              sample, ".overdispersed_genes.txt"),
                              sep = "\t", header = F, stringsAsFactors = F,
                          check.names=F)
  total_var_genes <- rbind(total_var_genes,
                           cbind.data.frame(sample = sample,
                                            var_genes = var_genes[,1]))
}

# Save top gene programs
top_programs <- top_programs[,-1] # Removing holder column
colnames(top_programs) <- paste0("GEP.", colnames(top_programs))
write.table(top_programs,
            file = paste0(output_nmf, "cNMF_top_correlated_genes_with_each_tumor_program.txt"),
            sep = "\t", row.names = F, quote = F)

genome_scores <- na.omit(genome_scores) # 12k remain
colnames(genome_scores)[2:ncol(genome_scores)] <- paste0("GEP.", colnames(genome_scores)[2:ncol(genome_scores)])
write.table(genome_scores,
            file = paste0(output_nmf, "cNMF_full_genome_correlation_within_each_sample.txt"),
            sep = "\t", row.names = F, quote = F)

# Save NMF-called programs 
barcode_score <- na.omit(barcode_score)
barcode_score$program_cNMF_call <- paste0("GEP.", barcode_score$program_cNMF_call)
write.table(barcode_score,
            file = paste0(output_nmf, "cNMF-calls_by_barcode.txt"),
            sep = "\t", row.names = F, quote = F)

# Inspect total variable genes
total_var_genes <- na.omit(total_var_genes)
table(total_var_genes$sample) # Number of variable genes per tumor (should be same)
freq_var_genes <- data.frame(unclass(table(total_var_genes$var_genes)),
                             check.names=F)
freq_var_genes <- cbind.data.frame(Gene = rownames(freq_var_genes),
                                   var_in_n_samples = freq_var_genes[,1])
freq_var_genes <- freq_var_genes[order(freq_var_genes$var_in_n_samples, decreasing = T),]
write.table(freq_var_genes,
            file = paste0(output_nmf, "cNMF_summary_variable_genes_per_sample.txt"),
            sep = "\t", row.names = F, quote = F)


## Downstream analysis using cell scores
output_nmf_merge_gene_scores <- paste0(output_nmf, "/merged_programs/")
dir.create(output_nmf_merge_gene_scores)

## Genes per program (used in next section)
top_genes_by_program <- read.table(file = paste0(output_nmf, "cNMF_top_correlated_genes_with_each_tumor_program.txt"),
                                      sep = "\t", header = T,
                                      check.names = F)
NMF_genes <- cbind.data.frame(Sample.NMF = NA,
                              Top_genes = NA)
for(i in 1:ncol(top_genes_by_program)){
  NMF_genes <- rbind(NMF_genes,
                     cbind.data.frame(Sample.NMF = colnames(top_genes_by_program)[i],
                                      Top_genes = top_genes_by_program[,i]))
}
NMF_genes <- na.omit(NMF_genes)

## Score each cell by the top gene loadings (FOR EACH INDIVIDUAL PROGRAM)
NMF_genes_list <- list()
for(NMF in unique(NMF_genes$Sample.NMF)){
  NMF_genes_list[[NMF]] <- NMF_genes$Top_genes[NMF_genes$Sample.NMF == NMF]
}

so <- AddModuleScore(
  so,
  features = NMF_genes_list,
  name = paste0(names(NMF_genes_list), "XXXX")
)
# Replace the '1' (and other numbers) added to the signature names
colnames(so@meta.data) <- gsub("XXXX\\d+$", "", colnames(so@meta.data))


## Score barcodes by only programs identified in that sample
# I labeled programs per cell within that sample’s called merged programs only. 
barcode_programs <- cbind.data.frame(barcode = NA,
                              top_sample_program = NA,
                              top_sample_program_score = NA) # Formally top_merge
for(sample in unique(so$Sample)){
  print(sample)
  
  so_sample <- subset(so, Sample == sample)
  sample_scores <- so_sample@meta.data[,grep(paste0("GEP.", sample), colnames(so_sample@meta.data))]
  # Z-score of each score
  for(i in 1:ncol(sample_scores)){
    sample_scores[,i] <- scale(sample_scores[,i], center = T, scale = T)
  }
  barcode_programs_sample <- cbind.data.frame(barcode = rownames(sample_scores),
                                       top_sample_program = NA,
                                       top_sample_program_score = NA)
  for(i in 1:nrow(sample_scores)){
    n <- as.numeric(which.max(sample_scores[i,]))
    barcode_programs_sample$top_sample_program[i] <- colnames(sample_scores)[n]
    barcode_programs_sample$top_sample_program_score[i] <- as.numeric(sample_scores[i,n])
  }
  barcode_programs <- rbind(barcode_programs,
                      barcode_programs_sample)
}
barcode_programs <- na.omit(barcode_programs)


## Cluster NMFs by cell scores

clusters <- 6 # I based this on the dendrogram and heatmap 

NMF_scores_cells <- so@meta.data[,grep("GEP", colnames(so@meta.data))]
colnames(NMF_scores_cells)

# Pairwise correlation between samples (columns)
cols.cor <- cor(NMF_scores_cells, use = "pairwise.complete.obs", method = "pearson")

# Plot the heatmap
ph <- pheatmap(cols.cor,
               fontsize = 8,
               clustering_distance_cols = as.dist(1 - cols.cor),
               clustering_distance_rows = as.dist(1 - cols.cor),
               cutree_cols = clusters,
               color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
)
save_pheatmap_pdf(ph,
                  paste0(output_nmf_merge_gene_scores, "Sample_program_correlation.pdf"),
                  width = 8)

## Merge gene loadings for each NMF group 
NMF_grouping <- data.frame(cutree(ph$tree_col, k=clusters))
NMF_grouping <- cbind.data.frame(merged_NMF = NMF_grouping[,1],
                                 original_NMF = rownames(NMF_grouping))
NMF_grouping <- data.frame(NMF_grouping[ph$tree_col$order,])
rownames(NMF_grouping) <- NMF_grouping$original_NMF
NMF_grouping$merged_NMF <- paste0("Merge.GEP.", NMF_grouping$merged_NMF)

# Only after further going through the script, I am re-labeling merged GEP names
gep_lookup <- read.table(file = paste0(output, "merge_GEP_lookup.txt"),
                         sep = "\t", header = T, stringsAsFactors = F)
for(i in 1:nrow(gep_lookup)){
  NMF_grouping$merged_NMF[NMF_grouping$merged_NMF == gep_lookup$Merge_orig[i]] <- gep_lookup$Merge_final[i]
}
table(NMF_grouping$merged_NMF)
write.table(NMF_grouping, paste0(output_nmf_merge_gene_scores, "Merged_NMF_groupings.txt"),
            sep = "\t", row.names = F, quote = F)

barcode_programs <- join(barcode_programs,
                         cbind.data.frame(top_sample_program = NMF_grouping$original_NMF,
                                          merged_program = NMF_grouping$merged_NMF))

merged_gene_loadings <- cbind.data.frame(NMF = NA,
                                         Top_genes = NA)
for(group in unique(NMF_grouping[,1])){

  group_samples <- rownames(NMF_grouping)[NMF_grouping[,1] == group]
  
  NMF_genes_group <- NMF_genes[paste0(NMF_genes$Sample.NMF) %in% group_samples,]
  
  # Printing number of genes in multiple samples over total genes
  print(paste0("Merge ", group, " has ", sum(duplicated(NMF_genes_group$Top_genes)), 
        " genes in multiple samples out of ", length(unique(NMF_genes_group$Top_genes)),
        " unique genes in program"))
  
  merged_gene_loadings <- rbind(merged_gene_loadings,
                                cbind.data.frame(NMF = group,
                                                 Top_genes = unique(NMF_genes_group$Top_genes))
  )
}
merged_gene_loadings <- na.omit(merged_gene_loadings)

# Check how many genes span merged loadings
sum(duplicated(merged_gene_loadings$Top_genes)) # 63 genes found in multiple merged programs
nrow(merged_gene_loadings) # 945
merged_gene_loadings <- merged_gene_loadings[!(merged_gene_loadings$Top_genes %in% 
                                                 unique(merged_gene_loadings$Top_genes[duplicated(merged_gene_loadings$Top_genes)])),]
nrow(merged_gene_loadings) # 822
# Reduced from this NMF grouping from 945 genes to 822
table(merged_gene_loadings$NMF)

write.table(merged_gene_loadings, paste0(output_nmf_merge_gene_scores, "Merged_gene_programs.txt"),
            sep = "\t", row.names = F, quote = F)

## ORA by the top gene loadings

all = msigdbr(species = "Homo sapiens")

gene_sets <- all[all$gs_cat %in% c("H") | all$gs_subcat %in% c("CP:KEGG", "GO:BP"),] %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
output_res_ora <- paste0(output_nmf_merge_gene_scores, "/ORA_merged/H_KEGG_GOBP/")
dir.create(output_res_ora, recursive = T)

gene_sets2 <- all[all$gs_cat %in% c("C3"),] %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
output_res_ora2 <- paste0(output_nmf_merge_gene_scores, "/ORA_merged/C3/")
dir.create(output_res_ora2, recursive = T)

gene_sets3 <- all[all$gs_cat %in% c("C6"),] %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
output_res_ora3 <- paste0(output_nmf_merge_gene_scores, "/ORA_merged/C6/")
dir.create(output_res_ora3, recursive = T)

gene_sets4 <- all[all$gs_cat %in% c("C8"),] %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
output_res_ora4 <- paste0(output_nmf_merge_gene_scores, "/ORA_merged/C8/")
dir.create(output_res_ora4, recursive = T)

gene_sets5 <- all[all$gs_cat %in% c("C2"),] %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
output_res_ora5 <- paste0(output_nmf_merge_gene_scores, "/ORA_merged/C2/")
dir.create(output_res_ora5, recursive = T)

gene_sets6 <- all %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
output_res_ora6 <- paste0(output_nmf_merge_gene_scores, "/ORA_merged/full/")
dir.create(output_res_ora6, recursive = T)

# Cycle through each cluster and perform ORA
clusters_final <- unique(merged_gene_loadings$NMF)
background <- rownames(so) # All gene names in dataset
for(specifNMF_cluster in clusters_final){
  print(specifNMF_cluster)
  
  genes <- merged_gene_loadings$Top_genes[merged_gene_loadings$NMF %in% specifNMF_cluster]
  
  # gene set 1
  ora_enricher <- enricher(genes, TERM2GENE=gene_sets, pAdjustMethod = "fdr", 
                           minGSSize = 5, maxGSSize = 1000, 
                           qvalueCutoff = 1, 
                           pvalueCutoff =1, # Usually don't but don't want to overwhelm
                           universe = background) 
  ora_result <- data.frame(ora_enricher)
  write.table(ora_result,
              file = paste0(output_res_ora,
                            specifNMF_cluster, ".txt"),
              row.names = F, sep = "\t", quote = F)
  
  # gene set 2
  ora_enricher <- enricher(genes, TERM2GENE=gene_sets2, pAdjustMethod = "fdr", 
                           minGSSize = 5, maxGSSize = 1000, 
                           qvalueCutoff = 1, 
                           pvalueCutoff =1, # Usually don't but don't want to overwhelm
                           universe = background) 
  ora_result <- data.frame(ora_enricher)
  write.table(ora_result,
              file = paste0(output_res_ora2,
                            specifNMF_cluster, ".txt"),
              row.names = F, sep = "\t", quote = F)
  
  # gene set 3
  ora_enricher <- enricher(genes, TERM2GENE=gene_sets3, pAdjustMethod = "fdr", 
                           minGSSize = 5, maxGSSize = 1000, 
                           qvalueCutoff = 1, 
                           pvalueCutoff =1, # Usually don't but don't want to overwhelm
                           universe = background) 
  ora_result <- data.frame(ora_enricher)
  write.table(ora_result,
              file = paste0(output_res_ora3,
                            specifNMF_cluster, ".txt"),
              row.names = F, sep = "\t", quote = F)
  
  # gene set 4
  ora_enricher <- enricher(genes, TERM2GENE=gene_sets4, pAdjustMethod = "fdr", 
                           minGSSize = 5, maxGSSize = 1000, 
                           qvalueCutoff = 1, 
                           pvalueCutoff =1, # Usually don't but don't want to overwhelm
                           universe = background) 
  ora_result <- data.frame(ora_enricher)
  write.table(ora_result,
              file = paste0(output_res_ora4,
                            specifNMF_cluster, ".txt"),
              row.names = F, sep = "\t", quote = F)
  
  # gene set 5
  ora_enricher <- enricher(genes, TERM2GENE=gene_sets5, pAdjustMethod = "fdr", 
                           minGSSize = 5, maxGSSize = 1000, 
                           qvalueCutoff = 1, 
                           pvalueCutoff =1, # Usually don't but don't want to overwhelm
                           universe = background) 
  ora_result <- data.frame(ora_enricher)
  write.table(ora_result,
              file = paste0(output_res_ora5,
                            specifNMF_cluster, ".txt"),
              row.names = F, sep = "\t", quote = F)
  
  # gene set 6
  ora_enricher <- enricher(genes, TERM2GENE=gene_sets6, pAdjustMethod = "fdr", 
                           minGSSize = 5, maxGSSize = 1000, 
                           qvalueCutoff = 1, 
                           pvalueCutoff =1, # Usually don't but don't want to overwhelm
                           universe = background) 
  ora_result <- data.frame(ora_enricher)
  write.table(ora_result,
              file = paste0(output_res_ora6,
                            specifNMF_cluster, ".txt"),
              row.names = F, sep = "\t", quote = F)
}

## ORA by the top gene loadings for each INDIVIDUAL program

# Mostly as a check to confirm that sample-specific programs I am merging
# have similar gene set analysis

all = msigdbr(species = "Homo sapiens")

gene_sets <- all[all$gs_cat %in% c("H") | all$gs_subcat %in% c("CP:KEGG", "GO:BP"),] %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
output_res_ora <- paste0(output_nmf, "/sample_program_ORA/H_KEGG_GOBP/")
dir.create(output_res_ora, recursive = T)

gene_sets2 <- all[all$gs_cat %in% c("C3"),] %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
output_res_ora2 <- paste0(output_nmf, "/sample_program_ORA/C3/")
dir.create(output_res_ora2, recursive = T)

gene_sets3 <- all[all$gs_cat %in% c("C6"),] %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
output_res_ora3 <- paste0(output_nmf, "/sample_program_ORA/C6/")
dir.create(output_res_ora3, recursive = T)

gene_sets4 <- all[all$gs_cat %in% c("C8"),] %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
output_res_ora4 <- paste0(output_nmf, "/sample_program_ORA/C8/")
dir.create(output_res_ora4, recursive = T)

gene_sets5 <- all[all$gs_cat %in% c("C2"),] %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
output_res_ora5 <- paste0(output_nmf, "/sample_program_ORA/C2/")
dir.create(output_res_ora5, recursive = T)

# Cycle through each cluster and perform ORA
background <- rownames(so) # All gene names in dataset
for(i in 1:length(NMF_genes_list)){
  print(i)
  
  genes <- NMF_genes_list[[i]]
  specifNMF_cluster <- names(NMF_genes_list)[i]
  
  # gene set 1
  ora_enricher <- enricher(genes, TERM2GENE=gene_sets, pAdjustMethod = "fdr", 
                           minGSSize = 5, maxGSSize = 1000, 
                           qvalueCutoff = 1, 
                           pvalueCutoff =1, # Usually don't but don't want to overwhelm
                           universe = background) 
  ora_result <- data.frame(ora_enricher)
  write.table(ora_result,
              file = paste0(output_res_ora,
                            specifNMF_cluster, ".txt"),
              row.names = F, sep = "\t", quote = F)
  
  # gene set 2
  ora_enricher <- enricher(genes, TERM2GENE=gene_sets2, pAdjustMethod = "fdr", 
                           minGSSize = 5, maxGSSize = 1000, 
                           qvalueCutoff = 1, 
                           pvalueCutoff =1, # Usually don't but don't want to overwhelm
                           universe = background) 
  ora_result <- data.frame(ora_enricher)
  write.table(ora_result,
              file = paste0(output_res_ora2,
                            specifNMF_cluster, ".txt"),
              row.names = F, sep = "\t", quote = F)
  
  # gene set 3
  ora_enricher <- enricher(genes, TERM2GENE=gene_sets3, pAdjustMethod = "fdr", 
                           minGSSize = 5, maxGSSize = 1000, 
                           qvalueCutoff = 1, 
                           pvalueCutoff =1, # Usually don't but don't want to overwhelm
                           universe = background) 
  ora_result <- data.frame(ora_enricher)
  write.table(ora_result,
              file = paste0(output_res_ora3,
                            specifNMF_cluster, ".txt"),
              row.names = F, sep = "\t", quote = F)
  
  # gene set 4
  ora_enricher <- enricher(genes, TERM2GENE=gene_sets4, pAdjustMethod = "fdr", 
                           minGSSize = 5, maxGSSize = 1000, 
                           qvalueCutoff = 1, 
                           pvalueCutoff =1, # Usually don't but don't want to overwhelm
                           universe = background) 
  ora_result <- data.frame(ora_enricher)
  write.table(ora_result,
              file = paste0(output_res_ora4,
                            specifNMF_cluster, ".txt"),
              row.names = F, sep = "\t", quote = F)
  
  # gene set 5
  ora_enricher <- enricher(genes, TERM2GENE=gene_sets5, pAdjustMethod = "fdr", 
                           minGSSize = 5, maxGSSize = 1000, 
                           qvalueCutoff = 1, 
                           pvalueCutoff =1, # Usually don't but don't want to overwhelm
                           universe = background) 
  ora_result <- data.frame(ora_enricher)
  write.table(ora_result,
              file = paste0(output_res_ora5,
                            specifNMF_cluster, ".txt"),
              row.names = F, sep = "\t", quote = F)
}


## Score cells with merged programs

merged_gene_loadings <- read.table(file = paste0(output_nmf_merge_gene_scores, 
                                            "Merged_gene_programs.txt"),
                                   sep = "\t", header = T, stringsAsFactors = F)
colnames(merged_gene_loadings)[1] <- "merge"

merge_genes_list <- list()
for(merge in unique(merged_gene_loadings$merge)){
  merge_genes_list[[merge]] <- merged_gene_loadings$Top_genes[merged_gene_loadings$merge == merge]
}

so@assays$RNA@data[10:20,10:20]

so <- AddModuleScore(
  so,
  features = merge_genes_list,
  name = paste0(names(merge_genes_list), "XXXX")
)
# Replace the '1' (and other numbers) added to the signature names
colnames(so@meta.data) <- gsub("XXXX\\d+$", "", colnames(so@meta.data))

# Z-score of each score
for(merge_gep in colnames(so@meta.data)[grep("Merge.GEP", colnames(so@meta.data))]){
  print(merge_gep)
  so@meta.data[,merge_gep] <- scale(so@meta.data[,merge_gep],
                                                  center = T, scale = T)
}

# Add program labels to so
barcode_programs <- barcode_programs[match(colnames(so),
                                           barcode_programs$barcode),]
so$sample_program <- barcode_programs$top_sample_program
so$merged_program <- barcode_programs$merged_program


## Merged scores in integrated UMAP

output_nmf_merge_gene_scores_umap <- paste0(output_nmf_merge_gene_scores,
                                       "UMAP_violin_merged_scores/")
dir.create(output_nmf_merge_gene_scores_umap)

for(mp in unique(so$merged_program)){
  ggf <- FeaturePlot(so,
                     features = mp,
                     order = T,
                     cols = c("lightgrey", "lightblue", "orange", "red"))
  ggv <- VlnPlot(so,
                 features = mp,
                 group.by = "merged_program")
  
  # Aggregate all plots
  gg_arr_gene <- grid.arrange(ggf,
                              ggv,
                              ncol = 2, nrow = 1)
  ggsave(paste0(output_nmf_merge_gene_scores_umap, mp, ".pdf"),
         gg_arr_gene)
}

## Merge call UMAP

gg_merged_call <- DimPlot(so,
                          group.by = "merged_program",
                          reduction = "umap")
gg_merged_call
ggsave(filename = paste0(output_nmf_merge_gene_scores, "UMAP_merged_call.pdf"),
       gg_merged_call)


####### -----            Save RDS

saveRDS(so,
        "Data/Seurat_objects/AggrNoNorm_filtered_primary_celltype_CNAsubs_integrateP_clustered_pathway_NMF.RDS")


####### -----            Additional figures
so <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_primary_celltype_CNAsubs_integrateP_clustered_pathway_NMF.RDS")

## Merged GEP scores by sample

output_nmf_merge_gene_scores_umap <- paste0(output_nmf_merge_gene_scores,
                                            "violin_merged_scores_by_sample/")
dir.create(output_nmf_merge_gene_scores_umap)
for(mp in unique(so$merged_program)){
  ggv <- VlnPlot(so,
                 features = mp,
                 group.by = "Sample")
  ggsave(paste0(output_nmf_merge_gene_scores_umap, mp, ".pdf"),
         ggv)
}

## Merged GEP scores by sample
output_nmf_merge_gene_scores_umap <- paste0(output_nmf_merge_gene_scores,
                                            "violin_QC_by_sample_gep_label/")
dir.create(output_nmf_merge_gene_scores_umap)
gg_umi <- VlnPlot(so,
              features = "nCount_RNA",
              split.by = "merged_program",
              group.by = "Sample") + xlab("Sample") + 
        ggtitle("Number of UMIs") + ylab("# of UMIs") +
        geom_vline(xintercept=(1:(length(unique(so$Sample))-1)+0.5))
gg_genes <- VlnPlot(so,
                  features = "nFeature_RNA",
                  split.by = "merged_program",
                  group.by = "Sample") + xlab("Sample") + 
  ggtitle("Number of Genes") + ylab("# of Genes") +
  geom_vline(xintercept=(1:(length(unique(so$Sample))-1)+0.5))
gg_mt <- VlnPlot(so,
                  features = "percent_mito",
                  split.by = "merged_program",
                  group.by = "Sample") + xlab("Sample") + 
  ggtitle("UMIs from mitochondria") + ylab("MT UMI (%)") +
  geom_vline(xintercept=(1:(length(unique(so$Sample))-1)+0.5))

ggarrange()
ggsave(paste0(output_nmf_merge_gene_scores_umap, mp, ".pdf"),
         ggv)

aggr_qc <- grid.arrange(gg_umi,
                        gg_genes,
                        gg_mt,
                        ncol = 1, nrow = 3)
ggsave(paste0(output_nmf_merge_gene_scores_umap, "/Aggregated_QC_by_GEP.pdf"),
       aggr_qc,
       width = 10, height = 10)



## Stacked barplot of merged calls and sample

prop_df_sample <- data.frame(prop.table(table(so@meta.data[,"merged_program"], 
                                              so@meta.data[,"Sample"]), margin = 2))
gg_stackbar_sample <- ggplot(prop_df_sample, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="stack", stat="identity") +ylab("Proportion") + xlab("Sample") +
  theme_classic()
gg_stackbar_sample
ggsave(paste0(output_nmf_merge_gene_scores, "/", "Sample_stacked_barplot_by_merged_program.pdf"),
       gg_stackbar_sample)

# Same but styled by proportion heatmap
source("functions/prop_hm.R")
dir.create(paste0(output_nmf_merge_gene_scores, "proportions_merged_groups_by_sample/"))
prop_hm(so = so,
        group1 = "Sample",
        group2 = "merged_program",
        outdir = paste0(output_nmf_merge_gene_scores, "proportions_merged_groups_by_sample/"),
        breaks = seq(0, 1, by = 0.01), # Sets max color to 1
        height = 4,
        width = 8
)

prop_df_sample <- data.frame(prop.table(table(so@meta.data[,"Sample"], 
                                              so@meta.data[,"merged_program"]), margin = 2))
gg_stackbar_sample <- ggplot(prop_df_sample, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="stack", stat="identity") +ylab("Proportion") + xlab("Merged IC") +
  theme_classic()
gg_stackbar_sample
ggsave(paste0(output_nmf_merge_gene_scores, "/", "Merged_program_stacked_barplot_by_Sample.pdf"),
       gg_stackbar_sample)


####### -----            Heatmap of the top gene loadings

merged_gene_loadings_order <- merged_gene_loadings[order(merged_gene_loadings$merge),]

DefaultAssay(so) <- "RNA"
so_scale <- ScaleData(so,
                      features = unique(merged_gene_loadings_order$Top_genes))
Idents(so_scale) <- so_scale$merged_program
table(Idents(so_scale))
gghm <- DoHeatmap(subset(so_scale, downsample = 100),
                  group.by=NULL,
                  label = F,
                  features = unique(merged_gene_loadings_order$Top_genes)
)
gghm
ggsave(paste0(output_nmf_merge_gene_scores, "Heatmap_downsample100.pdf"),
       gghm)

####### -----            Dotplot of select genes
ggd <- DotPlot(so,
               group.by = "merged_program",
               features = list("1" = c("JUN", "EGR3", "RHOB"),
                               "2" = c("KDSR", "PKP1", "SIRPA"),
                               "3" = c("TOP2A", "MKI67", "PCNA"),
                               "4" = c("HSP90B1", "HSPA5", "HYOU1"),
                               "5" = c( "RPL8", "RPS18", "RPL12"),
                               "6" = c("ATP5F1A", "NDUFB2", "COX6C")
               )) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggd
ggsave(paste0(output_nmf_merge_gene_scores, "Dotplot_top_select_genes.pdf"),
       ggd)


####### -----            Heatmap of select genes

select_genes <- unlist(list("1" = c("JUN", "EGR3", "RHOB"),
                            "2" = c("KDSR", "PKP1", "SIRPA"),
                            "3" = c("TOP2A", "MKI67", "PCNA"),
                            "4" = c("HSP90B1", "HSPA5", "HYOU1"),
                            "5" = c( "RPL8", "RPS18", "RPL12"),
                            "6" = c("ATP5F1A", "NDUFB2", "COX6C")))

so_dge_scale <- ScaleData(so,
                          features = select_genes
                ) 
Idents(so_dge_scale) <- so_dge_scale$merged_program
gghm_select <- DoHeatmap(subset(so_dge_scale, downsample = 100), 
                      features = select_genes) + NoLegend()
gghm_select
ggsave(paste0(output_nmf_merge_gene_scores, "Merged_program_heatmap_select_genes.pdf"),
       gghm_select)

gghm_select_w_legend <- DoHeatmap(subset(so_dge_scale, downsample = 100), 
                         features = select_genes)
gghm_select_w_legend
ggsave(paste0(output_nmf_merge_gene_scores, 
              "Merged_program_heatmap_select_genes_w_legend.pdf"),
       gghm_select_w_legend)

##########     ORA Heatmap

## Select ORA for merged program 

merged_gene_loadings <- read.table(file = paste0(output_nmf_merge_gene_scores, 
                                                 "Merged_gene_programs.txt"),
            sep = "\t", header = T)

all = msigdbr(species = "Homo sapiens")

selected_gene_sets <- c(
  # 1
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "PDGF_ERK_DN.V1_DN",
  "STK33_UP",
  
  # 2
  "RIGGI_EWING_SARCOMA_PROGENITOR_UP",
  "ZHANG_TARGETS_OF_EWSR1_FLI1_FUSION",
  "MIYAGAWA_TARGETS_OF_EWSR1_ETS_FUSIONS_UP",
  
  # 3 
  "HALLMARK_E2F_TARGETS",
  "GOBP_MITOTIC_CELL_CYCLE",

  # 4
  "REACTOME_UNFOLDED_PROTEIN_RESPONSE_UPR",
  "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
  
  # 4
  "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
  "WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS",
  "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",
  "KEGG_RIBOSOME",
  
  # 6
  "KEGG_OXIDATIVE_PHOSPHORYLATION",
  "GOBP_ATP_METABOLIC_PROCESS"

)

gene_sets <- all[all$gs_name %in% selected_gene_sets,] %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
output_res_ora <- paste0(output_nmf_merge_gene_scores, "/ORA_merged_select_gene_sets/")
dir.create(output_res_ora, recursive = T)

# Cycle through each cluster and perform ORA
clusters_final <- unique(merged_gene_loadings$NMF)
background <- rownames(so) # All gene names in dataset
combined_pval <- cbind.data.frame(GeneSet = selected_gene_sets,
                                  holder = NA)
for(specifNMF_cluster in clusters_final){
  print(specifNMF_cluster)
  
  genes <- merged_gene_loadings$Top_genes[merged_gene_loadings$NMF %in% specifNMF_cluster]
  
  # gene set 1
  ora_enricher <- enricher(genes, TERM2GENE=gene_sets, pAdjustMethod = "fdr", 
                           minGSSize = 1, maxGSSize = 10000, 
                           qvalueCutoff = 1, 
                           pvalueCutoff =1, # Usually don't but don't want to overwhelm
                           universe = background) 
  ora_result <- data.frame(ora_enricher)
  write.table(ora_result,
              file = paste0(output_res_ora,
                            specifNMF_cluster, ".txt"),
              row.names = F, sep = "\t", quote = F)
  
  combined_pval <- join(combined_pval,
                        cbind.data.frame(GeneSet = ora_result$ID,
                                         pval = -log10(ora_result$p.adjust)))
  colnames(combined_pval)[ncol(combined_pval)] <- specifNMF_cluster
  
}

rownames(combined_pval) <- combined_pval$GeneSet
combined_pval <- combined_pval[,-c(1,2)]
combined_pval[is.na(combined_pval)] <- 0 # If it wasn't tested, there were no overlapping genes
combined_pval <- combined_pval[,c(colnames(combined_pval)[1:ncol(combined_pval)][order(colnames(combined_pval)[1:ncol(combined_pval)])]
                                  )]
# REORDERING - CHANGE FOR SPECIFIC CLUSTERING
combined_pval <- combined_pval[,c(2,1,4,6,5,3)]

# Plot as heatmap
ph <- pheatmap(combined_pval,
               fontsize = 8,
               cluster_cols = F,
               cluster_rows = F,
               color = colorRampPalette(c("white", "orange", "darkorange", "red", "firebrick3", "darkred"))(100)
)
save_pheatmap_pdf(ph,
                  paste0(output_nmf_merge_gene_scores, "Merge_ORA_select_genesets.pdf"),
                  width = 8)


## Merge score grid featureplot
# https://github.com/satijalab/seurat/issues/1080
featureplot_celltype_final <- FeaturePlot(so,
                                          features = colnames(combined_pval),
                                          combine = F,
                                          order = T,
                                          cols = c("lightgrey", "lightblue", "orange", "red"))
for(i in 1:length(featureplot_celltype_final)) {
  featureplot_celltype_final[[i]] <- featureplot_celltype_final[[i]] + 
    NoLegend() + NoAxes() + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
}
featureplot_celltype_final_grid <- cowplot::plot_grid(plotlist = featureplot_celltype_final)
ggsave(paste0(output_nmf_merge_gene_scores, "UMAP_grid_merge_score.PDF"), 
       plot = print(featureplot_celltype_final_grid),
       height = 5,
       width = 7)

featureplot_celltype_final <- FeaturePlot(so,
                                          features = unique(so$merged_program),
                                          combine = F,
                                          order = T,
                                          cols = c("lightgrey", "lightblue", "orange", "red"))
for(i in 1:length(featureplot_celltype_final)) {
  featureplot_celltype_final[[i]] <- featureplot_celltype_final[[i]] + 
    # NoLegend() + 
    NoAxes() + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
}
featureplot_celltype_final_grid <- cowplot::plot_grid(plotlist = featureplot_celltype_final)
ggsave(paste0(output_nmf_merge_gene_scores, "UMAP_grid_merge_score_w_legend.PDF"), 
       plot = print(featureplot_celltype_final_grid),
       height = 5,
       width = 7)


#######    Calculate transcription factor activity across merged programs
# https://rdrr.io/github/christianholland/dorothea/f/vignettes/single_cell_vignette.Rmd#google_vignette

# Output directory
output_tfs <- paste0(output_nmf_merge_gene_scores, "merged_programs_TFs_by_VIPER/")
dir.create(output_tfs)

# Read Dorothea Regulons for Human:
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
length(unique(dorothea_regulon_human$tf)) # 1,333 TFs prior to filtering

# Obtain the regulons based on interactions with confidence level A, B and C
tf_confidence <- dorothea_regulon_human[!duplicated(dorothea_regulon_human$tf), c("tf", "confidence")]
tf_confidence$tf <- gsub("-", ".", tf_confidence$tf)
regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C"))
length(unique(regulon$tf))

# Compute Viper Scores 
so_viper <- so
DefaultAssay(so_viper) <- "RNA"
Idents(so_viper) <- so_viper$merged_program
so_viper <- run_viper(so_viper, regulon,
                  options = list(method = "scale", minsize = 4, 
                                 eset.filter = FALSE, cores = 1, 
                                 verbose = FALSE))
DefaultAssay(so_viper) <- "dorothea"
so_viper <- ScaleData(so_viper,
                      features = rownames(so_viper))
viper_scores_df <- GetAssayData(so_viper, slot = "scale.data", 
                                assay = "dorothea") %>%
  data.frame() %>%
  t()
# Create a data frame containing the cells and their clusters
CellsClusters <- data.frame(cell = names(Idents(so_viper)), 
                            cell_type = as.character(Idents(so_viper)),
                            stringsAsFactors = FALSE)
CellsClusters$cell <- gsub("-", ".", CellsClusters$cell) # The next step introduced periods
# Create a data frame with the Viper score per cell and its clusters
viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)
# Summarize the Viper scores by cell population
summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))
summarized_viper_scores <- join(summarized_viper_scores,
        tf_confidence) # Add confidence

# MAYBE THE BEST WAY TO USE THIS FILE IS TO IDENTIFY THE TOP SCORE PER PROGRAM
# AND CALCULATE THE DIFFERENCE BETWEEN THE MEDIAN AND 2 HIGHEST
summary_df <- summarized_viper_scores
summary_df$highest <- ""
summary_df$diff_to_2nd_highest <- ""
summary_df$diff_to_median <- "" # median of rest (doesn't include the max group)
summary_df$diff_to_mean <- "" # mean of rest 
for(tf in unique(summary_df$tf)){
  summary_df_temp <- summary_df[summary_df$tf == tf,]
  
  X <- summary_df_temp$avg
  n <- length(unique(X))
  n_high <- which(X == sort(unique(X),partial=n)[n])
  n_2nd_high <- which(X == sort(unique(X),partial=n-1)[n-1])
  
  summary_df$highest[summary_df$tf == tf & summary_df$cell_type == summary_df_temp$cell_type[n_high]] <- "Max"
  summary_df$diff_to_2nd_highest[summary_df$tf == tf & summary_df$cell_type == summary_df_temp$cell_type[n_high]] <- summary_df_temp$avg[n_high] - summary_df_temp$avg[n_2nd_high]
  summary_df$diff_to_median[summary_df$tf == tf & summary_df$cell_type == summary_df_temp$cell_type[n_high]] <- summary_df_temp$avg[n_high] - median(summary_df_temp$avg[-n_high])
  summary_df$diff_to_mean[summary_df$tf == tf & summary_df$cell_type == summary_df_temp$cell_type[n_high]] <- summary_df_temp$avg[n_high] - mean(summary_df_temp$avg[-n_high])
}

write.table(summary_df,
            file = paste0(output_tfs, "TF_activity_summary.txt"),
            sep = "\t", row.names = F, quote = F)

top_tfs_per_program <- summary_df %>%
  group_by(cell_type) %>%
  slice_max(n = 3, order_by = diff_to_mean) 
avg_for_hm <- dcast(summarized_viper_scores[summarized_viper_scores$tf %in% 
                                              unique(top_tfs_per_program$tf),],
      tf ~ cell_type, value.var = "avg")
rownames(avg_for_hm) <- avg_for_hm$tf
avg_for_hm <- avg_for_hm[,-1]
avg_for_hm <- avg_for_hm[match(unique(unique(top_tfs_per_program$tf)),
                               rownames(avg_for_hm)),]
pheatmap(avg_for_hm,
         cluster_cols = F,
         cluster_rows = F,
         scale = "none")
# Probably shouldn't scale because it is already the average of the scale
palette_length <- 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
my_breaks <- c(seq(min(avg_for_hm), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(avg_for_hm)/palette_length, 
                   max(avg_for_hm), 
                   length.out=floor(palette_length/2)))
ph <- pheatmap(avg_for_hm,
         cluster_cols = F,
         cluster_rows = F,
         scale = "none",
         color=my_color, breaks = my_breaks, 
         main = "DoRothEA (ABC)", angle_col = 45)
ph
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(ph,
                  paste0(output_tfs, "Heatmap_top_3_TFs_per_program.pdf"))


################################################################################
####### -------              Work with Aynaud gene sets          ------- #######


#######    Summarize Aynaud gene sets

aynaud_ic_summary <- read.table(file = "helper_data/IC_ngenes_phenotype.txt",
                                sep = "\t", header = T, stringsAsFactors = F)
aynaud_full_weights <- read.table(file = "helper_data/aggregated_1964_log_ica_S.txt",
                                  sep = "\t", header = T, stringsAsFactors = F)
aynaud_genesets <- cbind.data.frame(IC = NA,
                                    Gene = NA)
for(ic in unique(aynaud_ic_summary$IC)){
  ic_ngene <- aynaud_ic_summary$n_genes[aynaud_ic_summary$IC == ic]
  aynaud_full_weights_temp <- aynaud_full_weights[order(aynaud_full_weights[,ic], decreasing = T),]
  aynaud_genesets <- rbind(aynaud_genesets,
                           cbind.data.frame(IC = ic,
                                            Gene = aynaud_full_weights_temp$PROBE[1:ic_ngene]))
}
aynaud_genesets <- na.omit(aynaud_genesets)
aynaud_genesets <- join(aynaud_genesets, aynaud_ic_summary)
write.table(aynaud_genesets,
            file = "helper_data/aynaud_genesets.txt",
            sep = "\t", row.names = F, quote = F)


#######    Correlation between nmf and aynaud programs 

# Phenotype-labled and IC7 and IC13 

output_aynaud_cor <- paste0(output_nmf_merge_gene_scores, "/correlation_w_aynaud_v2/")
dir.create(output_aynaud_cor)

so <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_primary_celltype_CNAsubs_integrateP_clustered_pathway_NMF.RDS")
DefaultAssay(so)

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Gene lists
aynaud_genesets <- read.table(file = "Data/helper_data/aynaud/aynaud_genesets.txt", 
                              sep = "\t", header = T)
# Subsetting aynaud to those with pheontype (and also IC7 and IC13)
aynaud_genesets_sub <- aynaud_genesets
aynaud_genesets_sub$phenotype[aynaud_genesets_sub$IC %in% c("IC7", "IC13")] <- 
  aynaud_genesets_sub$IC[aynaud_genesets_sub$IC %in% c("IC7", "IC13")]
aynaud_genesets_sub <- aynaud_genesets_sub[!(aynaud_genesets_sub$phenotype == ""),]
aynaud_genesets_list <- split(aynaud_genesets_sub$Gene,
                              aynaud_genesets_sub$phenotype)

so <- AddModuleScore(
  so,
  features = aynaud_genesets_list,
  name = paste0(names(aynaud_genesets_list), "XXXX")
)
# Replace the '1' (and other numbers) added to the signature names
colnames(so@meta.data) <- gsub("XXXX\\d+$", "", colnames(so@meta.data))

# Z-score of each score
for(ic in names(aynaud_genesets_list)){
  print(ic)
  so@meta.data[,ic] <- scale(so@meta.data[,ic], center = T, scale = T)
}

# Extract all scores
scores <- so@meta.data[,c(names(nmf_programs_list), names(aynaud_genesets_list))]

myBreaks <- c(seq(-1, 0, length.out=ceiling(50/2) + 1), 
              seq(1/50, 1, length.out=floor(50/2)))

# Pairwise correlation - all gene sets
cols.cor <- cor(scores, use = "pairwise.complete.obs", method = "pearson")
ph <- pheatmap(cols.cor,
               fontsize = 8,
               clustering_distance_cols = as.dist(1 - cols.cor),
               clustering_distance_rows = as.dist(1 - cols.cor),
               color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
               breaks=myBreaks,
               display_numbers = round(cols.cor,2),
               fontsize_number = 8,
               number_color = "black"
)
save_pheatmap_pdf(ph, width = 7,
                  filename = paste0(output_aynaud_cor, "NMF_and_Aynaud_programs.pdf"))

# Pairwise correlation - NMF-only
cols.cor <- cor(scores[,names(nmf_programs_list)], 
                use = "pairwise.complete.obs", method = "pearson")
ph <- pheatmap(cols.cor,
               fontsize = 8,
               clustering_distance_cols = as.dist(1 - cols.cor),
               clustering_distance_rows = as.dist(1 - cols.cor),
               #cutree_cols = clusters,
               color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
               breaks=myBreaks,
               display_numbers = round(cols.cor,2),
               fontsize_number = 8,
               number_color = "black"
)
save_pheatmap_pdf(ph, width = 5, height = 5,
                  filename = paste0(output_aynaud_cor, "NMF-only.pdf"))

# Pairwise correlation - Aynaud-only
cols.cor <- cor(scores[,names(aynaud_genesets_list)], 
                use = "pairwise.complete.obs", method = "pearson")
ph <- pheatmap(cols.cor,
               fontsize = 8,
               clustering_distance_cols = as.dist(1 - cols.cor),
               clustering_distance_rows = as.dist(1 - cols.cor),
               #cutree_cols = clusters,
               color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
               breaks=myBreaks,
               display_numbers = round(cols.cor,2),
               fontsize_number = 8,
               number_color = "black"
)
save_pheatmap_pdf(ph, width = 5, height = 5,
                  filename = paste0(output_aynaud_cor, "Aynaud-only.pdf"))

# Correlate only NMF to Aynaud
# https://stackoverflow.com/questions/38631223/compare-only-pairwise-comparisons-between-two-matrices-using-rcorr-in-r
matrix1 <- scores[,names(nmf_programs_list)]
matrix2 <- scores[,names(aynaud_genesets_list)]
tmp <- with(expand.grid(seq(ncol(matrix1)), seq(ncol(matrix2))),
            mapply(function(i, j) cor.test(matrix1[, i], matrix2[, j]),
                   Var1, Var2))
cor_matrix <- matrix(unlist(tmp['estimate', ]), nrow=ncol(matrix1),
                dimnames=list(colnames(matrix1), colnames(matrix2)))
ph <- pheatmap(cor_matrix,
               fontsize = 8,
               #clustering_distance_cols = as.dist(1 - cor_matrix),
               #clustering_distance_rows = as.dist(1 - cor_matrix),
               #cutree_cols = clusters,
               color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
               breaks=myBreaks,
               display_numbers = round(cor_matrix,2),
               fontsize_number = 8,
               number_color = "black"
)
save_pheatmap_pdf(ph, width = 6, height = 5,
                  filename = paste0(output_aynaud_cor, "NMF_vs_Aynaud.pdf"))
