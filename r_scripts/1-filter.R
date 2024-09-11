
################################################################################
####### -------                   Goals                          ------- #######

# Explore QC
# Filter the Seurat object

################################################################################
####### -------        Notes on this analysis                    ------- #######
# Final filtering: 
    # # Remove cells with <200 detected genes
    # # Remove cells with >6000 detected genes (could be doublets)
    # # Remove cells with >10% mitochondrial reads
    # # Remove cells with > 50000 UMIs

# Dataset specific
# NA values will not show any cutoffs on final plots
umi_low_cut <- NA
umi_high_cut <- 50000
mt_high_cut <- 20
gene_low_cut <- 200
gene_high_cut <- 7000

# Vectors for plotting
umi_cut_vector <- c(umi_low_cut, umi_high_cut)
umi_cut_vector[is.na(umi_cut_vector)] <- 0

gene_cut_vector <- c(gene_low_cut, gene_high_cut)
gene_cut_vector[is.na(gene_cut_vector)] <- 0

# Other plotting params
gg_title_text <- 12
gg_axis_text <- 12


################################################################################
####### -------                Read in data                      ------- #######
library(dplyr)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(scDblFinder)
library(plotly)

source("functions/SeuratQC.R")

setwd("/Volumes/Partition 2/Core/Hayashi_Ewing_scRNAseq_Dec2022/")
sobj <- readRDS("Data/Seurat_objects/AggrNoNorm_minFiltered.RDS")

output <- "Analysis/preFilter_QC"
dir.create(output)

################################################################################
####### -------       REMOVING TROUBLESOME SAMPLES         ------- #######

# For this dataset, cellranger QC identified 5 troublesome samples

table(sobj$Sample)
remove_samples <- c("006_Primary",
                    "009_Primary",
                    "025_Primary",
                    "051_CTC",
                    "066_CTC") 
sobj <- subset(sobj,
                Sample %in% remove_samples,
                invert = T)
table(sobj$Sample) # Confirm desired samples are removed


################################################################################
####### -------              Prefiltering table                  ------- #######

# Note that very minimal filtering occurred within aggr_to_seurat.R

# Dataset level
table_prefilter <- cbind.data.frame(Metric = c("Number of barcodes", 
                                             "Number of genes",
                                             "Median UMI per cell",
                                             "Median genes per cell"),
                                  Value = c(ncol(sobj),
                                            nrow(sobj),
                                            median(sobj$nCount_RNA),
                                            median(sobj$nFeature_RNA)))
write.table(table_prefilter,
            file = paste0(output, "/prefilter_dataset-level_metrics.txt"),
            sep = "\t", row.names = F, quote = F)

# Sample-level
table_sample <- cbind.data.frame(Sample = c("Number of barcodes", 
                                            "Number of genes",
                                            "Median UMI per cell",
                                            "Median genes per cell"))
for(samp in unique(sobj$Sample)){
  print(samp)
  
  sobj_samp <- subset(sobj,
                      Sample == samp)
  sobj_samp <- sobj_samp[rowSums(sobj_samp) > 0,]
  
  table_sample$samp <- c(ncol(sobj_samp),
                         nrow(sobj_samp),
                         median(sobj_samp$nCount_RNA),
                         median(sobj_samp$nFeature_RNA))
  colnames(table_sample)[ncol(table_sample)] <- samp
}

table_sample <- data.frame(t(table_sample))
write.table(table_sample,
            file = paste0(output, "/prefilter_sample-level_metrics.txt"),
            sep = "\t", row.names = T, col.names = F, quote = F)

################################################################################
####### -------         Add metadata including % MT              ------- #######

rownames(sobj)[grep("MT-", rownames(sobj))]

sobj <- sobj %>%
  PercentageFeatureSet(
    pattern  = "^MT-", 
    col.name = "percent_mito" 
  ) 

sobj$Sample2 <- paste0(sapply(strsplit(sobj$Sample, split = "_"),
                              function(x){x[[2]]}),
                       "_",
                       sapply(strsplit(sobj$Sample, split = "_"),
                              function(x){x[[1]]}))
sobj$Sample_type <- sapply(strsplit(gsub("CTC", "CTC-enriched", sobj$Sample), split = "_"),
                              function(x){x[[2]]}) # Will use at later point


################################################################################
####### -------                  Plotting                        ------- #######

gg_umi <- sobj %>%
  VlnPlot(
    features = c("nCount_RNA"),
    ncol     = 1,
    pt.size  = 0.01,
    split.by = "Sample_type",
    group.by = "Sample2"#,
  ) + geom_hline(yintercept=umi_cut_vector,
                 linetype = 2) + 
  ylab("# of UMIs") + ggtitle("Number of UMIs") + xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text=element_text(size = gg_title_text),
        axis.text = element_text(size = gg_axis_text))
gg_umi

gg_genes <- sobj %>%
  VlnPlot(
    features = c("nFeature_RNA"),
    ncol     = 1,
    pt.size  = 0.01,
    split.by = "Sample_type",
    group.by = "Sample2"#,
  ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            text=element_text(size = gg_title_text),
            axis.text = element_text(size = gg_axis_text)) + 
  ylab("# of genes") + ggtitle("Number of genes") + xlab("") +
  geom_hline(yintercept=gene_cut_vector,
             linetype = 2)
gg_genes

gg_MT <- sobj %>%
  VlnPlot(
    features = c("percent_mito"),
    ncol     = 1,
    pt.size  = 0.01,
    split.by = "Sample_type",
    group.by = "Sample2"#,
  )+ geom_hline(yintercept=mt_high_cut,
                linetype = 2) + 
  ylab("UMI (%)") + ggtitle("UMI from mitochondria") + xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text=element_text(size = gg_title_text),
        axis.text = element_text(size = gg_axis_text))
gg_MT

gg_umi_vs_genes <- FeatureScatter(object = sobj, 
                                  feature1 = 'nCount_RNA', 
                                  feature2 = 'nFeature_RNA',
                                  group.by = "Sample2",
                                  pt.size = 0.1)
gg_umi_vs_genes


gg_MT_vs_UMI <- FeatureScatter(object = sobj, 
                               feature2 = 'percent_mito', 
                               feature1 = 'nCount_RNA',
                               group.by = "Sample2",
                               pt.size = 0.1)

gg_MT_vs_UMI 



#########
# https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html

metadata <- sobj@meta.data

metadata$log10GenesPerUMI <- log10(metadata$nFeature_RNA) / log10(metadata$nCount_RNA)
sobj$log10GenesPerUMI <- log10(sobj$nFeature_RNA) / log10(sobj$nCount_RNA)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
hist_GenesPerUMI <- metadata %>% 
  ggplot(aes(color=Sample2, x=log10GenesPerUMI, fill= Sample2)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  ylab("Log10 Cell density") 

# Visualize the number UMIs/transcripts per cell
hist_UMI <- metadata %>% 
  ggplot(aes(color=Sample2, x=nCount_RNA, fill= Sample2)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Log10 UMI density") +
  geom_vline(xintercept = umi_cut_vector)

# Visualize the distribution of genes detected per cell via histogram
hist_genes <- metadata %>% 
  ggplot(aes(color=Sample2, x=nFeature_RNA, fill= Sample2)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = gene_cut_vector)

#####################################################################
#####    Save pre-filtering plots

ggsave(paste0(output, "/Genes.pdf"),
       gg_genes)
ggsave(paste0(output, "/UMIs.pdf"),
       gg_umi)
ggsave(paste0(output, "/MT.pdf"),
       gg_MT)
ggsave(paste0(output, "/Hist_UMIs.pdf"),
       hist_UMI)
ggsave(paste0(output, "/Hist_MT.pdf"),
       hist_MT)
ggsave(paste0(output, "/Hist_genes.pdf"),
       hist_genes)
ggsave(paste0(output, "/Hist_GenesPerUMI.pdf"),
       hist_GenesPerUMI)

sobj_publication_prefilter <- sobj

gg_umi_publication_prefilter <- sobj_publication_prefilter %>%
  VlnPlot(
    features = c("nCount_RNA"),
    ncol     = 1,
    pt.size  = 0,
    split.by = "Sample_type",
    group.by = "Sample2"#,
  ) + geom_hline(yintercept=umi_cut_vector,
                 linetype = 2) + 
  ylab("# of UMIs") + ggtitle("Number of UMIs") + xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text=element_text(size = gg_title_text),
        axis.text = element_text(size = gg_axis_text))
gg_umi_publication_prefilter

gg_genes_publication_prefilter <- sobj_publication_prefilter %>%
  VlnPlot(
    features = c("nFeature_RNA"),
    ncol     = 1,
    pt.size  = 0,
    split.by = "Sample_type",
    group.by = "Sample2"#,
  ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            text=element_text(size = gg_title_text),
            axis.text = element_text(size = gg_axis_text)) + 
  ylab("# of genes") + ggtitle("Number of genes") + xlab("") +
  geom_hline(yintercept=gene_cut_vector,
             linetype = 2)
gg_genes_publication_prefilter

gg_MT_publication_prefilter <- sobj_publication_prefilter %>%
  VlnPlot(
    features = c("percent_mito"),
    ncol     = 1,
    pt.size  = 0,
    split.by = "Sample_type",
    group.by = "Sample2"#,
  )+ geom_hline(yintercept=mt_high_cut,
                linetype = 2) + 
  ylab("UMI (%)") + ggtitle("UMI from mitochondria") + xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text=element_text(size = gg_title_text),
        axis.text = element_text(size = gg_axis_text))
gg_MT_publication_prefilter

aggr_vio <- grid.arrange(gg_umi_publication_prefilter,
                         gg_genes_publication_prefilter,
                         gg_MT_publication_prefilter,
                         ncol = 3, nrow = 1)
ggsave(paste0(output, "/Aggregated_violin_plots_pub_samples_only.pdf"),
       aggr_vio,
       width = 16, height = 4)


################################################################################
####### -------                Doublet finder                    ------- #######

sobj_type.list <- SplitObject(sobj, split.by = "Sample_type")

scDblFinder_results <- cbind.data.frame(Barcode = NULL,
                                        scDblFinder.class = NULL,
                                        scDblFinder.score = NULL)
for(i in 1:length(sobj_type.list)){
  so.sce <- as.SingleCellExperiment(sobj_type.list[[i]])
  so.sce <- scDblFinder(so.sce)
  scDblFinder_results_temp <- cbind.data.frame(Barcode = colnames(so.sce),
                                               scDblFinder.class = so.sce$scDblFinder.class,
                                               scDblFinder.score = so.sce$scDblFinder.score)
  scDblFinder_results <- rbind(scDblFinder_results, scDblFinder_results_temp)
}
scDblFinder_results <- scDblFinder_results[match(colnames(sobj),
                                                 scDblFinder_results$Barcode),]
sobj$scDblFinder.class <- scDblFinder_results$scDblFinder.class
sobj$scDblFinder.score <- scDblFinder_results$scDblFinder.score

################################################################################
####### -------                Filter cells                      ------- #######

sobj_metrics_filtered <- sobj %>%
  subset(
        nCount_RNA <= umi_high_cut & 
        nFeature_RNA >= gene_low_cut &   
        nFeature_RNA <= gene_high_cut & 
        percent_mito <= mt_high_cut 
    )

VlnPlot(sobj_metrics_filtered,
        features = "scDblFinder.score",
        group.by = "scDblFinder.class")

# Filter by singlet as well
sobj <- subset(sobj_metrics_filtered,
               scDblFinder.class == "singlet")

sobj_filtered <- sobj

saveRDS(sobj_filtered, 
        "Data/Seurat_objects/AggrNoNorm_filtered.RDS")


################################################################################
####### -------            Post filter QC plot                   ------- #######

output_postfilter <- "Analysis/postFilter_QC"
dir.create(output_postfilter,
           recursive = T)

SeuratQC(so = sobj_filtered,
         category = "Sample2",
         outdir = output_postfilter,
         violin_pt_size = 0.01)

sobj_filtered_publication <- sobj_filtered

table(sobj_publication_prefilter$Sample)

gg_umi_filtered_publication <- sobj_filtered_publication %>%
  VlnPlot(
    features = c("nCount_RNA"),
    ncol     = 1,
    pt.size  = 0,
    split.by = "Sample_type",
    group.by = "Sample2"#,
  ) + # geom_hline(yintercept=umi_cut_vector) + 
  ylab("# of UMIs") + ggtitle("Number of UMIs") + xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text=element_text(size = gg_title_text),
        axis.text = element_text(size = gg_axis_text))
gg_umi_filtered_publication

gg_genes_filtered_publication <- sobj_filtered_publication %>%
  VlnPlot(
    features = c("nFeature_RNA"),
    ncol     = 1,
    pt.size  = 0,
    split.by = "Sample_type",
    group.by = "Sample2"#,
  ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            text=element_text(size = gg_title_text),
            axis.text = element_text(size = gg_axis_text)) + 
  ylab("# of genes") + ggtitle("Number of genes") + xlab("") 
gg_genes_filtered_publication

gg_MT_filtered_publication <- sobj_filtered_publication %>%
  VlnPlot(
    features = c("percent_mito"),
    ncol     = 1,
    pt.size  = 0,
    split.by = "Sample_type",
    group.by = "Sample2"#,
  ) + # geom_hline(yintercept=mt_high_cut) + 
  ylab("UMI (%)") + ggtitle("UMI from mitochondria") + xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text=element_text(size = gg_title_text),
        axis.text = element_text(size = gg_axis_text))
gg_MT_filtered_publication

# Save an aggregated plot for ease of presentation
aggr_vio_publication <- grid.arrange(gg_umi_filtered_publication,
                         gg_genes_filtered_publication,
                         gg_MT_filtered_publication,
                         ncol = 3, nrow = 1)
ggsave(paste0(output_postfilter, "/Aggregated_violin_plots_post_filter_pub_samples_only.pdf"),
       aggr_vio_publication,
       width = 16, height = 4)



################################################################################
####### -------            Post filter tables                    ------- #######

# Dataset level
table_filter <- cbind.data.frame(Metric = c("Number of cells", 
                                               "Number of genes",
                                               "Median UMI per cell",
                                               "Median genes per cell"),
                                    Value = c(ncol(sobj),
                                              nrow(sobj),
                                              median(sobj$nCount_RNA),
                                              median(sobj$nFeature_RNA)))
write.table(table_filter,
            file = paste0(output_postfilter, "/dataset-level_metrics.txt"),
            sep = "\t", row.names = F, quote = F)

# Sample-level
table_sample <- cbind.data.frame(Sample = c("Number of cells", 
                                            "Number of genes",
                                            "Median UMI per cell",
                                            "Median genes per cell"))
for(samp in unique(sobj$Sample)){
  print(samp)
  
  sobj_samp <- subset(sobj,
                      Sample == samp)
  sobj_samp <- sobj_samp[rowSums(sobj_samp) > 0,]
  
  table_sample$samp <- c(ncol(sobj_samp),
                         nrow(sobj_samp),
                         median(sobj_samp$nCount_RNA),
                         median(sobj_samp$nFeature_RNA))
  colnames(table_sample)[ncol(table_sample)] <- samp
}

table_sample <- data.frame(t(table_sample))
write.table(table_sample,
            file = paste0(output_postfilter, "/sample-level_metrics.txt"),
            sep = "\t", row.names = T, col.names = F, quote = F)

