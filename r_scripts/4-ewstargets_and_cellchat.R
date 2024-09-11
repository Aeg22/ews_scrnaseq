

################################################################################
####### -------                    Goals                         ------- #######

# Processing the primary samples only
# Subclustered labels added to the primary dataset

# Dotplots of proposed EWS targets
# CellChat

################################################################################
####### -------                 Load R packages                  ------- #######

library(Seurat)
library(ggplot2)
library(plyr)
library(CellChat)
library(patchwork)
library(openxlsx)
options(stringsAsFactors = FALSE)
library(scales)


################################################################################
####### -------                Update variables                  ------- #######

setwd("/Volumes/Partition 2/Core/Hayashi_Ewing_scRNAseq_Dec2022/")

out_cellchat <- paste0("Analysis/Annotation/Primary/Normalized/cellchat/")
out_cellchat_objects <- paste0(out_cellchat, "objects/")
dir.create(out_cellchat_objects, recursive = T)

# For the Seurat object in the CellChat section:
# samples should be under 'Sample'
# cell type labels should be under 'Sub.Cell.Type'
cell_of_interest <- "Ewing sarcoma" # The main cell type I am interested in
data_is_unsorted <- FALSE 
min_cells_to_test <- 10

################################################################################
####### -------             Add subcluster labels to main        ------- #######

so <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_primary_celltype.RDS")
so_t <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_primary_celltype_CNAsubs_sublusteredTcells_pathway.RDS")
so_m <- readRDS("Data/Seurat_objects/AggrNoNorm_filtered_primary_celltype_CNAsubs_sublusteredMyeloid_pathway.RDS")

so_meta <- so@meta.data
so_meta <- cbind.data.frame(barcode = rownames(so_meta),
                 so_meta)

so_t_meta <- so_t@meta.data
so_t_meta <- cbind.data.frame(barcode = rownames(so_t_meta),
                              so_t_meta)

so_m_meta <- so_m@meta.data
so_m_meta <- cbind.data.frame(barcode = rownames(so_m_meta),
                              so_m_meta)

# Combined M and T
so_m_t_meta <- rbind(so_m_meta[,c("barcode", "Cell.Type", "Broad.Cell.Type")],
                     so_t_meta[,c("barcode", "Cell.Type", "Broad.Cell.Type")])
colnames(so_m_t_meta)[2:3] <- paste0("Sub.", colnames(so_m_t_meta)[2:3])

# Add to so
so_meta <- join(so_meta, so_m_t_meta)
so_meta$Sub.Cell.Type[which(is.na(so_meta$Sub.Cell.Type))] <- as.character(so_meta$Cell.Type[which(is.na(so_meta$Sub.Cell.Type))])
table(so_meta$Sub.Cell.Type,
      so_meta$Cell.Type)
table(so_meta$Sub.Cell.Type)
so_meta <- so_meta[match(colnames(so),
                         so_meta$barcode),]
so$Sub.Cell.Type <- so_meta$Sub.Cell.Type

so$Sub.Cell.Type <- gsub("Ewing sarcoma proliferating", "Ewing sarcoma", so$Sub.Cell.Type)

################################################################################
####### -------       Dotplots of proposed EWS targets           ------- #######

output_targets <- "Analysis/Annotation/Primary/Normalized/targets_plots/"
dir.create(output_targets)

####### -------                 Process objects
so2 <- so
tme_ews_targets <- c("CD274", "PAX3", "STEAP1", "PAPPA", rownames(so2)[grep("^HLA", rownames(so2))])
# so3 where I include all cell types but separate EWS by patient as well
so3 <- so2
so3$patEWS_w_type <- as.character(droplevels(so3$Cell.Type))
so3$patEWS_w_type[so3$patEWS_w_type == "Ewing sarcoma proliferating"] <- "Ewing sarcoma"
n <- which(so3$patEWS_w_type == "Ewing sarcoma")
so3$patEWS_w_type[n] <- paste0("EWS - ", gsub("_Primary", "", so3$Sample[n]))

####### -------                Curated EWS targets
mas_targets <- read.xlsx("helper_data/EWS_targets 12.2023_AG.xlsx")
mas_targets$EWS_target_gene[!(mas_targets$EWS_target_gene %in% rownames(so2))]
mas_targets$Class <- gsub("Surface", "Immunotherapeutic", mas_targets$Class)
mas_targets$Class <- gsub("Tyrosine kinase", "Signaling", mas_targets$Class)

for(class in unique(mas_targets$Class)){
  genes <- mas_targets$EWS_target_gene[mas_targets$Class == class]
  
  ggdot_ews_targets <- DotPlot(so2,
                               features = genes,
                               group.by = "Sub.Cell.Type",
                               cluster.idents = F) + 
    theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5)) + 
    ggtitle(paste0(class, "-related EWS targets by cell type"))
  ggsave(paste0(output_targets, "Literature_targets_", gsub(" ", "_", class),
                "_by_cell_type.PDF"), 
         plot = ggdot_ews_targets)
  
  # Sub to just EWS and plot across patient
  ggdot_ews_targets_pats <- DotPlot(so2_ews,
                                    features = genes,
                                    group.by = "Sample",
                                    cluster.idents = F) + 
    theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5)) + 
    ggtitle(paste0(class, "-related EWS targets by tumor"))
  ggsave(paste0(output_targets, "Literature_targets_", gsub(" ", "_", class),
                "_by_sample.PDF"), 
         plot = ggdot_ews_targets_pats)
  
  gg_temp <- DotPlot(so3,
          features = genes,
          group.by = "patEWS_w_type",
          cluster.idents = F) + 
    theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5)) + 
    ggtitle(paste0(class, "-related EWS targets by tumor"))
  ggsave(paste0(output_targets, "Literature_targets_", gsub(" ", "_", class),
                "_by_pat_and_celltype.PDF"),
         plot = gg_temp)
}

# Single, combined plot with different genes
tme_ews_targets2 <- tme_ews_targets[-which(tme_ews_targets == "STEAP1")]
gg_temp <- DotPlot(so3,
                   features = list("TME" = tme_ews_targets2,
                                   "Immunotherapeutic" = mas_targets$EWS_target_gene[mas_targets$Class == "Immunotherapeutic"],
                      "Signaling" = mas_targets$EWS_target_gene[mas_targets$Class == "Signaling"]),
                   group.by = "patEWS_w_type",
                   cluster.idents = F) + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
gg_temp
ggsave(paste0(output_targets, "combined_gene_sets_by_pat_and_celltype.PDF"),
       plot = gg_temp, width = 16, height = 6, units = "in")

# Potential immunotherapeutics
# https://aacrjournals.org/clincancerres/article/30/5/1022/734306/Surface-and-Global-Proteome-Analyses-Identify
# The prioritization generated 39 candidates in Group 1 (Z-scores ≥ 1) and 56 candidates in Group 2 (Z-scores > 0; Fig. 3C and D; Supplementary Data S3). 
it_df <- read.xlsx("helper_data/mooney_surfaceome_immune_score_ranking_modified_from_supp_data_s3.xlsx")
group1 <- it_df$symbol[it_df$Group == "Zscore>1"]
group2 <- it_df$symbol[it_df$Group == "Zscore>0"]
gg_temp <- DotPlot(so3,
                   features = list("Group 1" = group1,
                                   "Group 2" = group2),
                   group.by = "patEWS_w_type",
                   cluster.idents = F) + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5))
gg_temp
ggsave(paste0(output_targets, "Mooney_et_al_immunotherapeutic_groups1and2_by_pat_and_celltype.PDF"),
       plot = gg_temp, width = 24, height = 6, units = "in")


################################################################################
####### -------              Run CellChat by sample              ------- #######

# https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html

# Generate cellchat objects for each sample (plotting later)
for(sample in unique(so$Sample)){
  print(sample)
  
  so_sample <- subset(so,
                      Sample == sample)
  
  # Create a CellChat object
  cellchat <- createCellChat(object = so_sample, 
                             group.by = "Sub.Cell.Type")
  
  # Set the ligand-receptor interaction database
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
  showDatabaseCategory(CellChatDB)
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  # set the used database in the object
  cellchat@DB <- CellChatDB.use
  
  # Preprocessing the expression data for cell-cell communication analysis
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  future::plan("multiprocess", workers = 4) 
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # Compute the communication probability and infer cellular communication network
  cellchat <- computeCommunProb(cellchat,
                                population.size = data_is_unsorted)
  
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = min_cells_to_test)
  
  # Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  
  # Calculate the aggregated cell-cell communication network
  cellchat <- aggregateNet(cellchat)
  
  # Save RDS object for future - mainly because this step takes a while
  saveRDS(cellchat,
          paste0(out_cellchat_objects, sample, ".RDS"))
  
}


################################################################################
### ---- Significant interactions and bubble plot with cell of interest ---- ###

df.nets.all <- data.frame()
for(sample in unique(so$Sample)){
  
  cellchat <- readRDS(paste0(out_cellchat_objects, sample, ".RDS"))
  
  # Extract the inferred cellular communication network as a data frame
  df.net <- subsetCommunication(cellchat,
                                thresh = 1) # Default outputs only pval < 0.05
  df.net$sample <- sample
  
  df.net.sig <- df.net[df.net$pval < 0.05,]
  write.table(df.net.sig,
              file = paste0(out_cellchat, sample, "_sig_interactions.txt"),
              sep = "\t", row.names = F, quote = F)
  write.table(df.net,
              file = paste0(out_cellchat, sample, "_all_interactions.txt"),
              sep = "\t", row.names = F, quote = F)
  
  print(paste0(sample, ": ", nrow(df.net), " total interacts, sig : ", nrow(df.net.sig)))
  
  # Cell of interest as a source
  bubble <- netVisual_bubble(cellchat, sources.use = cell_of_interest,
                             remove.isolate = FALSE,
                             thresh = 0.05)
  ggsave(paste0(out_cellchat, sample, "_EWS_as_source_bubble.pdf"),
         bubble)
  
  # Cell of interest as a target
  bubble <- netVisual_bubble(cellchat, targets.use = cell_of_interest,
                             remove.isolate = FALSE,
                             thresh = 0.05)
  ggsave(paste0(out_cellchat, sample, "_EWS_as_target_bubble.pdf"),
         bubble)
  
  # Heatmap showing number of significant interaction by cell type
  hm <- netVisual_heatmap(cellchat, color.heatmap = "Reds")
  pdf(paste0(out_cellchat, sample, "_heatmap_sig_interactions.pdf"))
  print(hm)
  dev.off()
  
  if(nrow(df.nets.all) == 0){
    df.nets.all <- df.net
  }else{
    df.nets.all <- rbind(df.nets.all, df.net)
  }
  
}

df.nets.all_save <- df.nets.all


################################################################################
##### ----- Summarize sig interactions for ALL cell types by patient ----- #####

# In the next section and most of the figures, I am summarizing the data
# to only significant interaction involving the cell type of interest (Ewing)

# Here, I am summarizing all significant interactions THAT DO NOT involve Ewing, 
# so I can focus more on immune-immune interactions

df.nets <- df.nets.all_save[df.nets.all_save$pval < 0.05,]
df.nets$source_target <- paste0(df.nets$source, "_", df.nets$target)
table(df.nets$source_target)
df.nets$source_target_interaction <- paste0(df.nets$source_target, "_",
                                            df.nets$interaction_name_2)
table(df.nets$source_target_interaction) 
df.nets_interest <- df.nets[!(df.nets$source == cell_of_interest | df.nets$target == cell_of_interest),]
df.nets_interest$source_target <- paste0(df.nets_interest$source, "_", df.nets_interest$target)
table(df.nets_interest$source_target)
df.nets_interest$source_target_interaction <- paste0(df.nets_interest$source_target, ".",
                                                     df.nets_interest$interaction_name_2)
interest_interactions <- data.frame(unclass(table(df.nets_interest$source_target_interaction))) 
interest_interactions <- cbind.data.frame(source_target = sapply(strsplit(rownames(interest_interactions), split = "\\."),
                                                                 function(x){x[[1]]}),
                                          interaction = sapply(strsplit(rownames(interest_interactions), split = "\\."),
                                                               function(x){x[[2]]}),
                                          number_patients = interest_interactions[,1])
# Since there is a 10 cell limit to test the interaction, I want to label
# the number of patients where the interaction could be tested (this is by cell type)
df.nets.all_save$source_target <- paste0(df.nets.all_save$source, "_", df.nets.all_save$target)
df.nets.all.sub <- df.nets.all_save[,c("sample", "source_target")]
df.nets.all.sub <- unique(df.nets.all.sub)
possible_interactions <- data.frame(unclass(table(df.nets.all.sub$source_target))) 
possible_interactions <- cbind.data.frame(source_target = rownames(possible_interactions),
                                          testable_samples = possible_interactions[,1])
interest_interactions <- join(interest_interactions,
                              possible_interactions)
# Order
interest_interactions <- interest_interactions[order(interest_interactions$testable_samples, decreasing = F),]
interest_interactions <- interest_interactions[order(interest_interactions$number_patients, decreasing = T),]
# Write sig interactions that contain cell of interest
write.table(interest_interactions,
            file = paste0(out_cellchat, "IMMUNE_interaction_summary.txt"),
            sep = "\t", row.names = F, quote = F)



###### ------ Stacked bar plot of select interactions

out_cellchat_violin <- paste0(out_cellchat, "violin_immune/")
dir.create(out_cellchat_violin)

int_df <- interest_interactions[interest_interactions$number_patients > 2,]
int_df <- rbind(cbind.data.frame(source_target = int_df$source_target,
                                 interaction = int_df$interaction,
                                 Testable = "Nonsignificant samples",
                                 samples = int_df$testable_samples - int_df$number_patients,
                                 mock = paste0("X", 1:nrow(int_df))),
                cbind.data.frame(source_target = int_df$source_target,
                                 interaction = int_df$interaction,
                                 Testable = "Significant samples",
                                 samples = int_df$number_patients,
                                 mock = paste0("X", 1:nrow(int_df)))
)
int_df$mock <- factor(int_df$mock,
                      levels = paste0("X", 1:(nrow(int_df)/2)))
gg_stacked <- ggplot(int_df, aes(fill=Testable, y=samples, x=mock)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_x_discrete(labels = gsub("_", " - ",int_df$source_target[1:(nrow(int_df)/2)])) + coord_flip() +
  theme_classic() + annotate("text", 
                             x = seq(1, (nrow(int_df)/2), length.out =  (nrow(int_df)/2)),
                             y = max(int_df$samples) + 1.25,
                             size = 2.5,
                             label = int_df$interaction[1:(nrow(int_df)/2)],
                             hjust = 0) + 
  expand_limits(y = c(1, 7.2)) + 
  scale_y_continuous(breaks= pretty_breaks()) + xlab("Source - Target") + 
  ylab("Samples (n = 7)") + guides(fill=guide_legend(title="Testable samples"))
gg_stacked
ggsave(filename = paste0(out_cellchat_violin, "immune_stacked_sample_interesting_only.pdf"),
       width = 9.5, height = 3.5)


################################################################################
### --- Summarize sig interaction (containing interesting cell) by patient -- ##

df.nets <- df.nets.all[df.nets.all$pval < 0.05,]
df.nets$source_target <- paste0(df.nets$source, "_", df.nets$target)
table(df.nets$source_target)
df.nets$source_target_interaction <- paste0(df.nets$source_target, "_",
                                            df.nets$interaction_name_2)
table(df.nets$source_target_interaction) 
df.nets_interest <- df.nets[df.nets$source == cell_of_interest | df.nets$target == cell_of_interest,]
df.nets_interest$source_target <- paste0(df.nets_interest$source, "_", df.nets_interest$target)
table(df.nets_interest$source_target)
df.nets_interest$source_target_interaction <- paste0(df.nets_interest$source_target, ".",
                                                df.nets_interest$interaction_name_2)
interest_interactions <- data.frame(unclass(table(df.nets_interest$source_target_interaction))) 
interest_interactions <- cbind.data.frame(source_target = sapply(strsplit(rownames(interest_interactions), split = "\\."),
                                                            function(x){x[[1]]}),
                                     interaction = sapply(strsplit(rownames(interest_interactions), split = "\\."),
                                                          function(x){x[[2]]}),
                                     number_patients = interest_interactions[,1])
# Since there is a 10 cell limit to test the interaction, I want to label
# the number of patients where the interaction could be tested (this is by cell type)
df.nets.all$source_target <- paste0(df.nets.all$source, "_", df.nets.all$target)
df.nets.all.sub <- df.nets.all[,c("sample", "source_target")]
df.nets.all.sub <- unique(df.nets.all.sub)
possible_interactions <- data.frame(unclass(table(df.nets.all.sub$source_target))) 
possible_interactions <- cbind.data.frame(source_target = rownames(possible_interactions),
                                          testable_samples = possible_interactions[,1])
interest_interactions <- join(interest_interactions,
                              possible_interactions)
# Order
interest_interactions <- interest_interactions[order(interest_interactions$testable_samples, decreasing = F),]
interest_interactions <- interest_interactions[order(interest_interactions$number_patients, decreasing = T),]
# Write sig interactions that contain cell of interest
write.table(interest_interactions,
            file = paste0(out_cellchat, "EWS_interaction_summary.txt"),
            sep = "\t", row.names = F, quote = F)


################################################################################
###### ------ Sample # of stacked plot for interesting interactions ------ #####

### Stuff to consider
# Should I show expression in patients with <10 cells? Would be tough to remove in the violin
# Probability is different than p-value
#     In some cases, there are very low probabilities with very sig p-values
#     p-values come from a permutation

interesting_interactions <- read.xlsx(paste0(out_cellchat, 
                                      "EWS_interaction_summary_to_highlight.xlsx"))

# Stacked barplot of highlighted interactions <- mock data
specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)

# Stacked barplot of highlighted interactions
# Goal format: https://r-graph-gallery.com/48-grouped-barplot-with-ggplot2
int_df <- interesting_interactions[n,]
int_df <- rbind(cbind.data.frame(source_target = int_df$source_target,
                           interaction = int_df$interaction,
                           Testable = "Nonsignificant samples",
                           samples = int_df$testable_samples - int_df$number_patients,
                           mock = paste0("X", 1:nrow(int_df))),
          cbind.data.frame(source_target = int_df$source_target,
                 interaction = int_df$interaction,
                 Testable = "Significant samples",
                 samples = int_df$number_patients,
                 mock = paste0("X", 1:nrow(int_df)))
)
int_df$mock <- factor(int_df$mock,
                      levels = paste0("X", 1:(nrow(int_df)/2)))

gg_stacked <- ggplot(int_df, aes(fill=Testable, y=samples, x=mock)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_x_discrete(labels = gsub("_", " - ",int_df$source_target[1:(nrow(int_df)/2)])) + coord_flip() +
  theme_classic() + annotate("text", 
                             x = seq(1, (nrow(int_df)/2), length.out =  (nrow(int_df)/2)),
                             y = max(int_df$samples) + 1.25,
                             size = 2.5,
                             label = int_df$interaction[1:(nrow(int_df)/2)],
                             hjust = 0) + 
  expand_limits(y = c(1, 9.2)) + 
  scale_y_continuous(breaks= pretty_breaks()) + xlab("Source - Target") + 
  ylab("Samples (n = 7)") + guides(fill=guide_legend(title="Testable samples"))
ggsave(filename = paste0(out_cellchat, "stacked_sample_interesting_only.pdf"),
       width = 9.5, height = 3.5)


################################################################################
####### ------- Violin plots for interesting IMMUNE interactions ------- #######

### Stuff to consider
# Should I show expression in patients with <10 cells? Would be tough to remove in the violin
# Probability is different than p-value
#     In some cases, there are very low probabilities with very sig p-values
#     p-values come from a permutation

out_cellchat_violin <- paste0(out_cellchat, "violin_immune/")
dir.create(out_cellchat_violin)

highlighted_interactions_ligand <- c("CD69", "CLEC2D", "MIF")
highlighted_interactions_receptor <- c("KLRB1", "CD74", "CXCR4")

# Combined patients, both ligand and receptor, all cell types
gg_CombPat_AllGene <- VlnPlot(so,
                              features = c(highlighted_interactions_ligand,
                                           highlighted_interactions_receptor),
                              stack = T,
                              flip=T,
                              group.by = "Sub.Cell.Type",
                              #split.by = "Sub.Cell.Type",
                              pt.size = 0.1) + NoLegend() + theme(axis.text.x = element_text(size = 5))

my_color_palette <- hue_pal()(length(unique(so$Sample)))
gg_by_sample_AllGene <- VlnPlot(so,
                      features = c(highlighted_interactions_ligand,
                                   highlighted_interactions_receptor),
                      stack = T,
                      flip=T,
                      group.by = "Sub.Cell.Type",
                      split.by = "Sample",
                      pt.size = 0.1) + 
                theme(axis.text.x = element_text(size = 5)) + 
                geom_vline(xintercept=(1:(length(unique(so$Sub.Cell.Type))-1)+0.5)) + 
                scale_fill_manual(values = my_color_palette)

so_T <- subset(so, Broad.Cell.Type == "T")
my_color_palette2 <- hue_pal()(length(unique(so_T$Sample)))
gg_by_sample_AllGene_T <- VlnPlot(so_T,
                                features = c(highlighted_interactions_ligand,
                                             highlighted_interactions_receptor),
                                stack = T,
                                flip=T,
                                group.by = "Sub.Cell.Type",
                                split.by = "Sample",
                                pt.size = 0.1) + 
  theme(axis.text.x = element_text(size = 5)) + 
  geom_vline(xintercept=(1:(length(unique(so_T$Sub.Cell.Type))-1)+0.5)) + 
  scale_fill_manual(values = my_color_palette2)

# Save plots
gg_CombPat_AllGene
ggsave(filename = paste0(out_cellchat_violin, "All_cell_types_CombinedPatient.pdf"),
       width = 10, height = 5)

gg_by_sample_AllGene
ggsave(filename = paste0(out_cellchat_violin, "All_cell_types_ByPatient.pdf"),
       width = 11, height = 6)

gg_by_sample_AllGene_T
ggsave(filename = paste0(out_cellchat_violin, "T_cells_only_ByPatient.pdf"),
       width = 11, height = 6)

t_dist <- data.frame(unclass(table(so_T$Sample,
                                   so_T$Sub.Cell.Type)))
t_dist <- cbind.data.frame(Sample = rownames(t_dist), t_dist)
write.table(t_dist,
            file = paste0(out_cellchat_violin, "T_cell_distribution_by_sample.txt"),
            sep = "\t", row.names = F, quote = F)


################################################################################
####### -------     Violin plots for general interactions        ------- #######

out_cellchat_violin_general <- paste0(out_cellchat_violin,
                                      "/two_main_interactions_v2/")
dir.create(out_cellchat_violin_general)

general_interactions <- read.xlsx(paste0(out_cellchat, "two_main_interactions_v2.xlsx"))

for(gen_interaction in unique(general_interactions$general_interaction)){
  print(gen_interaction)
  
  VlnPlot(so_vio_target,
          features = general_interactions$Gene[general_interactions$general_interaction == gen_interaction],
          stack = T,
          flip=T,
          group.by = "Sub.Cell.Type",
          split.by = "Sample",
          pt.size = 0.1) + 
    theme(axis.text.x = element_text(size = 5)) + 
    geom_vline(xintercept=(1:(length(unique(so_vio_target$Sub.Cell.Type))-1)+0.5)) + 
    scale_fill_manual(values = my_color_palette)
  
  ggsave(filename = paste0(out_cellchat_violin_general, gen_interaction, ".pdf"),
         width = 8, height = 6)
  
}



