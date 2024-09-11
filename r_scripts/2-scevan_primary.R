
################################################################################
####### -------                    Goals                         ------- #######

# Processing the primary samples only with SCEVAN

################################################################################
####### -------                 Load R packages                  ------- #######

library(Seurat)
library(ggplot2)
library(devtools)
library(SCEVAN)
library(gridExtra)
library(ggthemes)
library(extrafont)
library(remotes)
extrafont::font_import()

source("functions/SCEVAN/custom_multisample.R")
# I edited the SCEVAN function to skip the plotOncoHeatSubclones() plotting step
# If there aren't subclones called then the plot fails and the function would
# fail before the results data frame is written
source("functions/SCEVAN/custom_pipelineCNA.R")


################################################################################
####### -------                Update variables                  ------- #######

project <- "/Volumes/Partition 2/Core/Hayashi_Ewing_scRNAseq_Dec2022/"
setwd(project)

tumor_cells <- c("Ewing sarcoma", "Ewing sarcoma proliferating") # Used only for clone step
input_rds = "Data/Seurat_objects/AggrNoNorm_filtered_primary_celltype.RDS"

output_all_sample <- "Analysis/Annotation/Primary/Normalized/SCEVAN/multiSample/"
dir.create(output_all_sample, recursive = T)

output_per_sample <- "Analysis/Annotation/Primary/Normalized/SCEVAN/perSample_beta_2/"
dir.create(output_per_sample, recursive = T)


################################################################################
####### -------                 Read in data                     ------- #######

so <- readRDS(input_rds)

# Make list of matrices (one per sample)
listCountMtx <- list(Primary_001 = subset(so, Sample == "001_Primary")@assays$RNA@counts,
                     Primary_005 = subset(so, Sample == "005_Primary")@assays$RNA@counts,
                     Primary_010 = subset(so, Sample == "010_Primary")@assays$RNA@counts,
                     Primary_038 = subset(so, Sample == "038_Primary")@assays$RNA@counts,
                     Primary_051 = subset(so, Sample == "051_Primary")@assays$RNA@counts,
                     Primary_061 = subset(so, Sample == "061_Primary")@assays$RNA@counts,
                     Primary_066 = subset(so, Sample == "066_Primary")@assays$RNA@counts)


################################################################################
####### -------             Run SCEVAN via multisample           ------- #######

setwd(output_all_sample)
results_multiSample <- multiSampleComparisonClonalCN_AGreorder(listCountMtx, 
                                                               analysisName = "All_primary", 
                                                               organism = "human", 
                                                               par_cores = 10)

dev.off()

load("output/All_primary_outputAnalysis.RData")
calls_df_list <- outputAnalysis[[1]]
cell_calls <- bind_rows(calls_df_list)
cell_calls$confidentNormal <- gsub("yes", "confident ", cell_calls$confidentNormal)
cell_calls$class.detail <- paste0(cell_calls$confidentNormal,
                                  cell_calls$class)
cell_calls$class.detail <- gsub("NA", "", cell_calls$class.detail)
cell_calls <- cell_calls[,-2]

# Add SCEVAN tumor/normal calls to seurat object
sum(!(rownames(cell_calls) %in% colnames(so))) # All are present
cell_calls <- cell_calls[match(colnames(so), rownames(cell_calls)),]
so$SCEVAN_class <- cell_calls$class
so$SCEVAN_class.detail <- cell_calls$class.detail


################################################################################
####### -------               Plot SCEVAN calls                  ------- #######

# Write
so_meta <- so@meta.data
so_meta <- cbind.data.frame(barcodes = rownames(so_meta),
                            so_meta)
write.table(so_meta,
            file = "SCEVAN_metadata.pdf",
            sep = "\t", row.names = F)

###### First, at the dataset level

prop_df_sample <- data.frame(prop.table(table(so@meta.data[,"SCEVAN_class"], 
                                              so@meta.data[,"Cell.Type"]), margin = 2))
gg_stackbar_dataset <- ggplot(prop_df_sample, aes(fill=Var1, y=Freq, x=Var2)) + 
  geom_bar(position="stack", stat="identity") +ylab("Proportion") + xlab("Cell type") +
  theme_classic() + 
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.5, size = 10)) + 
  labs(fill='SCEVAN call') 
gg_stackbar_dataset
ggsave(paste0("CellType_SCEVAN_stacked_barplot.pdf"),
        gg_stackbar_dataset)
