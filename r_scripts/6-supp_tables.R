

################################################################################
####### -------                     Set paths                    ------- #######

setwd("/Volumes/Partition 2/Core/Hayashi_Ewing_scRNAseq_Dec2022/")
output <- "Analysis/Copied_figures/Figures/"
dir.create(output)

library(plyr)
library(openxlsx)

################################################################################
####### -------                   Read in data                   ------- #######

seq_qc <- read.table(file = "Analysis/postFilter_QC/sample-level_metrics.txt",
                     sep = "\t", header = T, stringsAsFactors = F,
                     check.names = F)
seq_qc <- seq_qc[!(seq_qc$Sample %in% c("051_CTC", "066_CTC")),] # Removing samples we removed from the study due to poor QC
seq_qc <- seq_qc[order(seq_qc$Sample),]
seq_qc <- rbind(seq_qc[grep("Primary", seq_qc$Sample),],
                seq_qc[grep("CTC", seq_qc$Sample),])

primary_dge <- read.table(file = "Analysis/Annotation/Primary/Normalized/ORA_cell_type/Cell_type_FindAllMarkers_sigOnly.txt",
                     sep = "\t", header = T, stringsAsFactors = F)
primary_dge <- primary_dge[,c(6, 7, 1, 5, 2, 3, 4)]
colnames(primary_dge) <- c("Cell type", "Gene", "P value", "Adjusted p value", 
                           "Log2 fold-change", "% cells in Cell type expressing Gene",
                           "% cells not in Cell type expressing Gene")

t_dge <- read.table(file = "Analysis/Annotation/Primary/subcluster_Tcells/ORA_cell_type/Cell_type_FindAllMarkers_sigOnly.txt",
                     sep = "\t", header = T, stringsAsFactors = F)
t_dge <- t_dge[,c(6, 7, 1, 5, 2, 3, 4)]
colnames(t_dge) <- c("Cell type", "Gene", "P value", "Adjusted p value", 
                           "Log2 fold-change", "% cells in Cell type expressing Gene",
                           "% cells not in Cell type expressing Gene")

myeloid_dge <- read.table(file = "Analysis/Annotation/Primary/subcluster_myeloid/ORA_cell_type/Cell_type_FindAllMarkers_sigOnly.txt",
                    sep = "\t", header = T, stringsAsFactors = F)
myeloid_dge <- myeloid_dge[,c(6, 7, 1, 5, 2, 3, 4)]
colnames(myeloid_dge) <- c("Cell type", "Gene", "P value", "Adjusted p value", 
                           "Log2 fold-change", "% cells in Cell type expressing Gene",
                           "% cells not in Cell type expressing Gene")

blood_dge <- read.table(file = "Analysis/Annotation/CTC/Normalized/ORA_cell_type/Cell_type_FindAllMarkers_sigOnly.txt",
                        sep = "\t", header = T, stringsAsFactors = F)
blood_dge <- blood_dge[,c(6, 7, 1, 5, 2, 3, 4)]
colnames(blood_dge) <- c("Cell type", "Gene", "P value", "Adjusted p value", 
                           "Log2 fold-change", "% cells in Cell type expressing Gene",
                           "% cells not in Cell type expressing Gene")

cellchat <- read.table(file = "Analysis/Annotation/Primary/Normalized/cellchat/EWS_interaction_summary.txt",
                     sep = "\t", header = T, stringsAsFactors = F)
colnames(cellchat) <- c("Source - Target", "Ligand - Receptor", 
                        "Samples with significance", "Testable samples")

sample_gep <- read.table(file = "Analysis/Annotation/Primary/integrate_P_RPCA_CPM/cNMF_v3/cNMF_top_correlated_genes_with_each_tumor_program.txt",
                         sep = "\t", header = T, stringsAsFactors = F)
merged_gep <- read.table(file = "Analysis/Annotation/Primary/integrate_P_RPCA_CPM/cNMF_v3/merged_programs/Merged_gene_programs.txt",
                     sep = "\t", header = T, stringsAsFactors = F)
colnames(merged_gep) <- c("Program", "Genes")


# Exploration
p <- read.table(file = "/Volumes/Partition 2/Core/Hayashi_Ewing_scRNAseq_Dec2022/Analysis/GEO/cell_metadata_primary.txt",
                     sep = "\t", header = T, stringsAsFactors = F)
b <- read.table(file = "/Volumes/Partition 2/Core/Hayashi_Ewing_scRNAseq_Dec2022/Analysis/GEO/cell_metadata_blood.txt",
                sep = "\t", header = T, stringsAsFactors = F)

################################################################################
####### -------                Prepare excel object              ------- #######

key <- cbind.data.frame(Table = c("Supplementary Table 1",
                                  "Supplementary Table 2",
                                  "Supplementary Table 3",
                                  "Supplementary Table 4",
                                  "Supplementary Table 5",
                                  "Supplementary Table 6",
                                  "Supplementary Table 7",
                                  "Supplementary Table 8"),
  Description = c("Sequencing information of each single-cell RNA-seq sample.",
                  "Significant gene enrichment by broad cell type in primary tumors.",
                  "Genes of each sample-level gene expression program.",
                  "Genes of each merged gene expression program.",
                  "Significant gene enrichment by T cell subtype in primary tumors.",
                  "Significant gene enrichment by myeloid cell subtype in primary tumors.",
                  "Significant cell communication interactions involving Ewing Sarcoma cells in primary tumors.",
                  "Significant gene enrichment by cell type in blood samples."
) )

excel_list <- list(
  key,
  seq_qc,
  primary_dge,
  sample_gep,
  merged_gep,
  t_dge,
  myeloid_dge,
  cellchat,
  blood_dge
)
names(excel_list) <- c("Key",
                       key$Table)


################################################################################
####### -------                Write excel file                  ------- #######

write.xlsx(excel_list, paste0(output, "Supplementary_tables.xlsx"),
           firstActiveRow = T, colWidths = "auto", overwrite = T)




