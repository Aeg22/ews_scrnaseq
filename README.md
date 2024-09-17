# Description

This repository contains the main scripts used in the analysis of the Ewing sarcoma scRNA-seq data from the publication *Single cell RNA-sequencing of Ewing sarcoma tumors demonstrates transcriptional heterogeneity and clonal evolution* (PMID: XXXX).

# Code overview

Code is categorized and described in this README. 

### Cell Ranger

Scripts utilized during Cell Ranger alignment and counting steps. Aggregation was performed without normalization. The main outputs of the aggregation step are deposited to GEO as GSE277083.

Shell scripts:

* `cellranger/cellranger-count-example.sh`
* `cellranger/cellranger-aggregation.sh`

### Seurat initialization, quality control, and filtering

Generation of initial Seurat object, exploration of data quality, and filtering genes, cells, and samples. 

R scripts:

* `r_scripts/1-aggr_to_seurat.R`
* `r_scripts/1-filter.R`

### Analysis of primary tumors

The cells of the primary tumors were extracted and analyzed separately from the peripheral blood. These scripts perform cell clustering, cell type naming, cell proportions, pathway analysis, copy number alterations (CNAs), integration, cNMF, and identification and characterization of CNA subclusters. 

R scripts:

* `r_scripts/2-norm_cluster_celltype_Primary.R`
* `r_scripts/2-scevan_primary.R`
* `r_scripts/2-inferCNV_primary.R`
* `r_scripts/2-integration_and_cNMF.R`
* `r_scripts/2-inferCNV_primary_subclusters.R`
* `r_scripts/2-CNAsub_heterogeneity.R`

### Subclustering immune cells in primary tumors

The T and myeloid subsets were subclustered and subtyped. 

R scripts:

* `r_scripts/3-primary_subcluster_T.R`
* `r_scripts/3-primary_subcluster_myeloid.R`

### Additional analysis of primary tumors

Primary tumors were further analyzed for expression of potential EWS targets and cell-cell communication.

R scripts:

* `r_scripts/4-ewstargets_and_cellchat.R`

### Analysis of CTC/peripheral blood samples

The cells of the peripheral blood were extracted and analyzed separately. These scripts cluster cells, name cell types, calculate cell proportions, and run pathway analysis and inferCNV.

R scripts:

* `r_scripts/5-norm_cluster_celltype_CTC.R`

### Data availability

Scripts related to generation of supplementary tables and processing data for GEO submission. 

R scripts:

* `r_scripts/6-supp_tables.R`
* `r_scripts/6-GEO_processed.R`

# Releases and changes

* 1.0
  * Initial Release
