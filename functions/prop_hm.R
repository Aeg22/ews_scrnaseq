#### Calculate relative proportions across two categories (ex. cell proportions across sample)

prop_hm <- function(...,
                    so = NULL,
                    group1 = NULL,
                    group2 = NULL,
                    order_by = TRUE,
                    outdir = NULL,
                    # Changed several pheatmap() defaults
                    fontsize = 12,
                    #color = colorRampPalette(c("white", "orange", "darkorange", "red", "firebrick3"))(100),
                    color = c("white", colorRampPalette(c("#F3F2FE", "#0E0AFA"))(100)),
                    decimals = 3,
                    # Size of heatmap
                    width = 4,
                    height = 2
                    ){
  
  # This function will output a table and heatmap where `group1` represents the columns
  # and each column's total proportion will be normalized to 1 (ex. Sample). 
  # `group2` are the groups `group1` are divided into rows (ex. cell types). 
  #
  # so = NULL              # Seurat object. Prefilter if necessary
  # group1 = NULL          # Column name for `group1` in `so`. Represents the columns, which are normalized to 1. 
  # group2 = NULL          # Column name for `group2` in `so`. Represents the rows. 
  # order_by = T           # Order rows by descending mean proportions (TRUE). FALSE will order by ascending.  
  #                           Use NULL to keep rows and columns are ordered by the factors in `so`
  # outdir = NULL          # Directory to save files
  # fontsize = 12           # Heatmap fontsize
  # color = c("white", colorRampPalette(c("#F3F2FE", "#0E0AFA"))(100))
  #                        # Heatmap color scheme
  # decimals = 3           # Number of decimals to display on heatmap
  # width = 4              # Heatmap pdf width
  # height = 2             # Heatmap pdf height
  # ...additional arguments for pheatmap()
  #
  # Output
  # 
  # Proportions table
  # Files to outdir
  #     proportion_heatmap.pdf  # pdf heatmap of cell proportions
  #     proportion_table.txt    # txt of cell proportions

  library(Seurat)
  library(dplyr)
  library(pheatmap)
  library(grid)
  
  if (is.null(so) | is.null(group1)  | is.null(group2)  | is.null(outdir) ){
    stop("The following parameters are required: so, group1, group2, outdir")
  }
  
  # Function to save heatmap
  save_pheatmap_pdf <- function(x, filename, width=width, height=height) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  # Prop table
  pt <- data.frame(unclass(prop.table(table(so@meta.data[,group2], 
                                            so@meta.data[,group1]), margin = 2)),
                   check.names = F)
  
  # Reorder rows by mean (skip if order_by = NULL)
  if(!is.null(order_by)){
    row_mean <- rowMeans(pt)
    pt <- pt[order(row_mean, decreasing = order_by),]
  }
  
  # Plot and save heatmap
  # To set max to 1 regardless of data, use the param: breaks = seq(0, 1, by = 0.01)
  ph <- pheatmap(pt,
           fontsize = fontsize,
           color = color,
           display_numbers = T,
           number_format = paste0("%.", decimals,"f"),
           cluster_rows = F,
           cluster_cols = F,
           number_color = "black",
           main = paste0("Proportion of ", group2, " in ", group1),
           ...)
  save_pheatmap_pdf(ph,
                    paste0(outdir, "/proportion_heatmap.pdf"),
                    width=width,height=height)
  
  # Write table
  pt_save <- pt
  pt_save <- cbind.data.frame(group2 = rownames(pt_save), pt_save)
  colnames(pt_save)[1] <- group2
  write.table(pt_save, file = paste0(outdir, "/proportion_table.txt"), 
                sep = "\t", row.names = F, quote = F)
  
  # Return table
  return(pt_save)
  
}
