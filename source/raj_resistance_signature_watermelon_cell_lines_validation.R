library(Seurat)
source("source/cole_functions.R")

watermelon_cell_lines <- list.files("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/watermelon_data")
watermelon_cell_lines <- watermelon_cell_lines[1:8]

#remove this cell line since theres no day 0 cells
watermelon_cell_lines <- watermelon_cell_lines[-5]

resistance_heatmap <- matrix(NA, nrow=length(watermelon_cell_lines), ncol = 2)

curr_cell_line <- watermelon_cell_lines[1]

for(curr_cell_line in watermelon_cell_lines){
  
  watermelon_data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/watermelon_data/",curr_cell_line))
  
  if("cell_line" %in% colnames(watermelon_data@meta.data)){
    watermelon_data <- AddMetaData(watermelon_data,metadata = ifelse(watermelon_data$sample_type == "naive", 0, 10), col.name = "time_point")
  }
  
  ################################################################################
  # Calculate FC for each time point
  ################################################################################
  
  Idents(watermelon_data) <- watermelon_data$time_point
  
  time_points <- levels(Idents(watermelon_data))
  
  time_points_fc <- list()
  for(curr_time_point in time_points){
    
    curr_fc_result <- FoldChange(watermelon_data, ident.1 = curr_time_point)
    
    time_points_fc <- append(time_points_fc, list(curr_fc_result))
    
  }
  
  
  names(time_points_fc) <- paste0("pc9_day",time_points,"_fc")
  
  ################################################################################
  # Run GSEA of raj resistance signature for each time point
  ################################################################################
  
  #Create ranked genelists for each time point
  genesets_with_ranks <- list()
  for(curr_fc_res in time_points_fc){
    
    ranks <- curr_fc_res$avg_log2FC
    
    names(ranks) <- rownames(curr_fc_res)
    
    ranks <- sort(ranks, decreasing = T)
    
    genesets_with_ranks <- append(genesets_with_ranks, list(ranks))
    
  }
  
  names(genesets_with_ranks) <- paste0("pc9_day",time_points,"_fc")
  
  #Read in resistance signature
  raj_resistant_cancer_type <- "common"
  
  # refined_raj_markers <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_", raj_resistant_cancer_type, "_resistance_refined_cluster_markers.rds"))
  
  resistant_vs_control_de_genes <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_", raj_resistant_cancer_type, "_resistant_vs_control_de_signature_500.rds"))
  resistant_vs_control_de_genes <- list("resistant_vs_control_de_genes" = resistant_vs_control_de_genes)
  
  genesets2 <- resistant_vs_control_de_genes
  
  # genesets2 <- list("raj_resistance_signature" = raj_resistance_signature)
  
  # Run GSEA
  
  result <- create_GSEA_matrix(genesets_with_ranks, genesets2)
  
  
  gsea_matrix <- result[[1]]
  
  i <- which(curr_cell_line == watermelon_cell_lines)
  
  resistance_heatmap[i,] <- gsea_matrix
}
 


Heatmap(gsea_matrix, cluster_columns = F, cluster_rows = F, column_title = plot_title, name="NES", 
        column_names_rot = 45,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.5f", gsea_matrix[i, j]), x, y, gp = gpar(fontsize = 20))})



colnames(resistance_heatmap) <- paste0("Day ", time_points)
rownames(resistance_heatmap) <- sapply(watermelon_cell_lines, FUN = function(x) gsub("_processed.rds","",x))

plot_title <- paste0("Raj Common Resistance Signature Enrichment Along\nResistance Development in Multiple Cell Lines")

# resistance_heatmap <- resistance_heatmap[-2,]

Heatmap(resistance_heatmap, cluster_columns = F, column_title = plot_title, name="NES", 
        column_names_rot = 45, column_title_gp = gpar(fontsize=20),column_names_gp = gpar(fontsize=16),
        row_names_gp = gpar(fontsize=14),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.5f", resistance_heatmap[i, j]), x, y, gp = gpar(fontsize = 15))})

