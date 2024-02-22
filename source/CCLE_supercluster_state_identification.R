library(Seurat)
library(AUCell)
source("/data/ruoffcj/projects/aucell_scoring/aucell_thresholding.R")

calculate_distances <- function(obj,names_x,names_y=NULL){
  
  if(is.null(names_y)){
    
    pc_mat <- obj@reductions$pca@cell.embeddings[names_x,1:20]
    
    
    distances <- stats::dist(pc_mat, diag = F)
    
  } else {
    pc_mat <- obj@reductions$pca@cell.embeddings[c(names_x,names_y),1:20]
    
    
    distances <- as.matrix(stats::dist(pc_mat, diag = F))
    
    distances <- as.numeric(distances[names_x,names_y])
    
  }
  
  return(distances)
}

#################################################################################

supercluster_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/type1_supercluster_signatures.rds")

CCLE_data <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/CCLE_data/CCLE_Object_normalized.rds")
cell_lines <- unique(CCLE_data@meta.data$Cell_line)

# cell_lines <- cell_lines[1:5]
curr_cell_line <- cell_lines[2]

final_df <- list()
for(curr_cell_line in cell_lines){
  cat(curr_cell_line, "\n")
  
  # Create and process Seurat object
  curr_obj <- CreateSeuratObject(counts = CCLE_data@assays$RNA@counts[,CCLE_data$Cell_line == curr_cell_line],
                                 min.cells = 3, 
                                 min.features = 200)
  
  curr_obj <- NormalizeData(curr_obj)
  
  curr_obj <- FindVariableFeatures(curr_obj, selection.method = "vst", nfeatures = 2000)
  
  curr_obj <- ScaleData(curr_obj, features = rownames(curr_obj))
  
  curr_obj <- RunPCA(curr_obj, features = VariableFeatures(object = curr_obj))

  # Score all cells and identify active cells
  # auc_obj <- compute_AUCell_scores(curr_obj, supercluster_signatures, compute_thresholds=F, nCores = 2, assay_to_use = "RNA")
  # 
  # #Compute separate thresholds across time-points
  # computed_thresholds_df <- compute_shuffled_gene_set_AUCell_scores(curr_obj, gene_sets=supercluster_signatures, nCores=2, do_sample_wise=F, q_thresh=0.95, num_controls=100, assay_to_use = "RNA")
  # 
  
  cells_AUC <- AUCell_run(curr_obj@assays$RNA$data, supercluster_signatures)

  cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, assign=TRUE)
  
  # For each signature
  final_df[["cell_line"]] <- append(final_df[["cell_line"]],curr_cell_line)
  for(curr_signature in names(supercluster_signatures)){
    cat(curr_signature,"\n")
    
    # active_cells <- rownames(auc_obj$auc_mat)[auc_obj$auc_mat[,curr_signature] > computed_thresholds_df$threshold[which(curr_signature == names(supercluster_signatures))]]
    # active_cells <- cells_assignment[[curr_signature]]$assignment
     
    top_75 <- quantile(cells_AUC@assays@data$AUC[curr_signature,], probs = .75)
    active_cells <- colnames(cells_AUC)[cells_AUC@assays@data$AUC[curr_signature,]>top_75]
    
    if(length(active_cells) == 0){
      active_cells <- colnames(cells_AUC)[cells_AUC@assays@data$AUC[curr_signature,] > cells_assignment[[curr_signature]]$aucThr$thresholds[2]]
    }
    
    if(length(active_cells) > 4){
      inactive_cells <- colnames(curr_obj)[!colnames(curr_obj) %in% active_cells]
      
      # Calculate distances between active cells and rest
      active_distances <- calculate_distances(curr_obj,active_cells)
      
      rest_distances <- calculate_distances(curr_obj,active_cells, inactive_cells)
      
            # boxplot(active_distances,rest_distances)
      
      # Test distributions for significance
      
      wilcox_res <- wilcox.test(active_distances,rest_distances)
      
      # Declare whether cells are one state or are stochastic
      if(wilcox_res$p.value < 0.05 & median(active_distances) < median(rest_distances)){
        final_df[[curr_signature]] <- append(final_df[[curr_signature]], 1)
        final_df[[paste0(curr_signature,"_active")]] <- append(final_df[[paste0(curr_signature,"_active")]] , length(active_cells))
        final_df[[paste0(curr_signature,"_inactive")]] <- append(final_df[[paste0(curr_signature,"_inactive")]] , length(inactive_cells))
        
      } else{
        final_df[[curr_signature]] <- append(final_df[[curr_signature]], 0)
        final_df[[paste0(curr_signature,"_active")]] <- append(final_df[[paste0(curr_signature,"_active")]] , length(active_cells))
        final_df[[paste0(curr_signature,"_inactive")]] <- append(final_df[[paste0(curr_signature,"_inactive")]] , length(inactive_cells))
      }
    } else{
      final_df[[curr_signature]] <- append(final_df[[curr_signature]], NA)
      final_df[[paste0(curr_signature,"_active")]] <- append(final_df[[paste0(curr_signature,"_active")]] , length(active_cells))
      final_df[[paste0(curr_signature,"_inactive")]] <- append(final_df[[paste0(curr_signature,"_inactive")]] , length(inactive_cells))
    }
  }
}



ifelse(colnames(curr_obj) %in% active_cells, 1, 0)

curr_obj <- AddMetaData(curr_obj, ifelse(colnames(curr_obj) %in% active_cells, 1, 0), col.name = "active")

DimPlot(curr_obj, group.by = "active")


final_df <- do.call(cbind, final_df)



saveRDS(final_df, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/CCLE_supercluster_state_results_75.rds")

final_df <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/CCLE_supercluster_state_results_75.rds")


sum(final_df[,2]==1)/nrow(final_df)

sum(final_df[,5]==1)/nrow(final_df)

# 
# 
# supercluster1_res <- final_df[,2][!is.na(final_df[,2])]
# supercluster2_res <- final_df[,3][!is.na(final_df[,3])]
# 
# sum(supercluster1_res == 1)/length(supercluster1_res)
# sum(supercluster2_res == 1)/length(supercluster2_res)
# 

