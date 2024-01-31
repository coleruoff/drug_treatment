library(ComplexHeatmap)
library(Seurat)
library(circlize)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/"
# dataDirectory <- "//hpcdrive.nih.gov/CDSL_hannenhalli/Cole/projects/drug_treatment/"

geneset_group_matrix <- function(data, ident_to_use, genesets_to_use){
  
  Idents(data) <- ident_to_use
  
  num_groups <- nlevels(data)
  
  groups <- levels(data)
  
  cat(genesets_to_use, "\n")
  
  scores <- readRDS(paste0(dataDirectory, "data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", genesets_to_use, "_aucell_scores.rds"))
  threshold <- readRDS(paste0(dataDirectory, "data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", genesets_to_use, "_aucell_thresholds.rds"))
  
  # scores <- scale(scores)
  
  geneset_names <- colnames(scores)
  
  score_matrix <- matrix(NA, ncol=length(groups), nrow=length(geneset_names))
  pvalue_matrix <- matrix(NA, ncol=length(groups), nrow=length(geneset_names))
  
  for(curr_group in groups){
    
    cat(curr_group, "\n")
    
    curr_cell_names <- colnames(data)[data[[ident_to_use]] == curr_group]
    rest_names <- colnames(data)[data[[ident_to_use]] != curr_group]
    
    for(curr_geneset in geneset_names){
  
      curr_scores <- scores[curr_cell_names, curr_geneset]
           
      rest_scores <- scores[rest_names, curr_geneset]
      
      wilcox_res <- wilcox.test(curr_scores,rest_scores)
      
      i <- which(curr_geneset == geneset_names)
      j <- which(curr_group == groups)
      
      score_matrix[i,j] <- mean(curr_scores)
      pvalue_matrix[i,j] <- wilcox_res$p.value
      
    }
  }
  
  scaled_matrix <- t(scale(t(score_matrix)))
  
  colnames(scaled_matrix) <- groups
  rownames(scaled_matrix) <- geneset_names
  
  colnames(score_matrix) <- groups
  rownames(score_matrix) <- geneset_names
  
  return(list(scaled_matrix, pvalue_matrix,score_matrix))
}

################################################################################

curr_cell_line <- "A549"

scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- RACs[[curr_cell_line]]

data <- readRDS(paste0(dataDirectory, "data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
data <- AddMetaData(data, metadata = ifelse(colnames(data) %in% active_cell_names, "active","inactive"), col.name = "resistant_active")
data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")

data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & data$resistant_active == "active", "resistant","nonresistant"), col.name = "resistant")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & data$resistant_active == "active", paste0(data$Cluster, "_resistant"),data$Cluster), col.name = "resistant_cluster")


genesets_name <- "hallmarks"

ret_mat <- geneset_group_matrix(data, "resistant", genesets_name)

#Set all non-significant heatmap cells to 0
heatmap_matrix <- ifelse(ret_mat[[2]]<.05,ret_mat[[1]], 0)

colnames(heatmap_matrix) <- paste0("Cluster ",colnames(ret_mat[[1]]))
rownames(heatmap_matrix) <- rownames(ret_mat[[1]])


RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- paste0("Cluster ",RACs[[curr_cell_line]])

rac_ha <- HeatmapAnnotation(RAC = c(ifelse(colnames(heatmap_matrix) %in% clusters_of_interest,"RAC","Non-RAC")),
                                col = list(RAC = c("RAC" = "orange", "Non-RAC" = "lightblue")))


ht <- Heatmap(heatmap_matrix, name="Z-Score", cluster_rows = F, cluster_columns = T,
              bottom_annotation = rac_ha, column_title = "Clusters", column_title_side = "bottom",
              row_title = "", column_names_rot = 45)


draw(ht, column_title = paste0("Cancer Meta-Programs Mean AUCell Score in Each Cluster (", curr_cell_line, ")\n"), 
     column_title_gp = gpar(fontsize = 26),  padding = unit(c(2, 2, 2, 60), "mm"),
     heatmap_legend_side = "left", annotation_legend_side = "left")

######
# active vs inactive heatmap

heatmap_matrix <- ifelse(ret_mat[[2]]<.05,ret_mat[[3]], 0)

colnames(heatmap_matrix) <- paste0("Cluster ",colnames(ret_mat[[1]]))
rownames(heatmap_matrix) <- rownames(ret_mat[[1]])

col_fun <-colorRamp2(c(0, median(heatmap_matrix), max(heatmap_matrix)), c("blue","white", "red"))
Heatmap(heatmap_matrix, cluster_rows = F, col = col_fun)

summary(heatmap_matrix)

###################

rac_values <- as.vector(heatmap_matrix[,clusters_of_interest])

rest_values <- as.vector(heatmap_matrix[,!colnames(heatmap_matrix) %in% clusters_of_interest])

boxplot(rac_values,rest_values, names=c("RACs","Rest"), main=paste0(curr_cell_line, " ", genesets_name, " AUCell Scores"))
wilcox.test(rac_values,rest_values)

wilcox_res <- wilcox.test(rac_values,rest_values)
boxplot(rac_values,rest_values, names=c("RACs","Rest"), main=paste0(genesets_name, " AUCell Scores (", curr_cell_line, ") pval: ", sprintf("%.4f", wilcox_res$p.value)))



heatmap_matrix_up <- heatmap_matrix[!grepl("down",rownames(heatmap_matrix)),]


ht <- Heatmap(heatmap_matrix_up, name="Z-Score", cluster_rows = F, cluster_columns = T,
              bottom_annotation = rac_ha, column_title = "Clusters", column_title_side = "bottom",
              row_title = "Cancer Meta-Programs", column_names_rot = 45)


draw(ht, column_title = paste0(genesets_name, " Up-Genes Mean AUCell Score in Each Cluster (", curr_cell_line, ")\n"), 
     column_title_gp = gpar(fontsize = 26),  padding = unit(c(2, 2, 2, 60), "mm"),
     heatmap_legend_side = "left", annotation_legend_side = "left")



rac_values <- as.vector(heatmap_matrix_up[,clusters_of_interest])

rest_values <- as.vector(heatmap_matrix_up[,!colnames(heatmap_matrix_up) %in% clusters_of_interest])


boxplot(rac_values,rest_values, names=c("RACs","Rest"), main=paste0(curr_cell_line, " ", genesets_name, " Upregulated Genes AUCell Scores"))
wilcox.test(rac_values,rest_values)


heatmap_matrix_down <- heatmap_matrix[!grepl("up",rownames(heatmap_matrix)),]


ht <- Heatmap(heatmap_matrix_down, name="Z-Score", cluster_rows = F, cluster_columns = T,
              bottom_annotation = rac_ha, column_title = "Clusters", column_title_side = "bottom",
              row_title = "Cancer Meta-Programs", column_names_rot = 45)


draw(ht, column_title = paste0(genesets_name, " Up-Genes Mean AUCell Score in Each Cluster (", curr_cell_line, ")\n"), 
     column_title_gp = gpar(fontsize = 26),  padding = unit(c(2, 2, 2, 60), "mm"),
     heatmap_legend_side = "left", annotation_legend_side = "left")



rac_values <- as.vector(heatmap_matrix_down[,clusters_of_interest])

rest_values <- as.vector(heatmap_matrix_down[,!colnames(heatmap_matrix_down) %in% clusters_of_interest])


boxplot(rac_values,rest_values, names=c("RACs","Rest"), main=paste0(curr_cell_line, " ", genesets_name, " Downregulated Genes AUCell Scores"))
wilcox.test(rac_values,rest_values)



temp <- as.matrix(dist(t(heatmap_matrix)))

Heatmap(temp,bottom_annotation = rac_ha)


Heatmap(cor(heatmap_matrix),bottom_annotation = rac_ha)




score_matrix <- ret_mat[[3]]
colnames(score_matrix) <- paste0("Cluster ",colnames(ret_mat[[1]]))
rownames(score_matrix) <- rownames(ret_mat[[1]])

Heatmap(score_matrix)

temp <- as.matrix(dist(t(score_matrix)))

Heatmap(temp,bottom_annotation = rac_ha, column_title = paste0(curr_cell_line," ", genesets_name, " distances"), name = "distance")
Heatmap(cor(score_matrix,method = "spearman"),bottom_annotation = rac_ha,column_title = paste0(curr_cell_line," ", genesets_name, " correlation"), name = "spearman\ncorrelation")


plot(density(score_matrix[,9]))

cor(heatmap_matrix[,9],heatmap_matrix[,13])

dist(heatmap_matrix[,9],heatmap_matrix[,13])

