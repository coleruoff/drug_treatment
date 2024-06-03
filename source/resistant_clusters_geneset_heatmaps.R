dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/"
# dataDirectory <- "//hpcdrive.nih.gov/CDSL_hannenhalli/Cole/projects/drug_treatment/"

geneset_group_matrix <- function(data, ident_to_use, genesets_to_use)

curr_cell_line <- "A549"

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

RACs <- list(c(4,9,12,13,14,16,18),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- RACs[[curr_cell_line]]


scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
clusters_of_interest <- RACs[[curr_cell_line]]

data <- AddMetaData(data, metadata = ifelse(colnames(data) %in% active_cell_names, "active","inactive"), col.name = "resistant_active")
data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & data$resistant_active == "active", "resistant","nonresistant"), col.name = "resistant")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & data$resistant_active == "active", paste0(data$Cluster, "_resistant"),data$Cluster), col.name = "resistant_cluster")


genesets_name <- "ITH_meta_programs"
ret_mat <- geneset_group_matrix(data, "resistant_cluster", genesets_name)

#Set all non-significant heatmap cells to 0
heatmap_matrix <- ifelse(ret_mat[[2]]<.05,ret_mat[[1]], 0)

colnames(heatmap_matrix) <- paste0("Cluster ",colnames(ret_mat[[1]]))
rownames(heatmap_matrix) <- rownames(ret_mat[[1]])


RACs <- list(c(4,9,12,13,14,16,18),c(4,5,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- paste0("Cluster ",RACs[[curr_cell_line]],"_resistant")

rac_ha <- HeatmapAnnotation(RAC = c(ifelse(colnames(heatmap_matrix) %in% clusters_of_interest,"RAC","Non-RAC")),
                            col = list(RAC = c("RAC" = "darkgreen", "Non-RAC" = "lightblue")))


ht <- Heatmap(heatmap_matrix, name="Z-Score", cluster_rows = T, cluster_columns = T,
              bottom_annotation = rac_ha, column_title = "Clusters", column_title_side = "bottom",
              row_title = "Cancer Meta-Programs", column_names_rot = 45)


draw(ht, column_title = paste0(str_to_title(genesets_name)," Mean AUCell Score in Each Cluster (", curr_cell_line, ")\n"), 
     column_title_gp = gpar(fontsize = 26),  padding = unit(c(2, 2, 2, 60), "mm"),
     heatmap_legend_side = "left", annotation_legend_side = "left")





