library(tidyverse)
library(Seurat)

curr_cell_line <- "A549"

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

RACs <- list(c(4,9,12,13,14,16,18),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- RACs[[curr_cell_line]]


scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
clusters_of_interest <- RACs[[curr_cell_line]]


data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, data$Cluster, "Non-RAC"), col.name = "RAC")
#Add metadata for RAC and Cell Group
data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, "1", ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, paste0(data$Cluster, "_1"), ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), paste0(data$Cluster, "_2"), paste0(data$Cluster, "_0"))), col.name = "cell_cluster_group")

#Supercluster scores and thresholds
scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line,"_processed_filtered_type1_supercluster_signatures_aucell_scores.rds"))
thresholds <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line,"_processed_filtered_type1_supercluster_signatures_aucell_thresholds.rds"))



supercluster1_active_cells <- rownames(scores)[scores[,1] > scores[,2]]
supercluster2_active_cells <- rownames(scores)[scores[,1] < scores[,2]]

# supercluster1_active_cells <- rownames(scores)[scores[,1] > thresholds$threshold[1]]
# supercluster2_active_cells <- rownames(scores)[scores[,2] > thresholds$threshold[2]]

both_active <- intersect(supercluster1_active_cells,supercluster2_active_cells)

supercluster1_active_cells <- supercluster1_active_cells[!supercluster1_active_cells %in% both_active]
supercluster2_active_cells <- supercluster2_active_cells[!supercluster2_active_cells %in% both_active]

if(length(supercluster1_active_cells) > length(supercluster2_active_cells)){
  supercluster1_proliferation_scores <- data@meta.data %>% 
    filter(cell %in% sample(supercluster1_active_cells, length(supercluster2_active_cells))) %>% 
    pull(proliferation_index)
  
  supercluster2_proliferation_scores <- data@meta.data %>% 
    filter(cell %in% supercluster2_active_cells) %>% 
    pull(proliferation_index)
  
} else {
  supercluster1_proliferation_scores <- data@meta.data %>% 
    filter(cell %in% supercluster1_active_cells) %>% 
    pull(proliferation_index)
  
  supercluster2_proliferation_scores <- data@meta.data %>% 
    filter(cell %in% sample(supercluster2_active_cells, length(supercluster1_active_cells))) %>% 
    pull(proliferation_index)
}


wilcox_res <- wilcox.test(supercluster1_proliferation_scores,supercluster2_proliferation_scores)
plot_title <- paste0(curr_cell_line, " (pval: ", sprintf("%.5f",wilcox_res$p.value),")")
boxplot(supercluster1_proliferation_scores,supercluster2_proliferation_scores, main=plot_title,
        xlab='',ylab="Proliferation Index",
        names = c("Supercluster 1 Active",'Supercluster 2 Active'))








