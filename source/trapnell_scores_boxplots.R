library(Seurat)
library(tidyverse)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")

#Set cell line
curr_cell_line <- "MCF7"
#Read in cell line data
data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
clusters_of_interest <- RACs[[curr_cell_line]]

data <- AddMetaData(data, metadata = ifelse(colnames(data) %in% active_cell_names, "active","inactive"), col.name = "resistant_active")
data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & data$resistant_active == "active", "resistant","nonresistant"), col.name = "resistant")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & data$resistant_active == "active", paste0(data$Cluster, "_resistant"),data$Cluster), col.name = "resistant_cluster")

df <- data@meta.data %>% 
  select(proliferation_index, Cluster,resistant,rac,resistant_cluster,dose,treatment_stage) %>% 
  filter(treatment_stage == "post")


df$dose <- as.character(df$dose)

p <- ggboxplot(df, x = "Cluster", y = "proliferation_index",fill="resistant",palette = "jco", short.panel.labs = FALSE)

# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format", method = "wilcox")+ggtitle(paste0(curr_cell_line, " Proliferation Scores"))

