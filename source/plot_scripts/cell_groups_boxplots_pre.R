library(ggpubr)
library(Seurat)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)


# curr_cell_line <- "A549"
# ident_to_use <- "Cluster"
# geneset_to_use <- "raj_watermelon_resistance_signature"


create_df_for_boxplots <- function(data, ident_to_use, geneset_to_use){
  
  Idents(data) <- ident_to_use
  
  num_groups <- nlevels(data)
  
  groups <- levels(data)
  
  cat(geneset_to_use, "\n")
  
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", geneset_to_use, "_aucell_scores.rds"))
  threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", geneset_to_use, "_aucell_thresholds.rds"))
  
  # Create list of lists containing active cells for each geneset
  all_active_cells <- list()
  
  curr_geneset <- colnames(scores)[1]
  # for(curr_geneset in colnames(scores)){
  
  cat(curr_geneset, "\n")
  active_cells <- rownames(scores)[scores[,curr_geneset] > threshold$threshold[threshold$gene_set == curr_geneset]]
  
  # all_active_cells <- append(all_active_cells, list(curr_active_cells))
  
  data <- AddMetaData(data, metadata = scores[,curr_geneset], col.name = "curr_geneset")
  
  
  df <- data@meta.data %>% 
    dplyr::select(ident_to_use, curr_geneset, dose, treatment_stage) %>% 
    group_by(eval(parse(text=ident_to_use))) %>% 
    mutate("group_median_score" = median(curr_geneset))
  # }
  
  return(df)    
}

###################################################################################

curr_cell_line <- "K562"

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

RACs <- list(c(4,9,12,13,14,16,18),c(4,5,11),c(5,8,12,13,17))
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


data$cell_group <- ifelse(data$treatment_stage == "pre", "pre", data$cell_group)

df <- create_df_for_boxplots(data, "cell_group", "raj_watermelon_resistance_signature")

plot_title <- paste0("Resistance Signature AUCell Scores in ", curr_cell_line)

df$cell_group <- ifelse(df$cell_group == 0, "Non-RAC", paste0("RAC Type ", df$cell_group))

df$cell_group <- ifelse(df$cell_group == "RAC Type pre", "pre", df$cell_group)

df$cell_group <- factor(df$cell_group, levels = c("RAC Type 1","RAC Type 2","Non-RAC","pre"))


p <- ggboxplot(df, x = "cell_group", y = "curr_geneset",fill = "cell_group")

my_comparisons <- list( c("Non-RAC", "RAC Type 2"), c("Non-RAC", "RAC Type 1"),c("RAC Type 1","RAC Type 2"),c("RAC Type 2","pre"),c("Non-RAC","pre"))


# png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/resistance_score_boxplots/", curr_cell_line,"_cell_group_boxplots.png"),
# width=800, height=600)

p + stat_compare_means(comparisons = my_comparisons,label = "p.format", method = "wilcox", label.x = 2.2, size=8)+
  ggtitle(plot_title)+
  xlab("")+
  ylab("Score")+
  scale_fill_manual(values=c("red","orangered", "lightblue","gray"),name = "Cell Groups")+
  theme(legend.position="right",
        title = element_text(size=20, face = "bold"),
        axis.text = element_text(size=20),
        legend.text = element_text(size=20),
        legend.title = element_text(size=24))

dev.off()
