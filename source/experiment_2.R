library(tidyverse)
library(ggpubr)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
cell_lines <- c("A549","K562","MCF7")

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")

curr_cell_line <- cell_lines[1]

cat(curr_cell_line,"\n")

#Read in cell line data
data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

#read in DR signature scores and set active cells
scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
clusters_of_interest <- RACs[[curr_cell_line]]

#Add metadata for RAC and Cell Group
data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, "1", ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, paste0(data$Cluster, "_1"), ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), paste0(data$Cluster, "_2"), paste0(data$Cluster, "_0"))), col.name = "cell_cluster_group")

Idents(data) <- data$rac

df <- data@meta.data

p <- ggboxplot(df, x = "rac", y = "g1s_score",fill = "rac")

my_comparisons <- list( c("rac", "nonrac"))

plot_title <- "G1S Score Distributions"

p <- p + stat_compare_means(comparisons = my_comparisons,label = "p.format", method = "wilcox", label.x = 2.2, size=8)+
  ggtitle(plot_title)+
  xlab("")+
  ylab("G1S Score")+
  scale_fill_manual(values=c("lightblue", "pink"),name = "Cell Groups")+
  theme(legend.position="right",
        title = element_text(size=20, face = "bold"),
        axis.text = element_text(size=20),
        legend.text = element_text(size=24),
        legend.title = element_text(size=26),
        legend.key.height = unit(1.5,"cm"),
        legend.key.width = unit(1.5,"cm"))

p
