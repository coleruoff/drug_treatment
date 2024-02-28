library(Seurat)
library(presto)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(GSVA)
library(tidyverse)


cell_lines <- c("A549","K562","MCF7")

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- cell_lines

emergent <- list(c(14:19),c(9,11),c(13,15,18))
names(emergent) <- cell_lines

curr_cell_line <- "MCF7"

all_data <- list()
for(curr_cell_line in cell_lines){
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/",curr_cell_line, "_processed_filtered.rds"))
  
  #read in DR signature scores and set active cells
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
  threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
  active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
  clusters_of_interest <- RACs[[curr_cell_line]]
  
  #Add metadata for RAC and Cell Group
  data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, "1", ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, paste0(data$Cluster, "_1"), ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), paste0(data$Cluster, "_2"), paste0(data$Cluster, "_0"))), col.name = "cell_cluster_group")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & data$Cluster %in% emergent[[curr_cell_line]], "emergent_rac", ifelse(data$rac == "rac" & (!data$Cluster %in% emergent[[curr_cell_line]]), "non_emergent_rac", "non_rac")), col.name = "emergent_rac")
  
  
  all_data[curr_cell_line] <- data
}


# Idents(data) <- data$emergent_rac

####################


supercluster1_proliferation_scores <- c(all_data[["A549"]]$proliferation_index[all_data[["A549"]]$Cluster == 9],
all_data[["K562"]]$proliferation_index[all_data[["K562"]]$Cluster == 5],
all_data[["MCF7"]]$proliferation_index[all_data[["MCF7"]]$Cluster == 8])


supercluster2_proliferation_scores <- c(all_data[["A549"]]$proliferation_index[all_data[["A549"]]$Cluster == 14],
  all_data[["K562"]]$proliferation_index[all_data[["K562"]]$Cluster == 9],
  all_data[["MCF7"]]$proliferation_index[all_data[["MCF7"]]$Cluster == 13])


boxplot(supercluster1_proliferation_scores,supercluster2_proliferation_scores,names=c("Supercluster 1","Supercluster 2"))


######################


ggboxplot(data@meta.data, x="Cluster",y="proliferation_index")


data <- AddMetaData(data, scores, col.name ="resistance_score")


ggboxplot(data@meta.data, x="emergent_rac",y="resistance_score",fill="emergent_rac")


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

ggboxplot(data@meta.data, x="emergent_rac",y="percent.mt",fill="emergent_rac")


temp <- data@meta.data %>% 
  dplyr::count(emergent_rac,dose)


temp$dose <- factor(temp$dose)

ggplot(temp)+
  geom_col(aes(x=emergent_rac,y=n,fill=dose),position = "dodge")


scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line,"_processed_filtered_kegg_endocytosis_aucell_scores.rds"))
data <- AddMetaData(data, scores,col.name = "endocytosis_score")

ggboxplot(data@meta.data, x="emergent_rac",y="endocytosis_score",fill="emergent_rac")













