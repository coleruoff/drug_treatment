library(tidyverse)
library(ComplexHeatmap)
library(AUCell)
source("source/cole_functions.R")

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")


rac_type1_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_type1_signatures.rds")
rac_type2_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_type2_signatures.rds")


Heatmap(calc_jaccard_matrix(rac_type1_signatures,rac_type2_signatures), cluster_rows = F,cluster_columns = F)

curr_cell_line <- "A549"

cat(curr_cell_line,"\n")

#Read in cell line data
data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
clusters_of_interest <- RACs[[curr_cell_line]]

#Add metadata for RAC and Cell Group
data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, "1", ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, paste0(data$Cluster, "_1"), ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), paste0(data$Cluster, "_2"), paste0(data$Cluster, "_0"))), col.name = "cell_cluster_group")



i <- 9
temp <- data[,data$Cluster == i]

signature1_name <- paste0(curr_cell_line,"_rac_subcluster",i,"_1_signature")
signature2_name <- paste0(curr_cell_line,"_rac_subcluster",i,"_2_signature")

signatures_to_use <- c(rac_type1_signatures[signature1_name],rac_type2_signatures[signature2_name])

cells_AUC <- AUCell_run(temp@assays$RNA@data, signatures_to_use)

temp <- AddMetaData(temp, t(cells_AUC@assays@data$AUC), col.name = colnames((t(cells_AUC@assays@data$AUC))))

df <- temp@meta.data %>% 
  select(Cluster,cell_group,signature1_name,signature2_name)

df <- df %>% 
  pivot_longer(!c(Cluster,cell_group),names_to = "type_geneset", values_to = "score")

ggplot(df)+
  geom_boxplot(aes(x=cell_group, y=score,fill=type_geneset))



###################################################

Idents(temp) <- temp$cell_group
de_res <- FindAllMarkers(temp)

type1_markers <- de_res %>% 
  filter(p_val_adj <0.05,cluster==1) %>% 
  arrange(avg_log2FC) %>% 
  pull(gene)


type2_markers <- de_res %>% 
  filter(p_val_adj <0.05,cluster==2) %>% 
  arrange(avg_log2FC) %>% 
  pull(gene)


shared_markers <- intersect(type1_markers,type2_markers)

type1_markers <- type1_markers[!type1_markers %in% shared_markers]
type2_markers <- type2_markers[!type2_markers %in% shared_markers]



shared_genes <- intersect(signatures_to_use[[1]],signatures_to_use[[2]])

signatures_to_use[[1]] <- signatures_to_use[[1]][!signatures_to_use[[1]] %in% shared_genes]
signatures_to_use[[2]] <- signatures_to_use[[2]][!signatures_to_use[[2]] %in% shared_genes]


signatures_to_use[[1]] <- unique(c(signatures_to_use[[1]],type1_markers))
signatures_to_use[[2]] <- unique(c(signatures_to_use[[2]],type2_markers))


cells_AUC <- AUCell_run(temp@assays$RNA@data, signatures_to_use)

temp <- AddMetaData(temp, t(cells_AUC@assays@data$AUC), col.name = colnames((t(cells_AUC@assays@data$AUC))))

df <- temp@meta.data %>% 
  select(Cluster,cell_group,signature1_name,signature2_name)

df <- df %>% 
  pivot_longer(!c(Cluster,cell_group),names_to = "type_geneset", values_to = "score")

ggplot(df)+
  geom_boxplot(aes(x=cell_group, y=score,fill=type_geneset))

