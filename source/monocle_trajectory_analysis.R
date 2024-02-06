library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(tidyverse)

cell_lines <- c("A549", "K562", "MCF7")
RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")

curr_cell_line <- "A549"

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/",curr_cell_line,"_processed_filtered.rds"))

scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
clusters_of_interest <- RACs[[curr_cell_line]]

#Add metadata for RAC and Cell Group
data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, "1", ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, paste0(data$Cluster, "_1"), ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), paste0(data$Cluster, "_2"), paste0(data$Cluster, "_0"))), col.name = "cell_cluster_group")



all_type1 <- c()
all_type2 <- c()

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

for(curr_cluster in RACs[[curr_cell_line]]){
  
  cat(curr_cluster, "\n")
  cluster_data <- data[,data$Cluster  == curr_cluster]
  
  
  cluster_data <- cluster_data[,1:1000]
  
  cds <- as.cell_data_set(cluster_data)
  
  fData(cds)$gene_short_name <- rownames(fData(cds))
  
  recreate.partition <- c(rep(1,length(cds@colData@rownames)))
  names(recreate.partition) <- cds@colData@rownames
  recreate.partition <- as.factor(recreate.partition)
  
  cds@clusters$UMAP$partitions <- recreate.partition
  
  list_cluster <- cluster_data$cell_group
  cds@clusters$UMAP$clusters <- list_cluster
  
  cds@int_colData@listData$reducedDims$UMAP <- cluster_data@reductions$UMAP@cell.embeddings
  
  cds <- learn_graph(cds,use_partition = F)
  
  cds <- order_cells(cds, reduction_method='UMAP', root_cells = colnames(cds)[clusters(cds) == 1])
  
  cds$monocle3_pseudotime <- range01(pseudotime(cds))
  data.pseudo <- as.data.frame(colData(cds))
  
  
  
  
  
  type1_times <- data.pseudo$monocle3_pseudotime[data.pseudo$cell_group==1]
  type2_times <- data.pseudo$monocle3_pseudotime[data.pseudo$cell_group==2]
  
  all_type1 <- append(all_type1,type1_times)
  all_type2 <- append(all_type2,type2_times)
  
}


boxplot(all_type1,all_type2)


min(all_type1)
max(all_type1)



cluster_data <- data[,data$Cluster == 4]

cds <- as.cell_data_set(cluster_data)

fData(cds)$gene_short_name <- rownames(fData(cds))

recreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds@clusters$UMAP$partitions <- recreate.partition

list_cluster <- cluster_data$cell_group
cds@clusters$UMAP$clusters <- list_cluster

cds@int_colData@listData$reducedDims$UMAP <- cluster_data@reductions$UMAP@cell.embeddings



cds <- learn_graph(cds,use_partition = F)


plot_cells(cds)


cds <- order_cells(cds, reduction_method='UMAP', root_cells = colnames(cds)[clusters(cds) == 1])


plot_cells(cds,
           color_cells_by = 'pseudotime')


cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))


ggplot(data.pseudo, aes(monocle3_pseudotime,resistance_stage))+
  geom_boxplot()




type1_times <- data.pseudo$monocle3_pseudotime[data.pseudo$cell_group==1]
type2_times <- data.pseudo$monocle3_pseudotime[data.pseudo$cell_group==2]


wilcox.test(type1_times,type2_times)


deg_traj <- graph_test(cds, neighbor_graph = 'principal_graph', cores=4)


traj_genes <- deg_traj %>% 
  filter(q_value < 0.05 & status == "OK") %>% 
  arrange(q_value) %>% 
  pull(gene_short_name)

traj_genes <- traj_genes[1:500]







