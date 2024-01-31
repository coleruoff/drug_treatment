library(Seurat)
library(tidyverse)
library(fgsea)
library(clusterProfiler)
source("source/cole_functions.R")

################################################################################

curr_cell_line <- "A549"

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/",curr_cell_line, "_processed_filtered.rds"))
Idents(data) <- data$treatment_stage

non_emergent <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/trapnell_cluster_groups/nonemergent_clusters.rds")
non_emergent <- non_emergent[[curr_cell_line]]

################################################################################

cluster_treatment_de_genes <- list()

for(i in non_emergent){
  cat("Cluster ", i, ":\n")
  
  curr_de_res <- FindMarkers(data[,data$Cluster == i], ident.1="post")
  
  curr_cluster_de_genes <- curr_de_res %>% 
    rownames_to_column() %>%
    arrange(desc(avg_log2FC)) %>% 
    pull(avg_log2FC)
  
  names(curr_cluster_de_genes) <- curr_de_res %>% 
    rownames_to_column() %>%
    arrange(desc(avg_log2FC)) %>% 
    pull(rowname)
  
  
  cluster_treatment_de_genes <- append(cluster_treatment_de_genes, list(curr_cluster_de_genes))
}

names(cluster_treatment_de_genes) <- paste0(curr_cell_line,"_cluster_", non_emergent)

genesets_with_ranks <- cluster_treatment_de_genes

emergent_vs_pre_de <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_emergent_vs_pre_nonemergent_de.rds"))

get_ranks <- function(x){
  
  curr_cluster_de_genes <- x %>% 
    rownames_to_column() %>%
    arrange(desc(avg_log2FC)) %>% 
    pull(avg_log2FC)
  
  names(curr_cluster_de_genes) <- x %>% 
    rownames_to_column() %>%
    arrange(desc(avg_log2FC)) %>% 
    pull(rowname)
  
  return(curr_cluster_de_genes)
  
}

emergent_vs_pre_de_ranks <- sapply(emergent_vs_pre_de, FUN = function(x) get_ranks(x))


genesets_with_ranks <- append(genesets_with_ranks,emergent_vs_pre_de_ranks)
  
################################################################################

raj_de_genesets <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_breast_clusters_vs_control_de_genesets_500.rds")

raj_cluster_markers <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_breast_cluster_de_markers.rds")

genesets2 <- raj_de_genesets

################################################################################


genesets1 <- sapply(genesets_with_ranks, FUN = function(x) names(x[x > 0]))


result <- compare_genesets_fisher(genesets1, genesets2, background_genes = rownames(data))


plot_pretty_heatmap(log(result), plot_title,legend_title="log(OR)",col_fun=col_fun)






