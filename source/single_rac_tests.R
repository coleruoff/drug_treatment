library(tidyverse)
library(Seurat)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)

curr_cell_line <- "MCF7"

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
# RACs <- list(c(4,9,12,13,14,16,18),c(4,5,11),c(5,8,12,13,17))
RACs <- list(c(12),c(11),c(8))
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


DimPlot(data, group.by='Cluster',split.by = "rac")


Idents(data) <- data$rac
de_res <- FindAllMarkers(data)



rac_genes_up <- de_res %>% 
  filter(cluster == "rac" & p_val_adj < 0.05 & avg_log2FC > 2) %>% 
  arrange(desc(avg_log2FC)) %>% 
  pull(gene)

rac_genes_up <- de_res %>% 
  filter(cluster == "rac" & p_val_adj < 0.05 & avg_log2FC < -2) %>% 
  arrange(desc(avg_log2FC)) %>% 
  pull(gene)



m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, human_gene_symbol)
# hallmark_uni <- unique(c(unlist(all_signatures), m_t2g$human_gene_symbol))

mp_t2g <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ith_meta_programs_t2g.rds")


curr_geneset <- rac_genes

cell_line_universes <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cell_line_universes.rds")
universe_to_use <- cell_line_universes[[curr_cell_line]]

# Hallmarks Enrichment
hallmark_enrichment_results <- enricher(curr_geneset, TERM2GENE=m_t2g,
                                        universe = universe_to_use)

dotplot(hallmark_enrichment_results)

hallmark_enrichment_results@result %>% 
  filter(qvalue < 0.05) %>% 
  pull(Description)

# MPs Enrichment
mp_enrichment_results <- enricher(curr_geneset, TERM2GENE=mp_t2g,
                                  universe = universe_to_use)

dotplot(mp_enrichment_results)

mp_enrichment_results@result %>% 
  filter(qvalue < 0.05) %>% 
  pull(Description)

# GO Functional Enrichment
ego <- enrichGO(gene          = curr_geneset,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = .01,
                qvalueCutoff  = .05,
                readable      = TRUE,
                keyType = "SYMBOL",
                universe = universe_to_use)


dotplot(ego)






