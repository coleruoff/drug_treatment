setwd("/data/ruoffcj/projects/drug_treatment/")
library(Seurat)
library(ComplexHeatmap)
library(tidyverse)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

A549.data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/A549_processed_filtered.rds"))
K562.data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/K562_processed_filtered.rds"))
MCF7.data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/MCF7_processed_filtered.rds"))

#################################################################################
# Cluster based on mean expression of variable genes
#################################################################################

A549.data <- FindVariableFeatures(A549.data)
A549_variable_genes <- VariableFeatures(A549.data)

K562.data <- FindVariableFeatures(K562.data)
K562_variable_genes <- VariableFeatures(K562.data)

MCF7.data <- FindVariableFeatures(MCF7.data)
MCF7_variable_genes <- VariableFeatures(MCF7.data)

all_variable_genes <- unique(c(A549_variable_genes, K562_variable_genes, MCF7_variable_genes))

all_variable_genes <- all_variable_genes[all_variable_genes %in% rownames(A549.data)]
all_variable_genes <- all_variable_genes[all_variable_genes %in% rownames(K562.data)]
all_variable_genes <- all_variable_genes[all_variable_genes %in% rownames(MCF7.data)]

# length(all_variable_genes)

# Remove cell line specific markers
# cell_line_de <- readRDS("processed_data/cell_line_de.rds")
# 
# A549_cell_line_markers <- cell_line_de %>% 
#   filter(cluster == "A549" & p_val_adj < 0.05) %>% 
#   arrange(desc(avg_log2FC)) %>% 
#   pull(gene) %>% 
#   head(200)
# 
# K562_cell_line_markers <- cell_line_de %>% 
#   filter(cluster == "K562" & p_val_adj < 0.05) %>% 
#   arrange(desc(avg_log2FC)) %>% 
#   pull(gene) %>% 
#   head(200)
# 
# 
# MCF7_cell_line_markers <- cell_line_de %>% 
#   filter(cluster == "MCF7" & p_val_adj < 0.05) %>% 
#   arrange(desc(avg_log2FC)) %>% 
#   pull(gene) %>% 
#   head(200)
# 
# all_variable_genes <- all_variable_genes[!all_variable_genes %in% A549_cell_line_markers]
# 
# all_variable_genes <- all_variable_genes[!all_variable_genes %in% K562_cell_line_markers]
# 
# all_variable_genes <- all_variable_genes[!all_variable_genes %in% MCF7_cell_line_markers]

cell_lines <- c("A549", "K562", "MCF7")

heatmap <- matrix(NA, ncol=49, nrow=length(all_variable_genes))
colnames(heatmap) <- rep("", 49)

j <- 1
for(cell_line in cell_lines){
  
  cat(cell_line,"\n")
  data <- eval(parse(text=paste0(cell_line, ".data")))
  
  #total_mean <- rowMeans(data@assays$originalexp@data[all_variable_genes,])
  
  clusters <- 1:nlevels(data$Cluster)
  
  for(curr_cluster in clusters){
    cat(curr_cluster,"\n")
    
    cluster_mean <- rowMeans(data@assays$RNA$data[all_variable_genes, data$Cluster == curr_cluster])
    
    total_mean <- rowMeans(data@assays$RNA$data[all_variable_genes, data$Cluster != curr_cluster])
    
    #heatmap[,j] <- cluster_mean-total_mean
    heatmap[,j] <- total_mean-cluster_mean
    colnames(heatmap)[j] <- paste0(cell_line,"_",curr_cluster)
    
    j <- j + 1
  }
}

# active_clusters <- c("A549_4","A549_9","A549_13", "K562_5","K562_11","MCF7_5","MCF7_8","MCF7_13","MCF7_16","MCF7_17")
# 
# emergent_cluster <- c("A549_14","A549_15","A549_16", "A549_17","A549_18","A549_19","K562_9","MCF7_13","MCF7_15","MCF7_18")
# 
# clusters_of_interest <- emergent_cluster

cor_heatmap <- cor(heatmap, method = "spearman")

# row_idx <- which(rownames(cor_heatmap) %in% clusters_of_interest)
# fontsizes <- rep(10, nrow(cor_heatmap))
# fontsizes[row_idx] <- 18figure
# fontcolors <- rep('black', nrow(cor_heatmap))
# fontcolors[row_idx] <- 'red'
# fontfaces <- rep('plain',nrow(cor_heatmap))
# fontfaces[row_idx] <- 'bold'
# 
# rowAnno <- rowAnnotation(rows = anno_text(rownames(cor_heatmap), gp = gpar(fontsize = fontsizes, fontface = fontfaces, col = fontcolors)))

##################################################################################
# Global RAC supercluster plot
##################################################################################

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- paste0("A549_",RACs[["A549"]])
clusters_of_interest <- append(clusters_of_interest, paste0("K562_",RACs[["K562"]]))
clusters_of_interest <- append(clusters_of_interest, paste0("MCF7_",RACs[["MCF7"]]))

rac_ha <- rowAnnotation(RAC = c(ifelse(colnames(cor_heatmap) %in% clusters_of_interest,"RAC","Non-RAC")),
                        col = list(RAC = c("RAC" = "orange", "Non-RAC" = "lightblue")),show_annotation_name = F,
                        annotation_legend_param = list(title_gp=gpar(fontsize=30), grid_height=unit(2,"cm"),grid_width=unit(2,"cm"),
                                                       title="Cluster Type", labels_gp = gpar(fontsize = 16)))


ht <- Heatmap(cor_heatmap, name="Spearman\nCorrelation", cluster_rows = T, cluster_columns = T,
              right_annotation = rac_ha,column_title = "", column_title_side = "bottom",
              row_title_side = "left", row_title_gp = gpar(fontsize=20),
              column_names_rot = 45, 
              row_names_gp = gpar(fontsize=20),
              column_names_gp = gpar(fontsize=18),
              heatmap_legend_param = list(title_gp = gpar(fontsize = 30),legend_height = unit(5, "cm"), grid_width=unit(2,"cm"),
                                          labels_gp = gpar(fontsize = 16)))

png(paste0("/data/ruoffcj/projects/drug_treatment/final_figures/figure_3a.png"),
    width = 25,height=25, units = 'in',res = 300)

draw(ht, column_title="Correlations of Cluster Differential Mean Expression of Most Variable Genes", 
     column_title_gp = gpar(fontsize = 35, fontface = "bold"),  padding = unit(c(6, 10, 10, 2), "mm"),
     heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend=T)

dev.off()


