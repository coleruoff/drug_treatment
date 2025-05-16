args = commandArgs(trailingOnly=TRUE)
# dataDirectory <- paste0(args[1],"final_data/")
# plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("revision_scripts/drug_treatment_functions.R")

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_figures/"

################################################################################

cell_lines <- c("A549", "K562", "MCF7")

# all_data <- readRDS(paste0(dataDirectory, "processed_data/sciplex_data/all_cell_lines_data.rds"))
# 
# A549.data <- all_data[["A549"]]
# K562.data <- all_data[["K562"]]
# MCF7.data <- all_data[["MCF7"]]

all_racs <-  readRDS(paste0(dataDirectory, "processed_data/all_racs.rds"))

A549.data <- readRDS(paste0(dataDirectory,"processed_data/sciplex_data/A549_processed_filtered.rds"))
K562.data <- readRDS(paste0(dataDirectory,"processed_data/sciplex_data/K562_processed_filtered.rds"))
MCF7.data <- readRDS(paste0(dataDirectory,"processed_data/sciplex_data/MCF7_processed_filtered.rds"))

clusters_of_interest <- all_racs[["A549"]]
A549.data <- AddMetaData(A549.data, metadata = ifelse(A549.data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")

clusters_of_interest <- all_racs[["K562"]]
K562.data <- AddMetaData(K562.data, metadata = ifelse(K562.data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")

clusters_of_interest <- all_racs[["MCF7"]]
MCF7.data <- AddMetaData(MCF7.data, metadata = ifelse(MCF7.data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")


A549.data <- A549.data[,A549.data$rac == "rac" & A549.data$treatment_stage == "post"]

K562.data <- K562.data[,K562.data$rac == "rac" & K562.data$treatment_stage == "post"]

MCF7.data <- MCF7.data[,MCF7.data$rac == "rac" & MCF7.data$treatment_stage == "post"]

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

length(all_variable_genes)
################################################################################


num_cols <- nlevels(A549.data)+nlevels(K562.data)+nlevels(MCF7.data)

heatmap <- matrix(NA, ncol=num_cols, nrow=length(all_variable_genes))
colnames(heatmap) <- rep("", num_cols)

j <- 1
for(cell_line in cell_lines){
  
  cat(cell_line,"\n")
  data <- eval(parse(text=paste0(cell_line, ".data")))
  
  clusters <- unique(data$Cluster)
  
  for(curr_cluster in clusters){
    cat(curr_cluster,"\n")
    
    
    cluster_mean <- rowMeans(data[["RNA"]]$data[all_variable_genes, data$Cluster == curr_cluster])
    
    total_mean <- rowMeans(data[["RNA"]]$data[all_variable_genes,])
    
    heatmap[,j] <- total_mean-cluster_mean
    colnames(heatmap)[j] <- paste0(cell_line,"_",curr_cluster)
    
    j <- j + 1
  }
}

cor_heatmap <- cor(heatmap, method = "spearman")

sil_data <- fviz_nbclust(cor_heatmap, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")


jpeg(paste0(plotDirectory,"figure_S2a.jpg"), width=100, height = 80, units = "mm", res = 1000)
print(sil_data)
dev.off()

num_superclusters <- which(sil_data$data$y == max(sil_data$data$y))

# Heatmap(cor_heatmap, name="Spearman\nCorrelation", cluster_rows = T, cluster_columns = T,
#         column_title_side = "top",
#         row_title_side = "left", row_title_gp = gpar(fontsize=20),
#         column_names_rot = 90, row_title = c(),
#         row_names_gp = gpar(fontsize=8),
#         column_names_gp = gpar(fontsize=8),
#         heatmap_legend_param = list(title_gp = gpar(fontsize = 8),legend_height = unit(4, "mm"), grid_width=unit(4,"mm"),
#                                     labels_gp = gpar(fontsize = 6)))


ht <- Heatmap(cor_heatmap, name="Spearman\nCorrelation", cluster_rows = T, cluster_columns = T,
        column_title_side = "top",row_split = num_superclusters,column_split = num_superclusters,
        row_title_side = "left", row_title_gp = gpar(fontsize=6),
        column_names_rot = 90, row_title = c(),
        # column_title = c("Supercluster 1","Supercluster 2","Supercluster 3"),
        column_title_gp = gpar(fontsize=6),
        row_names_gp = gpar(fontsize=6),
        column_names_gp = gpar(fontsize=6),
        heatmap_legend_param = list(title_gp = gpar(fontsize = 6),legend_height = unit(2, "mm"), grid_width=unit(2,"mm"),
                                    labels_gp = gpar(fontsize = 4),legend_gp = gpar(lwd = .5)))


ht <- draw(ht,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend=T)






jpeg(paste0(plotDirectory,"figure_2a.jpg"), width=100, height = 80, units = "mm", res = 1000)
print(ht)
dev.off()



# colnames(cor_heatmap)


cluster_components <- row_order(ht)


# all(cor_heatmap[cluster_components[[1]],cluster_components[[1]]]>0)
# all(cor_heatmap[cluster_components[[2]],cluster_components[[2]]]>0)
# all(cor_heatmap[cluster_components[[3]],cluster_components[[3]]]>0)
# all(cor_heatmap[cluster_components[[4]],cluster_components[[4]]]>0)
# all(cor_heatmap[cluster_components[[5]],cluster_components[[5]]]>0)


supercluster_components <- list()
for(curr_cluster in cluster_components){
  
  
  supercluster_components <- append(supercluster_components, list(colnames(cor_heatmap)[curr_cluster]))
}


# supercluster_components <- supercluster_components[c(2)]

names(supercluster_components) <- paste0("supercluster", 1:length(supercluster_components))


new_components <- list()
for(i in 1:length(supercluster_components)){
  
  curr_supercluster <- list()
  
  for(curr_cell_line in cell_lines){
    
    curr_clusters <- supercluster_components[[i]][grepl(curr_cell_line,supercluster_components[[i]])]
    
    curr_supercluster[[curr_cell_line]] <- as.numeric(gsub(paste0(curr_cell_line,"_"),"",curr_clusters))
  }
  
  new_components[[paste0("supercluster",i)]] <- curr_supercluster
}

supercluster_components <- new_components


saveRDS(supercluster_components, paste0(dataDirectory, "processed_data/supercluster_components.rds"))

