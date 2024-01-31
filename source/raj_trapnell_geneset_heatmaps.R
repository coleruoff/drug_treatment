library(ComplexHeatmap)
library(Seurat)
library(circlize)
library(tidyverse)

geneset_group_matrix <- function(data, data_name, ident_to_use, genesets_to_use){
  
  Idents(data) <- ident_to_use
  
  num_groups <- nlevels(data)
  
  groups <- levels(data)
  
  cat(genesets_to_use, "\n")
  
  scores <- readRDS(paste0(dataDirectory, "data/aucell_score_objects/", data_name, "_", genesets_to_use, "_aucell_scores.rds"))
  # threshold <- readRDS(paste0(dataDirectory, "data/aucell_score_objects/",curr_cell_line, "_processed_filtered_", genesets_to_use, "_aucell_thresholds.rds"))
  
  # scores <- scale(scores)
  
  geneset_names <- colnames(scores)
  
  score_matrix <- matrix(NA, ncol=length(groups), nrow=length(geneset_names))
  pvalue_matrix <- matrix(NA, ncol=length(groups), nrow=length(geneset_names))
  
  for(curr_group in groups){
    
    cat(curr_group, "\n")
    
    curr_cell_names <- colnames(data)[data[[ident_to_use]] == curr_group]
    rest_names <- colnames(data)[data[[ident_to_use]] != curr_group]
    
    for(curr_geneset in geneset_names){
      
      curr_scores <- scores[curr_cell_names, curr_geneset]
      
      rest_scores <- scores[rest_names, curr_geneset]
      
      wilcox_res <- wilcox.test(curr_scores,rest_scores)
      
      i <- which(curr_geneset == geneset_names)
      j <- which(curr_group == groups)
      
      score_matrix[i,j] <- mean(curr_scores)
      pvalue_matrix[i,j] <- wilcox_res$p.value
      
    }
  }
  
  scaled_matrix <- t(scale(t(score_matrix)))
  
  colnames(scaled_matrix) <- groups
  rownames(scaled_matrix) <- geneset_names
  
  colnames(score_matrix) <- groups
  rownames(score_matrix) <- geneset_names
  
  return(list(scaled_matrix, pvalue_matrix,score_matrix))
}

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/"
################################################################################
# Set genesets to use as expression space for trapnell data
genesets_name <- "MPs_without_resistance"

curr_cell_line <- "A549"
data <- readRDS(paste0(dataDirectory, "data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")

################################################################################
# Create trapnell geneset heatmap/expression space values

trapnell_mat <- geneset_group_matrix(data, paste0(curr_cell_line,"_processed_filtered"), "Cluster", genesets_name)

Heatmap(trapnell_mat[[1]])

################################################################################
# Read in raj data
data2 <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/diverse_clonal_fates_data/raj_resistant_breast_processed.rds"))

# Set genesets to use as expression space for raj data
genesets_name <- "ITH_meta_programs"

# Create trapnell geneset heatmap/expression space values
raj_mat <- geneset_group_matrix(data2, "raj_resistant_breast_processed", "seurat_clusters", genesets_name)

################################################################################
# Z-score heatmaps within each cluster

trapnell_matrix <- scale(trapnell_mat[[3]])
colnames(trapnell_matrix) <- paste0("trapnell_",colnames(trapnell_matrix))
# Heatmap(trapnell_matrix)

raj_matrix <- scale(raj_mat[[3]])
colnames(raj_matrix) <- paste0("raj_",1:ncol(raj_matrix))
# Heatmap(raj_matrix)

################################################################################
# Calculate distances between all raj and trapnell clusters

dist_mat <- as.matrix(dist(rbind(t(raj_matrix),t(trapnell_matrix)), method="euclidean"))
subset_dist_mat <- dist_mat[rownames(dist_mat) %in% colnames(raj_matrix),colnames(dist_mat) %in% colnames(trapnell_matrix)]

# Select current RACs 
# clusters_of_interest <- paste0("trapnell_",RACs[[curr_cell_line]])

# Select RAC resistant subpops
clusters_of_interest <- paste0("trapnell_",RACs[[curr_cell_line]])

# Select distances for RACs and non-RACs
distances_to_RACs <- as.vector(subset_dist_mat[,clusters_of_interest])
rest_distances <- as.vector(subset_dist_mat[,!colnames(subset_dist_mat) %in% clusters_of_interest])

################################################################################
# Create dataframe for box plot
df <- data.frame(cbind(c(distances_to_RACs,rest_distances),c(rep("RACs", length(distances_to_RACs)),rep("non-RACs",length(rest_distances)))))
colnames(df) <- c("value","group")
df$group <- factor(df$group, levels = c("RACs","non-RACs"))

#Convert distances to similarity score
df$value <- as.numeric(df$value)
df$value <- 1/(exp(df$value))

################################################################################
# wilcoxon test for pvalue
RAC_similarity <- df %>% 
  filter(group=="RACs") %>% 
  pull(value)

rest_similarity <- df %>% 
  filter(group=="non-RACs") %>% 
  pull(value)

wilcox_res <- wilcox.test(RAC_similarity,rest_similarity)

################################################################################
# Plot and save figure as file

clean_geneset_name <- "ITH Meta-programs"

plot_title <- paste0(curr_cell_line, " Clusters Distances to Raj Resistant Clusters in ", clean_geneset_name, " Expression Space (p value: ", sprintf("%.6f", wilcox_res$p.value),")")

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/distance_to_raj_boxplots/",curr_cell_line,"_",genesets_name,".png"), width = 1500,height = 1000)
plot(ggplot(df)+
  geom_boxplot(aes(x=group,y=value, fill=group))+
  xlab("Cluster Group")+
  ylab("Distances")+
  ggtitle(plot_title)+
  NoLegend()+
  theme(axis.text.x = element_text(size=20),
        axis.title.y  = element_text(size=20),
        title = element_text(size=20),
        axis.text.y= element_text(size=20)))

dev.off()


################################################################################
# Read in resistance scores and set up meta data for RACs and resistant subpops
scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))

active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- RACs[[curr_cell_line]]

data <- AddMetaData(data, metadata = ifelse(colnames(data) %in% active_cell_names, "active","inactive"), col.name = "resistant_active")
data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")

data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & data$resistant_active == "active", "resistant","nonresistant"), col.name = "resistant")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & data$resistant_active == "active", paste0(data$Cluster, "_resistant"),data$Cluster), col.name = "resistant_cluster")

################################################################################
# Create trapnell geneset heatmap/expression space values

trapnell_mat <- geneset_group_matrix(data, paste0(curr_cell_line,"_processed_filtered"), "resistant_cluster", genesets_name)

Heatmap(trapnell_mat[[1]])

################################################################################
# Z-score heatmaps within each cluster

trapnell_matrix <- scale(trapnell_mat[[3]])
colnames(trapnell_matrix) <- paste0("trapnell_",colnames(trapnell_matrix))
# Heatmap(trapnell_matrix)

################################################################################
# Calculate distances between all raj and trapnell clusters

dist_mat <- as.matrix(dist(rbind(t(raj_matrix),t(trapnell_matrix)), method="euclidean"))
subset_dist_mat <- dist_mat[rownames(dist_mat) %in% colnames(raj_matrix),colnames(dist_mat) %in% colnames(trapnell_matrix)]

# Select RAC resistant subpops
clusters_of_interest <- paste0("trapnell_",RACs[[curr_cell_line]],"_resistant")

# Select distances for RACs and non-RACs
distances_to_RACs <- as.vector(subset_dist_mat[,clusters_of_interest])
rest_distances <- as.vector(subset_dist_mat[,!colnames(subset_dist_mat) %in% clusters_of_interest])

################################################################################
# Create dataframe for box plot
df <- data.frame(cbind(c(distances_to_RACs,rest_distances),c(rep("Resistant Cells", length(distances_to_RACs)),rep("Rest",length(rest_distances)))))
colnames(df) <- c("value","group")
df$group <- factor(df$group, levels = c("Resistant Cells","Rest"))

#Convert distances to similarity score
df$value <- as.numeric(df$value)
df$value <- 1/(exp(df$value))

################################################################################
# wilcoxon test for pvalue
resistant_similarity <- df %>% 
  filter(group=="Resistant Cells") %>% 
  pull(value)

rest_similarity <- df %>% 
  filter(group=="Rest") %>% 
  pull(value)

wilcox_res <- wilcox.test(resistant_similarity,rest_similarity)

################################################################################
# Plot and save figure as file

clean_geneset_name <- "ITH Meta-programs"

plot_title <- paste0(curr_cell_line, " Resistant Subpopulations Similarity to Raj Resistant Clusters in ", clean_geneset_name, " Expression Space (p value: ", sprintf("%.6f", wilcox_res$p.value),")")

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/figures/distance_to_raj_boxplots/",curr_cell_line,"_",genesets_name,"_resistant_cells.png"), width = 1500,height = 1000)
plot(ggplot(df)+
       geom_boxplot(aes(x=group,y=value, fill=group))+
       xlab("")+
       ylab("Similarity Score")+
       ggtitle(plot_title)+
       NoLegend()+
       theme(axis.text.x = element_text(size=20),
             axis.title.y  = element_text(size=20),
             title = element_text(size=20),
             axis.text.y= element_text(size=20)))

dev.off()

################################################################################
################################################################################
# Gaussian Mixture Model

library(mclust)


i <- 1
for(i in 1:19){
  gmm <- Mclust(t(subset_dist_mat),G=i)
  
  cat(i, ": ", max(gmm$BIC[!is.na(gmm$BIC)]),"\n")
  
  
}



gmm <- Mclust(t(subset_dist_mat),G=2)

summary(gmm)

group1 <- names(gmm$classification)[gmm$classification == 1]
group2 <- names(gmm$classification)[gmm$classification == 2]
group3 <- names(gmm$classification)[gmm$classification == 3]
group4 <- names(gmm$classification)[gmm$classification == 4]

group1_distances <- as.vector(subset_dist_mat[,group1])
group2_distances <- as.vector(subset_dist_mat[,group2])
group3_distances <- as.vector(subset_dist_mat[,group3])
group4_distances <- as.vector(subset_dist_mat[,group4])

boxplot(group1_distances,group2_distances,group3_distances,group4_distances)

gmm$classification


################################################################################











dim(raj_matrix[1:2,1])

dev.off()
# combined_mat <- cbind(trapnell_matrix,raj_matrix)



raj_scores <- readRDS(paste0(dataDirectory, "data/aucell_score_objects/raj_resistant_breast_processed_", genesets_name, "_aucell_scores.rds"))

combined_mat <- cbind(trapnell_matrix,scale(colMeans(raj_scores)))



RACs <- list(c(4,9,12,13,14,16,18),c(4,5,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- paste0("trapnell_",RACs[[curr_cell_line]])

rac_ha <- HeatmapAnnotation(RAC = c(ifelse(colnames(combined_mat) %in% clusters_of_interest,"RAC",ifelse(colnames(combined_mat) %in% colnames(trapnell_matrix),"Non-RAC","Raj"))),
                            col = list(RAC = c("RAC" = "darkgreen", "Non-RAC" = "lightblue", "Raj"="red")))

Heatmap(combined_mat, bottom_annotation = rac_ha)

col_fun = colorRamp2(c(1, .9), c("red", "white"))
Heatmap(1-(as.matrix(dist(t(combined_mat))))/100, bottom_annotation = rac_ha,col=col_fun)
