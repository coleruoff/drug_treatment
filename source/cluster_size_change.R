library(Seurat)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(fgsea)
library(monocle3)

#Set cell line
curr_cell_line <- "A549"

#Read in cell line data
all_data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

clusters <- 1:nlevels(all_data)

#################################################################################
#PRE TO POST

classes <- as.character(unique(all_data$pathway_level_1))
classes <- classes[-which(classes=="Vehicle")]
# classes <- classes[-which(classes=="Others")]

heatmap <- matrix(1, nrow=length(classes), ncol = length(clusters))
colnames(heatmap) <- paste0("Cluster ", clusters)
rownames(heatmap) <- classes
pvals <- matrix(1, nrow=length(classes), ncol = length(clusters))

for(curr_class in classes){
  
  
  cat(curr_class, "\n")
  
  i <- which(curr_class == classes)
  
  total_pre <- sum(all_data$treatment_stage == 'pre')
  total_post <- sum(all_data$treatment_stage == 'post' & all_data$pathway_level_1 == curr_class)
  
  for(curr_cluster in clusters){
    
    cluster_pre <- sum(all_data$Cluster==curr_cluster & all_data$treatment_stage == 'pre')
    cluster_post <- sum(all_data$Cluster==curr_cluster & all_data$treatment_stage == 'post' & all_data$pathway_level_1 == curr_class)
    
    
    
    fisher_table <- matrix(c(cluster_post,total_post,cluster_pre,total_pre), nrow = 2,
                           dimnames =
                             list(c("curr cluster", "total"),
                                  c("post", "pre")))
    
    
    fisher_table <- fisher_table+1
    fisher_results <- fisher.test(fisher_table, alternative = "two.sided")
    
    
    
    
    j <- which(curr_cluster == clusters)
    
    heatmap[i,j] <- fisher_results$estimate
    pvals[i,j] <- fisher_results$p.value
    
    # if(cluster_post < .5*sum(all_data$pathway_level_1==curr_class)){
    #   pvals[i,j] <- fisher_results$p.value
    # } else{
    #   pvals[i,j] <- 1
    # }
    
  }
}

################################################################################
#PLOT HEATMAP
################################################################################


sig_heatmap <- ifelse(pvals>.05, 1, heatmap)
colnames(sig_heatmap) <- colnames(heatmap)
rownames(sig_heatmap) <- rownames(heatmap)

col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

ht <- Heatmap(log(sig_heatmap), name="log(Odds Ratio)\n",cluster_rows = F, cluster_columns =F, col = col_fun,
              column_names_rot = 45, row_names_max_width = unit(12, "cm"), column_names_max_height = unit(12, "cm"),
              column_names_gp = gpar(fontsize = 18),row_names_gp = gpar(fontsize = 16),
              
              heatmap_legend_param = list(legend_height=unit(4,"cm"), grid_width=unit(1,"cm"), title_gp=gpar(fontsize=20), labels_gp=gpar(fontsize=15)))



# file_name <- paste0("documents/figures/cluster_size_fc/", curr_cell_line, "/",curr_cell_line, "_cluster_size_fc.png")
# png(file_name, width= 1200, height= 1300, units="px")

draw(ht, padding = unit(c(2, 20, 2, 2), "mm"),column_title = paste0("Relative Cluster Size Change Between Treatment Stages Across Classes (", curr_cell_line, ")"), column_title_gp = gpar(fontsize = 26))

# dev.off()



################################################################################

pan_drug_heatmap <- matrix(NA, ncol=length(clusters),nrow=1)
pan_drug_pvals <- matrix(NA, ncol=length(clusters),nrow=1)

total_pre <- sum(all_data$treatment_stage == 'pre')
total_post <- sum(all_data$treatment_stage == 'post')

for(curr_cluster in clusters){
  
  cluster_pre <- sum(all_data$Cluster==curr_cluster & all_data$treatment_stage == 'pre')
  cluster_post <- sum(all_data$Cluster==curr_cluster & all_data$treatment_stage == 'post')
  
  
  
  fisher_table <- matrix(c(cluster_post,total_post,cluster_pre,total_pre), nrow = 2,
                         dimnames =
                           list(c("curr cluster", "total"),
                                c("post", "pre")))
  
  
  fisher_table <- fisher_table+1
  fisher_results <- fisher.test(fisher_table, alternative = "two.sided")
  
  
  
  
  j <- which(curr_cluster == clusters)
  
  pan_drug_heatmap[,j] <- fisher_results$estimate
  pan_drug_pvals[,j] <- fisher_results$p.value
}


colnames(pan_drug_heatmap) <- colnames(heatmap)
rownames(pan_drug_heatmap) <- "All Drugs"


pan_drug_ht <- Heatmap(log(pan_drug_heatmap), name="log(Odds Ratio)\n",cluster_rows = F, cluster_columns =F, col = col_fun,
                       column_names_rot = 45, row_names_max_width = unit(12, "cm"), column_names_max_height = unit(12, "cm"),
                       column_names_gp = gpar(fontsize = 16),row_names_gp = gpar(fontsize = 16),
                       
                       heatmap_legend_param = list(legend_height=unit(4,"cm"), grid_width=unit(1,"cm"), title_gp=gpar(fontsize=20), labels_gp=gpar(fontsize=15)))



ht_list <- ht %v% pan_drug_ht



# file_name <- paste0("results/cluster_size_change_heatmaps/by_class/",curr_cell_line, "_cluster_size_fc_class.png")
# png(file_name, width= 1200, height= 1300, units="px")

draw(ht_list, padding = unit(c(2, 20, 2, 2), "mm"),column_title = paste0("Relative Cluster Size Change Between Treatment Stages Across Drug Classes (", curr_cell_line, ")"), column_title_gp = gpar(fontsize = 26))


dev.off()



################################################################################


if(curr_cell_line == "A549"){
  
  A549_emergent_clusters <- list()
  for(i in 1:nrow(sig_heatmap)){
    A549_emergent_clusters <- append(A549_emergent_clusters,list((1:ncol(sig_heatmap))[sig_heatmap[i,] > 1.5]))
    
  }
  names(A549_emergent_clusters) <- rownames(sig_heatmap)
  
} else if(curr_cell_line == "K562"){
  
  K562_emergent_clusters <- list()
  for(i in 1:nrow(sig_heatmap)){
    K562_emergent_clusters <- append(K562_emergent_clusters,list((1:ncol(sig_heatmap))[sig_heatmap[i,] > 1.5]))
    
  }
  names(K562_emergent_clusters) <- rownames(sig_heatmap)
  
} else if (curr_cell_line == "MCF7"){
  MCF7_emergent_clusters <- list()
  for(i in 1:nrow(sig_heatmap)){
    MCF7_emergent_clusters <- append(MCF7_emergent_clusters,list((1:ncol(sig_heatmap))[sig_heatmap[i,] > 1.5]))
    
  }
  names(MCF7_emergent_clusters) <- rownames(sig_heatmap)
  
}

save(MCF7_emergent_clusters, file="processed_data/emergent_clusters_by_class/MCF7_emergent_clusters.RData")


##################################################################################

#Growing clusters
MCF7_growing_clusters <- list()

for(i in 1:nrow(heatmap)){
  
  curr_shrink <- (1:ncol(sig_heatmap))[sig_heatmap[i,] > 1.5]
  MCF7_growing_clusters <- append(MCF7_growing_clusters,list(curr_shrink))
  
}

names(MCF7_growing_clusters) <- rownames(heatmap)

saveRDS(MCF7_growing_clusters, "processed_data/cluster_info/MCF7_growing_clusters.rds")


#Shrinking clusters
MCF7_shrinking_clusters <- list()

for(i in 1:nrow(heatmap)){
  
  curr_shrink <- (1:ncol(sig_heatmap))[sig_heatmap[i,] < 1]
  MCF7_shrinking_clusters <- append(MCF7_shrinking_clusters,list(curr_shrink))
  
}

names(MCF7_shrinking_clusters) <- rownames(heatmap)

saveRDS(MCF7_shrinking_clusters, "processed_data/cluster_info/MCF7_shrinking_clusters.rds")



load("processed_data/emergent_clusters_list.Rdata")
