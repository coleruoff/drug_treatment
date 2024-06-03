source("source/cole_functions.R")

cell_lines <- c("A549","K562","MCF7")

#################################################################################
# Read in Trapnell data
#################################################################################
curr_cell_line <- cell_lines[1]

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

num_clusters <- nlevels(data$Cluster)

#################################################################################
# Read in AUCell score results and thresholds
#################################################################################

resistant_cancer_type <- "common"

scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line, "_processed_filtered_raj_",resistant_cancer_type,"_resistance_signature_aucell_scores.rds"))

threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line,"_processed_filtered_raj_",resistant_cancer_type,"_resistance_signature_aucell_thesholds.rds"))

all_active_cells <- rownames(scores)[scores > threshold$threshold]

all_active_cells <- rownames(scores)[scores > quantile(scores, probs = .9)]

#################################################################################
# Populate heatmap
#################################################################################

Idents(data) <- data$Cluster
cluster_dose_combos <- c()

for(i in 1:nlevels(data)){
  
  cluster_dose_combos <- append(cluster_dose_combos, paste(i,sort(unique(data$dose)), sep = "_"))
}


cluster_dose_heatmap <- matrix(NA, nrow=1, ncol=length(cluster_dose_combos))

heatmap_column_names <- c()

for(j in 1:length(cluster_dose_combos)){
  cat(j, "\n")
  
  
  split_string <- strsplit(cluster_dose_combos[j], "_")
  curr_cluster <- split_string[[1]][1]
  curr_dose <- split_string[[1]][2]
  
  heatmap_column_names <- append(heatmap_column_names, paste0("Cluster ", curr_cluster, ": Dose ", curr_dose))
  
  curr_cells_names <- rownames(data@meta.data)[data$Cluster == curr_cluster & data$dose == curr_dose]
  
  curr_total_cells <- length(curr_cells_names)
  
  if(curr_total_cells < 10){
    
    cluster_dose_heatmap[1, j] <- 0
    
  } else {
    curr_active_cells <- sum(curr_cells_names %in% all_active_cells)
    
    cluster_dose_heatmap[1, j] <- curr_active_cells/curr_total_cells
  }
}


colnames(cluster_dose_heatmap) <- heatmap_column_names



df <- matrix(NA, ncol=3,nrow=0)

for(i in 1:length(cluster_dose_combos)){
  
  df <- rbind(df,c(strsplit(cluster_dose_combos[i],"_")[[1]][1],strsplit(cluster_dose_combos[i],"_")[[1]][2],t(cluster_dose_heatmap)[i,]))
  
}

colnames(df) <- c("cluster","dose","value")
df <- data.frame(df)
df$value <- as.numeric(df$value)
df$cluster <- as.factor(df$cluster)




ggplot(df,aes(x=cluster,y=value, fill=dose))+
  geom_bar(stat="identity",position=position_dodge())

#################################################################################
# Plot heatmap
#################################################################################


plot_title <- paste0("Percentage of Raj ", str_to_title(resistant_cancer_type), " Resistant Signature Active Cells in ",curr_cell_line, " Clusters")

col_fun <-colorRamp2(c(0, max(cluster_dose_heatmap)), c("white", "red"))

Heatmap(cluster_dose_heatmap, name="Percentage", column_title = plot_title, cluster_rows = F, cluster_columns = F, 
        row_names_side = "left", col = col_fun, column_names_rot = 45,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", cluster_dose_heatmap[i, j]), x, y, gp = gpar(fontsize = 5))})


#################################################################################

#################################################################################









