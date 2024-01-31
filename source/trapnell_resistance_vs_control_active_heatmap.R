source("source/cole_functions.R")

cell_lines <- c("A549","K562","MCF7")

A549.data <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/A549_processed_filtered.rds")

K562.data <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/K562_processed_filtered.rds")

MCF7.data <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/MCF7_processed_filtered.rds")

cell_lines_data <- c(A549.data,K562.data,MCF7.data)
#################################################################################
# Read in Trapnell data
#################################################################################
i <- 1
data <- cell_lines_data[[i]]
curr_cell_line <- cell_lines[[i]]

num_clusters <- nlevels(data$Cluster)

growing <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/trapnell_cluster_groups/growing_clusters.rds")
growing <- growing[[curr_cell_line]]


emergent <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/trapnell_cluster_groups/emergent_clusters.rds")
emergent <- emergent[[curr_cell_line]]

#################################################################################
# Read in AUCell score results and thresholds
#################################################################################

scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_resistance_signatures_aucell_scores.rds"))

thresholds <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line,"_processed_filtered_raj_resistance_signatures_aucell_thresholds.rds"))

#################################################################################
# Populate heatmap
#################################################################################

curr_geneset <- "raj_common_200"

raj_resistant_cancer_type <- strsplit(curr_geneset,"_")[[1]][[2]]

all_active_cells <- rownames(scores)[scores[,curr_geneset] > thresholds$threshold[thresholds$gene_set == curr_geneset]]

percentage_heatmap <- matrix(NA, nrow=2, ncol=(num_clusters+1))

for(curr_stage in c("pre","post")){
  cat(curr_stage, "\n")
  
  curr_row <- ifelse(curr_stage == "pre", 1, 2)
  
  for(j in 1:num_clusters){
    cat(j, "\n")
    
    curr_cells_names <- rownames(data@meta.data)[data$Cluster == j & data$treatment_stage == curr_stage]
    
    curr_total_cells <- length(curr_cells_names)
    
    if(curr_total_cells < 10){
      
      percentage_heatmap[curr_row,j] <- 0
      
    } else {
      curr_active_cells <- sum(curr_cells_names %in% all_active_cells)
      
      percentage_heatmap[curr_row,j]  <- curr_active_cells/curr_total_cells
    }
  }
}

#Add total pre percentage in last column
pre_cells_names <- rownames(data@meta.data)[data$treatment_stage == "pre"]
active_pre_cells <- sum(pre_cells_names %in% all_active_cells)

percentage_heatmap[1, ncol(percentage_heatmap)] <- active_pre_cells/length(pre_cells_names)

#Add total post percentage in last column
post_cells_names <- rownames(data@meta.data)[data$treatment_stage == "post"]
active_post_cells <- sum(post_cells_names %in% all_active_cells)

percentage_heatmap[2, ncol(percentage_heatmap)] <- active_post_cells/length(post_cells_names)


colnames(percentage_heatmap)[1:ncol(percentage_heatmap)] <- paste0("Cluster ", 1:ncol(percentage_heatmap))
colnames(percentage_heatmap)[ncol(percentage_heatmap)] <- "Total"
rownames(percentage_heatmap) <- c("Pre-Treatment","Post-Treatment")

#################################################################################
# Plot heatmap
#################################################################################

plot_title <- paste0("Percentage of Raj ", str_to_title(raj_resistant_cancer_type), " Resistant Signature Active Cells in ",curr_cell_line, " Clusters")

col_fun <-colorRamp2(c(0, 1), c("white", "red"))

ha <- HeatmapAnnotation(Emergent = c(ifelse(1:num_clusters %in% emergent,"Emergent","Non-emergent"),"NA"),
                        col = list(Emergent = c("Emergent" = "red", "Non-emergent" = "blue", "NA" = "black")))

ht <- Heatmap(percentage_heatmap, name="Percentage", column_title = plot_title, cluster_rows = F, cluster_columns = F, 
              row_names_side = "left", col = col_fun, column_names_rot = 45, column_title_gp = gpar(fontsize=22),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", percentage_heatmap[i, j]), x, y, gp = gpar(fontsize = 12))},
              bottom_annotation = ha)


draw(ht)




# Active Percentage Distributions in Emergent Clusters
df <- data.frame(cbind(percentage_heatmap[2,1:(ncol(percentage_heatmap)-1)],ifelse(1:num_clusters %in% emergent,"Emergent","Non-emergent")))
colnames(df) <- c("value","emergent")
df$value <- as.numeric(df$value)

vec1 <- as.vector(percentage_heatmap[2,emergent])
vec2 <- as.vector(percentage_heatmap[2,non_emergent])
wilcox_res <- wilcox.test(vec1, vec2)

pval <- wilcox_res$p.value

boxplot_title <- paste0("Percentage of Raj ", str_to_title(raj_resistant_cancer_type), " Resistant Signature Active Cells in ",curr_cell_line, " Clusters (p-value: ", sprintf("%.4f",pval),")")

ggplot(df)+
  geom_boxplot(aes(x=emergent,y=value, fill=emergent))+
  xlab("")+
  ylab("Percent of Active Cells")+
  ggtitle(boxplot_title)+
  NoLegend()


# 
percentages <- c()
for(j in 1:num_clusters){
  cat(j, "\n")
  
  curr_cells_names <- rownames(data@meta.data)[data$Cluster == j]
  
  curr_total_cells <- length(curr_cells_names)
  
  if(curr_total_cells < 10){
    
    percentages  <- append(percentages,0)
    
  } else {
    curr_active_cells <- sum(curr_cells_names %in% all_active_cells)
    
    percentages  <- append(percentages,curr_active_cells/curr_total_cells)
  }
}

plot_title <- paste0("Distribution of Active Percentages of ",curr_cell_line, " Clusters")

plot(density(percentages), main=plot_title, xlab="Percent of Active Cells", ylab="")+
  abline(v=.25, col="blue")




(1:nlevels(data))[percentages>.045]








