library(ComplexHeatmap)
library(circlize)

cell_lines <- c("A549","K562","MCF7")

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- cell_lines

pre_RACs <- list(c(4,9,12,13),c(4,5,11),c(5,8,12,17))
names(pre_RACs) <- cell_lines

curr_cell_line <- cell_lines[3]

data <- all_data[[curr_cell_line]]

drug_classes <- as.character(unique(data$pathway_level_1))
drug_classes <- drug_classes[-which(drug_classes == "Vehicle")]
drug_classes <- drug_classes[-which(drug_classes == "Other")]
drug_classes <- sort(drug_classes)
clusters <- as.numeric(levels(data))
#########################################################

pre_cluster_cell_counts <- c()
for(curr_cluster in clusters){
  cat(curr_cluster,"\n")
  
  curr_num <- sum(data$pathway_level_1=="Vehicle" & data$Cluster == curr_cluster)
  
  if(curr_num < 10){
    curr_num <- 0
  } 
  
  pre_cluster_cell_counts <- append(pre_cluster_cell_counts, curr_num)
}

# pre_cluster_cell_counts[pre_cluster_cell_counts < 10] <- 0

total_pre_cells <- sum(data$pathway_level_1=="Vehicle")       


heatmap <- matrix(NA, ncol=length(clusters),nrow=length(drug_classes))
colnames(heatmap) <- clusters
rownames(heatmap) <- drug_classes

curr_drug_class <- drug_classes[3]
curr_cluster <- 14

for(curr_drug_class in drug_classes){
  cat(curr_drug_class, "\n")
  
  total_curr_cells <- sum(data$pathway_level_1==curr_drug_class)           
  
  for(curr_cluster in clusters){
    
    curr_num <- sum(data$pathway_level_1==curr_drug_class & data$Cluster == curr_cluster)
    
    if(curr_num < 10){
      curr_num <- 0
    } 
    
    contin_table <- matrix(c(curr_num,pre_cluster_cell_counts[curr_cluster],total_curr_cells,total_pre_cells), ncol=2)
    
    if(pre_cluster_cell_counts[curr_cluster] == 0 & curr_num == 0){
      heatmap[curr_drug_class,curr_cluster] <- 1
      
    } else {
      res <- fisher.test(contin_table+1)
      if(res$p.value < 0.05){
        heatmap[curr_drug_class,curr_cluster] <- res$estimate
      } else{
        heatmap[curr_drug_class,curr_cluster] <- res$estimate
      }
    }
  }
}





# heatmap <- heatmap[,colnames(heatmap) %in% pre_RACs[[curr_cell_line]]]

heatmap <- heatmap[,colnames(heatmap) %in% RACs[[curr_cell_line]]]

heatmap <- ifelse(heatmap > 10, 10,heatmap)

col_fun = colorRamp2(c(-1, 0, max(log(heatmap))), c("blue", "white", "red"))
Heatmap(log(heatmap),cluster_rows = F,cluster_columns = F,column_title = curr_cell_line, col=col_fun)




