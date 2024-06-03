

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

doses <- sort(unique(data$dose))
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


total_pre_cells <- sum(data$pathway_level_1=="Vehicle")       


df <- list()

for(curr_drug_class in drug_classes){
  cat(curr_drug_class, "\n")
  
  for(curr_dose in doses){
    
    if(curr_dose == 0){
      dose_total <- sum(data$pathway_level_1=="Vehicle" & data$dose == curr_dose)
    } else{
      dose_total <- sum(data$pathway_level_1==curr_drug_class & data$dose == curr_dose)
    }
    
    
    for(curr_cluster in clusters){
      
      if(curr_dose == 0){
        curr_num <- sum(data$pathway_level_1=="Vehicle" & data$Cluster == curr_cluster & data$dose == curr_dose)
      } else{
        curr_num <- sum(data$pathway_level_1==curr_drug_class & data$Cluster == curr_cluster & data$dose == curr_dose)
      }
      
      
      value <- curr_num/dose_total
        
      df[["drug_class"]] <- append(df[["drug_class"]], curr_drug_class)
      df[["cluster"]] <- append(df[["cluster"]], curr_cluster)
      df[["dose"]] <- append(df[["dose"]], curr_dose)
      df[["value"]] <- append(df[["value"]], value)
      
    }
  }
}

df <- data.frame(df)

df$dose <- factor(df$dose)



ggplot(df %>% 
         filter(cluster %in% RACs[[curr_cell_line]]))+
  geom_point(aes(x=dose,y=value))+
  geom_line(aes(x=dose,y=value,group=1))+
  facet_grid(cluster~drug_class)+
  ggtitle(curr_cell_line)




