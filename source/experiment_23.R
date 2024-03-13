source("source/read_in_all_cell_lines.R")

cell_lines <- c("A549","K562","MCF7")

supercluster1_components <- c(9,5,8)
supercluster2_components <- c(19,11,5)
supercluster3_components <- c(14,9,13)

names(supercluster1_components) <- cell_lines
names(supercluster2_components) <- cell_lines
names(supercluster3_components) <- cell_lines


total_sc1_counts <- list()
total_sc2_counts <- list()
total_sc3_counts <- list()

for(curr_cell_line in cell_lines){
  data <- all_data[[curr_cell_line]]
  
  data <- AddMetaData(data, ifelse(data$Cluster == supercluster1_components[curr_cell_line], 1,0), "supercluster1")
  data <- AddMetaData(data, ifelse(data$Cluster == supercluster2_components[curr_cell_line], 1,0), "supercluster2")
  data <- AddMetaData(data, ifelse(data$Cluster == supercluster3_components[curr_cell_line], 1,0), "supercluster3")
  
  
  df <- data@meta.data
  
  curr_sc1 <- df %>% 
    count(supercluster1, Phase) %>% 
    filter(supercluster1 == 1)
  
  total_sc1_counts[[curr_cell_line]] <- df %>% 
    filter(supercluster1 == 1) %>% 
    dplyr::select(cell, Phase)
  
  curr_sc2 <- df %>% 
    count(supercluster2, Phase) %>% 
    filter(supercluster2 == 1)
  
  total_sc2_counts[[curr_cell_line]] <- df %>% 
    filter(supercluster2 == 1) %>% 
    dplyr::select(cell, Phase)
  
  curr_sc3 <- df %>% 
    count(supercluster3, Phase) %>% 
    filter(supercluster3 == 1)
  
  total_sc3_counts[[curr_cell_line]] <- df %>% 
    filter(supercluster3 == 1) %>% 
    dplyr::select(cell, Phase)
  
}


total_sc1_counts <- do.call("rbind",total_sc1_counts)
total_sc1_counts <- total_sc1_counts %>% 
  count(Phase)
total_sc1_counts$n <- total_sc1_counts$n/sum(total_sc1_counts$n)

total_sc2_counts <- do.call("rbind",total_sc2_counts)
total_sc2_counts <- total_sc2_counts %>% 
  count(Phase)
total_sc2_counts$n <- total_sc2_counts$n/sum(total_sc2_counts$n)

total_sc3_counts <- do.call("rbind",total_sc3_counts)
total_sc3_counts <- total_sc3_counts %>% 
  count(Phase)
total_sc3_counts$n <- total_sc3_counts$n/sum(total_sc3_counts$n)

total_sc1_counts$supercluster <- "Supercluster 1"
total_sc2_counts$supercluster <- "Supercluster 2"
total_sc3_counts$supercluster <- "Supercluster 3"

final_df <- rbind(total_sc1_counts,total_sc2_counts,total_sc3_counts)

ggplot(final_df)+
  geom_col(aes(x=Phase,y=n,fill=Phase), position = "dodge")+
  facet_wrap(~supercluster)+
  ggtitle("Superclusters Cell Cycle Phase Percentages")+
  xlab("Phase")+
  ylab("Percentage")


  
  
  
  
  
  
  
  