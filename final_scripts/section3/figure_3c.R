args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(tidyverse)
library(Seurat)
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
# plotDirectory <- "/data/ruoffcj/projects/drug_treatment/"

#################################################################################

all_data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/all_cell_lines_data.rds"))

cell_lines <- c("A549","K562","MCF7")

supercluster_components <- readRDS(paste0(dataDirectory, "processed_data/supercluster_components.rds"))

total_sc1_counts <- list()
total_sc2_counts <- list()

for(curr_cell_line in cell_lines){
  data <- all_data[[curr_cell_line]]
  
  data <- AddMetaData(data, ifelse(data$Cluster == supercluster_components[[1]][[curr_cell_line]], 1,0), "supercluster1")
  data <- AddMetaData(data, ifelse(data$Cluster == supercluster_components[[2]][[curr_cell_line]], 1,0), "supercluster2")
  
  df <- data@meta.data
  
  curr_sc1 <- df %>% 
    dplyr::count(supercluster1, Phase) %>% 
    filter(supercluster1 == 1)
  
  total_sc1_counts[[curr_cell_line]] <- df %>% 
    filter(supercluster1 == 1) %>% 
    dplyr::select(cell, Phase)
  
  curr_sc2 <- df %>% 
    dplyr::count(supercluster2, Phase) %>% 
    filter(supercluster2 == 1)
  
  total_sc2_counts[[curr_cell_line]] <- df %>% 
    filter(supercluster2 == 1) %>% 
    dplyr::select(cell, Phase)
  
  
  
}


total_sc1_counts <- do.call("rbind",total_sc1_counts)
total_sc1_counts <- total_sc1_counts %>% 
  dplyr::count(Phase)
total_sc1_counts$n <- total_sc1_counts$n/sum(total_sc1_counts$n)

total_sc2_counts <- do.call("rbind",total_sc2_counts)
total_sc2_counts <- total_sc2_counts %>% 
  dplyr::count(Phase)
total_sc2_counts$n <- total_sc2_counts$n/sum(total_sc2_counts$n)

total_sc1_counts$supercluster <- "Supercluster 1"
total_sc2_counts$supercluster <- "Supercluster 2"

final_df <- rbind(total_sc1_counts,total_sc2_counts)

p <- ggplot(final_df)+
  geom_col(aes(x=Phase,y=n,fill=Phase), position = "dodge")+
  facet_wrap(~supercluster)+
  ggtitle("Superclusters Cell Cycle Phase Percentages")+
  xlab("Phase")+
  ylab("Percentage")


png(paste0(plotDirectory, "figure_3c.png"),
    width = 10,height=5, units = 'in',res = 300)

print(p)

dev.off()







