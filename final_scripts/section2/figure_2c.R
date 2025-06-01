args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/"

#################################################################################

all_data <- readRDS(paste0(dataDirectory, "processed_data/sciplex_data/all_cell_lines_data.rds"))

cell_lines <- c("A549","K562","MCF7")

supercluster_components <- readRDS(paste0(dataDirectory, "processed_data/supercluster_components.rds"))


total_sc1_counts <- list()
total_sc2_counts <- list()
total_sc3_counts <- list()

for(curr_cell_line in cell_lines){
  cat(curr_cell_line, "\n")
  data <- all_data[[curr_cell_line]]
  
  data <- AddMetaData(data, ifelse(data$Cluster %in% supercluster_components[[1]][[curr_cell_line]], 1,0), "supercluster1")
  data <- AddMetaData(data, ifelse(data$Cluster %in% supercluster_components[[2]][[curr_cell_line]], 1,0), "supercluster2")
  data <- AddMetaData(data, ifelse(data$Cluster %in% supercluster_components[[3]][[curr_cell_line]], 1,0), "supercluster3")
  
  data <- data[,data$treatment_stage=="post"]
  
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

  curr_sc3 <- df %>%
    dplyr::count(supercluster3, Phase) %>%
    filter(supercluster3 == 1)

  total_sc3_counts[[curr_cell_line]] <- df %>%
    filter(supercluster3 == 1) %>%
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

total_sc3_counts <- do.call("rbind",total_sc3_counts)
total_sc3_counts <- total_sc3_counts %>%
  dplyr::count(Phase)
total_sc3_counts$n <- total_sc3_counts$n/sum(total_sc3_counts$n)


total_sc1_counts$supercluster <- "Supercluster 1"
total_sc2_counts$supercluster <- "Supercluster 2"
total_sc3_counts$supercluster <- "Supercluster 3"

final_df <- rbind(total_sc1_counts,total_sc2_counts,total_sc3_counts)

final_df$n <- final_df$n * 100

p <- ggplot(final_df)+
  geom_col(aes(x=Phase,y=n,fill=Phase), position = "dodge")+
  facet_wrap(~supercluster)+
  xlab("Phase")+
  ylab("Percentage")+
  theme_classic()+
  theme(strip.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 4),
        legend.key.height = unit(3,"mm"),
        legend.key.width = unit(3,"mm"),
        title = element_text(size=8),
        axis.line = element_line(linewidth=.2),
        axis.ticks = element_line(linewidth = .2))

p

jpeg(paste0(plotDirectory,"figure_2c.jpg"), width=100, height = 60, units = "mm", res = 600)
print(p)
dev.off()






