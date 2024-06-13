args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(tidyverse)
library(ggpubr)
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
# plotDirectory <- "/data/ruoffcj/projects/drug_treatment/final_figures/"

##############################################################################

cell_lines <- c("A549","K562","MCF7")

RACs <- readRDS(paste0(dataDirectory, "processed_data/all_RACs.rds"))

supercluster_componenets <- readRDS(paste0(dataDirectory, "processed_data/supercluster_components.rds"))

all_data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/all_cell_lines_data.rds"))

curr_cell_line <- cell_lines[2]

plots <- list()
all_plot_data <- list()
for(curr_cell_line in cell_lines){
  data <- all_data[[curr_cell_line]]
  
  drug_classes <- as.character(unique(data$pathway_level_1))
  drug_classes <- drug_classes[-which(drug_classes == "Vehicle")]
  drug_classes <- drug_classes[-which(drug_classes == "Other")]
  drug_classes <- sort(drug_classes)
  
  # clusters <- as.numeric(levels(data))
  clusters <- c(supercluster_componenets[[1]][[curr_cell_line]],supercluster_componenets[[2]][[curr_cell_line]])
  
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
        
        curr_sc <- which(curr_cluster == clusters)
        
        df[["drug_class"]] <- append(df[["drug_class"]], curr_drug_class)
        df[["cluster"]] <- append(df[["cluster"]], curr_cell_line)
        df[["dose"]] <- append(df[["dose"]], curr_dose)
        df[["value"]] <- append(df[["value"]], value)
        df[["supercluster"]] <- append(df[["supercluster"]], curr_sc)
        
      }
    }
  }
  
  plot_df <- data.frame(df)
  
  all_plot_data <- append(all_plot_data, list(plot_df))
  
  # plot_df$dose <- factor(plot_df$dose)
  # 
  # plot_df <- plot_df %>% 
  #   filter(cluster %in% RACs[[curr_cell_line]])
  # 
  # cluster_order <- sort(unique(plot_df$cluster))
  # plot_df$cluster <- paste0("Cluster ", plot_df$cluster)
  # plot_df$cluster <- factor(plot_df$cluster, levels = paste0("Cluster ", cluster_order))
  # 
  # plot_title <- curr_cell_line
  # 
  # p <- ggplot(plot_df)+
  #   geom_point(aes(x=dose,y=value))+
  #   geom_line(aes(x=dose,y=value,group=1))+
  #   facet_grid(cluster~drug_class)+
  #   ggtitle(plot_title)+
  #   xlab("")+
  #   ylab("")+
  #   theme(plot.title = element_text(size=30, face="bold"),
  #         axis.text.x = element_text(size=10))
  # 
  # 
  # plots <- append(plots, list(p))
}



all_plot_data <- do.call(rbind, all_plot_data)

plots <- list()
for(curr_sc in 1:2){
  
  df <- all_plot_data %>% 
    filter(supercluster == curr_sc)
  
  plot_title <- paste0("Supercluster ", curr_sc)
  
  df$dose <- factor(df$dose)

  
  p <- ggplot(df)+
      geom_point(aes(x=dose,y=value))+
      geom_line(aes(x=dose,y=value,group=1))+
      facet_grid(drug_class~cluster)+
      ggtitle(plot_title)+
      xlab("")+
      ylab("")+
      theme(plot.title = element_text(size=30, face="bold"),
            axis.text.x = element_text(size=10),
            strip.text.x.top = element_text(size=20))
  
  
  if(curr_sc == 1){
    p <- p+theme(strip.text.y.right  = element_blank())
  } else{ 
    p <- p+theme(strip.text.y.right  = element_text(size=16, angle = 0))
  }
  
  plots <- append(plots, list(p))
  
}



figure <- ggarrange(plotlist = plots, widths = c(1,1.5),ncol=2, common.legend = T,legend=c("right"))

p <- annotate_figure(figure, left = text_grob("Percentage", rot = 90, vjust = 1, size=35, face="bold"),
                     bottom = text_grob("Dose", size=35, face="bold"))



p

png(paste0(plotDirectory, "figure_S2d.png"),
    width = 20,height=20, units = 'in',res = 300)

print(p)

dev.off()
