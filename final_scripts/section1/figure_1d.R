args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

library(tidyverse)
library(ggpubr)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"

################################################################################

cell_lines <- c("A549","K562","MCF7")

plots <- list()
for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line,"\n")
  
  data <- readRDS(paste0(dataDirectory,"processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  scores <- readRDS(paste0(dataDirectory, "aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
  threshold <- readRDS(paste0(dataDirectory, "aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
  
  active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
  active_cell_names <- active_cell_names[active_cell_names %in% colnames(data)]
  
  all_doses <- sort(unique(data$dose))
  all_clusters <- sort(unique(data$Cluster))
  
  ORs <- c()
  pvals <- c()
  
  total_df <- data.frame(matrix(NA,nrow=0,ncol=7))
  
  for(i in all_clusters){
    
    for(j in all_doses){
      curr_dose_names <- colnames(data)[data$dose == j]
      
      total_num_active <- sum(active_cell_names %in% curr_dose_names)
      total_num_inactive <- length(curr_dose_names)-total_num_active
      
      curr_cluster_names <- colnames(data)[data$Cluster == i & data$dose == j]
      
      curr_num_active <- sum(curr_cluster_names %in% active_cell_names)
      curr_num_inactive <- length(curr_cluster_names) - curr_num_active
      
      contin_table <- matrix(c(curr_num_active,total_num_active,curr_num_inactive,total_num_inactive), ncol=2)
      
      res <- fisher.test(contin_table)
      
      ORs <- append(ORs, res$estimate)
      pvals <- append(pvals, res$p.value)
      
      active_percent <- curr_num_active/(curr_num_inactive+curr_num_active)
      inactive_percent <- 1-active_percent
      
      cluster_size_percent <- length(curr_cluster_names)/length(curr_dose_names)
      
      prolif_index <- log(mean(data$proliferation_index[colnames(data) %in% curr_cluster_names])/mean(data$proliferation_index))
      
      
      total_df <- rbind(total_df,c(i,j,res$estimate,active_percent,inactive_percent,cluster_size_percent,prolif_index))
    }
  }
  
  colnames(total_df) <- c("cluster","dose","OR","percent","inactive_percent","cluster_size_percent","prolif_index")
  total_df$OR <- as.numeric(total_df$OR)
  total_df$percent <- as.numeric(total_df$percent)
  total_df$inactive_percent <- as.numeric(total_df$inactive_percent)
  total_df$cluster_size_percent <- as.numeric(total_df$cluster_size_percent)
  total_df$prolif_index <- as.numeric(total_df$prolif_index)
  
  total_df$cluster <- paste0("Cluster ", total_df$cluster)
  all_clusters <- paste0("Cluster ", all_clusters)
  
  plot_title <- curr_cell_line
  p <- ggplot(total_df,aes(x=dose, y=OR, group=1))+
    geom_line()+
    geom_hline(yintercept=1, linetype="dashed", color = "red", linewidth=1)+
    geom_point()+
    xlab("")+
    ylab("")+
    ggtitle(plot_title)+
    facet_wrap(~factor(cluster, levels=paste0("Cluster ", 1:19)))+
    theme(axis.text.x = element_text(size=14,angle=45, hjust=1),
          legend.position="right",
          title = element_text(size=28, face="bold"),
          axis.text.y = element_text(size=18),
          legend.text = element_text(size=14),
          legend.title = element_text(size=12),
          strip.text = element_text(size=20))
  
  plots <- append(plots,list(p))
}

figure <- ggarrange(plotlist = plots, ncol=3, nrow=1,common.legend = F)

p <- annotate_figure(figure, left = text_grob("Odds Ratio", rot = 90, vjust = 1, size=35, face="bold"),
                     bottom = text_grob("Dose", size=35, face="bold"))

png(paste0(plotDirectory, "figure_1d.png"),
    width=30, height=12, units = "in", res=300)

print(p)

dev.off()
