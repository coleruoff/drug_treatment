args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

library(tidyverse)
library(Seurat)
library(ggpubr)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/"

################################################################################

cell_lines <- c("A549","K562","MCF7")

curr_cell_line <- "A549"

plots <- list()
all_RACs <- list()
for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line,"\n")
  
  data <- readRDS(paste0(dataDirectory,"processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  data <- data[,data$treatment_stage=="post"]
  
  scores <- readRDS(paste0(dataDirectory, "aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
  threshold <- readRDS(paste0(dataDirectory, "aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
  
  active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
  active_cell_names <- active_cell_names[active_cell_names %in% colnames(data)]
  
  all_clusters <- as.numeric(levels(data))
  
  total_num_active <- length(active_cell_names)
  total_num_inactive <- ncol(data)-total_num_active
  
  ORs <- c()
  pvals <- c()
  
  for(i in all_clusters){
    curr_cluster_names <- colnames(data)[data$Cluster == i]
    
    curr_num_active <- sum(curr_cluster_names %in% active_cell_names)
    curr_num_inactive <- length(curr_cluster_names) - curr_num_active
    
    
    contin_table <- matrix(c(curr_num_active,total_num_active,curr_num_inactive,total_num_inactive), ncol=2)
    
    res <- fisher.test(contin_table)
    
    ORs <- append(ORs, res$estimate)
    pvals <- append(pvals, res$p.value)
  }
  
  df <- data.frame(cbind(paste0("", all_clusters),ORs))
  colnames(df) <- c("cluster","or")
  df$cluster <- factor(df$cluster, levels = all_clusters)
  df$or <- as.numeric(df$or)
  
  df$color <- ifelse(df$or > 1.5 & pvals < 0.05, "RAC","Non-RAC")
  
  curr_RACs <- df %>% 
    filter(color == "RAC") %>% 
    pull(cluster) %>% 
    as.numeric()
  
  all_RACs[[curr_cell_line]] <- curr_RACs
  
  plot_title <- curr_cell_line
  
  p <- ggplot(df)+
    geom_col(aes(x=cluster, y=or, fill=color))+
    scale_fill_manual(name="Cluster Type",values=c("lightblue","orange"))+
    geom_hline(yintercept=1.5, linetype="dashed", color = "red", linewidth=1)+
    xlab("")+
    ylab("")+
    ggtitle(plot_title)+
    theme(legend.position="right",
          axis.text = element_text(size=10),
          legend.text = element_text(size=10),
          legend.title = element_text(size=10),
          legend.key.height = unit(5,"mm"),
          legend.key.width = unit(5,"mm"))
  
  plots <- append(plots,list(p))
}

figure <- ggarrange(plotlist = plots, nrow=3, common.legend = T,legend=c("right"))

p <- annotate_figure(figure, left = text_grob("Odds Ratio", rot = 90, vjust = 1, size=20, face="bold"),
                bottom = text_grob("Clusters", size=20, face="bold"))

tiff(paste0(plotDirectory,"figure_1c.tiff"),width=200, height = 140, units = "mm", res = 1000)

print(p)

dev.off()

# png(paste0(plotDirectory,"figure_1c.png"),
#      width=20, height=20, units= "in", res=300)
# 
# print(p)
# 
# dev.off()

saveRDS(all_RACs, paste0(dataDirectory, "processed_data/all_RACs.rds"))


