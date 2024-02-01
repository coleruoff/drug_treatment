setwd("/data/ruoffcj/projects/drug_treatment/")
library(tidyverse)
library(Seurat)
library(ggpubr)

curr_cell_line <- "MCF7"
use_pre <- F

cell_lines <- c("A549","K562","MCF7")

plots <- list()
for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line,"\n")
  
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  if(use_pre){
    pre_clusters <- data@meta.data %>%
      count(treatment_stage,Cluster) %>%
      filter(treatment_stage == "pre" & n > 10) %>%
      pull(Cluster)
    
    data <- data[,data$treatment_stage=='pre' & data$Cluster %in% pre_clusters]
  }
  
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
  threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
  
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
  
  df$color <- ifelse(df$or > 1.5, "RAC","Non-RAC")
  
  # plot_title <- paste0("Resistant Active/Inactive Odds Ratio (", curr_cell_line, " Pre-Treatment)")
  plot_title <- curr_cell_line
  # 
  # if(use_pre){
  #   png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/active_inactive_or_figures/",curr_cell_line,"_barplots_pre.png"),
  #       width=1000, height = 500)
  # } else {
  #   png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/active_inactive_or_figures/",curr_cell_line,"_barplots.png"),
  #       width=1000, height = 500)
  # }
  
  
  p <- ggplot(df)+
    geom_col(aes(x=cluster, y=or, fill=color))+
    scale_fill_manual(name="Cluster Type",values=c("lightblue","orange"))+
    geom_hline(yintercept=1.5, linetype="dashed", color = "red", size=1)+
    xlab("")+
    ylab("")+
    ggtitle(plot_title)+
    theme(legend.position="right",
          title = element_text(size=28, face="bold"),
          axis.text = element_text(size=30),
          legend.text = element_text(size=24),
          legend.title = element_text(size=26),
          legend.key.height = unit(1.5,"cm"),
          legend.key.width = unit(1.5,"cm"))
  
  plots <- append(plots,list(p))
}


png(paste0("/data/ruoffcj/projects/drug_treatment/final_figures/figure_1c.png"),
     width=20, height=20, units= "in", res=300)

figure <- ggarrange(plotlist = plots, nrow=3, common.legend = T,legend=c("right"))

p <- annotate_figure(figure, left = text_grob("Odds Ratio", rot = 90, vjust = 1, size=35, face="bold"),
                bottom = text_grob("Clusters", size=35, face="bold"),
                top=text_grob("Resistant Active/Inactive Odds Ratio", size=40, face="bold"))

print(p)

dev.off()


