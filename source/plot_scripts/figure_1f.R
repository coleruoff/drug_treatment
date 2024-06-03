setwd("/data/ruoffcj/projects/drug_treatment/")
source("source/read_in_all_cell_lines.R")
library(ggpubr)

cell_lines <- c("A549","K562","MCF7")

geneset_to_use <- "hallmarks"
clean_geneset_name <- "Cancer Hallmarks"

plots <- list()
for(curr_cell_line in cell_lines){
  data <- all_data[[curr_cell_line]]
  
  rac_cell_names <- data@meta.data %>% 
    filter(rac == "rac") %>% 
    pull(cell)
  
  nonrac_cell_names <- data@meta.data %>% 
    filter(rac != "rac") %>% 
    pull(cell)
  
  trapnell_scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line,"_processed_filtered_",geneset_to_use,"_aucell_scores.rds"))
  raj_scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/raj_resistant_breast_processed_",geneset_to_use,"_aucell_scores.rds"))
  all_scores <- rbind(trapnell_scores,raj_scores)
  
  dist_mat <- as.matrix(dist(all_scores, method="euclidean"))
  
  rac_dist <- dist_mat[rownames(dist_mat) %in% rac_cell_names,!(colnames(dist_mat) %in% rac_cell_names | colnames(dist_mat) %in% nonrac_cell_names)]
  nonrac_dist <- dist_mat[rownames(dist_mat) %in% nonrac_cell_names,!(colnames(dist_mat) %in% rac_cell_names | colnames(dist_mat) %in% nonrac_cell_names)]
  
  df <- list()
  df[["group"]] <- c(rep("rac",nrow(rac_dist)),rep("nonrac",nrow(nonrac_dist)))
  df[["value"]] <- c(as.vector(rac_dist),as.vector(nonrac_dist))
  df <- data.frame(df)
  
  df$group <- ifelse(df$group == "nonrac", "Non-RAC", "RAC")
  df$group <- factor(df$group, levels = c("RAC","Non-RAC"))
  
  my_comparisons <- list(c("Non-RAC", "RAC"))
  p <- ggboxplot(df, x="group",y="value",fill="group")
  
  plot_title <- curr_cell_line
  p <- p + stat_compare_means(comparisons = my_comparisons,label = "p.format", method = "wilcox", label.x = 2.2, size=8)+
    ggtitle(plot_title)+
    xlab("")+
    ylab("")+
    scale_fill_manual(values=c("orange", "lightblue"),name = "Cell Groups")+
    theme(legend.position="right",
          title = element_text(size=20, face = "bold"),
          axis.text = element_text(size=20),
          legend.text = element_text(size=24),
          legend.title = element_text(size=26),
          legend.key.height = unit(1.5,"cm"),
          legend.key.width = unit(1.5,"cm"))
  
  
  plots <- append(plots, list(p))
}


png(paste0("/data/ruoffcj/projects/drug_treatment/final_figures/figure_1f.png"),
    width=30, height=12, units="in",res = 300)


figure <- ggarrange(plotlist = plots, ncol=3, common.legend = T, legend=c("right"))

p <- annotate_figure(figure, left = text_grob("Distance", rot = 90, vjust = 1, size=35, face="bold"),
                     bottom = text_grob("", size=35, face="bold"),
                     top=text_grob(paste0("Cell Group Distances to Resistant Clusters in ",clean_geneset_name, " Expression Space"), size=40, face="bold"))


print(p)

dev.off()
