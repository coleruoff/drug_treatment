setwd("/data/ruoffcj/projects/drug_treatment/")
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
source("source/cole_functions.R")

cell_lines <- c("A549","K562","MCF7")


curr_cell_line <- cell_lines[1]
plots <- list()
for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line,"\n")
  
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
  threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
  active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
  
  RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
  names(RACs) <- c("A549","K562","MCF7")
  clusters_of_interest <- RACs[[curr_cell_line]]
  
  data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, data$Cluster, "Non-RAC"), col.name = "RAC")
  #Add metadata for RAC and Cell Group
  data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, "1", ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, paste0(data$Cluster, "_1"), ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), paste0(data$Cluster, "_2"), paste0(data$Cluster, "_0"))), col.name = "cell_cluster_group")
  
  
  genesets_name <- "yeast_human_orthologs_up"
  
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_",genesets_name,"_aucell_scores.rds"))
  data <- AddMetaData(data, metadata=scores,col.name="yeast_score")
  
  boxplot_df <- data@meta.data %>% 
    dplyr::select(yeast_score,rac)
  
  colnames(boxplot_df) <- c("value","group")
  boxplot_df$value <- as.numeric(boxplot_df$value)

  boxplot_df$group <- ifelse(boxplot_df$group == "rac","RAC","Non-RAC")
  
  my_comparisons <- list(c("Non-RAC", "RAC"))

  
  p <- ggboxplot(boxplot_df, x = "group", y = "value", fill="group",short.panel.labs = FALSE)
    
  # Use only p.format as label. Remove method name.
  
  p <- p + stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox", size=8)+
    ggtitle(curr_cell_line)+
    xlab("")+
    ylab("")+
    scale_fill_manual(values=c("lightblue","orange"),name = "Cell Groups")+
    theme(legend.position="right",
          title = element_text(size=20, face = "bold"),
          axis.text = element_text(size=20),
          legend.text = element_text(size=24),
          legend.title = element_text(size=26),
          legend.key.height = unit(1.5,"cm"),
          legend.key.width = unit(1.5,"cm"))
  
  plots <- append(plots,list(p))
  
}


png(paste0("/data/ruoffcj/projects/drug_treatment/final_figures/figure_5a.png"),
    width=30, height=12, units="in",res = 300)


figure <- ggarrange(plotlist = plots, ncol=3, common.legend = T, legend=c("right"))

p <- annotate_figure(figure, left = text_grob("AUCell Score", rot = 90, vjust = 1, size=35, face="bold"),
                     bottom = text_grob("", size=35, face="bold"),
                     top=text_grob("Yeast Antifungal Resistance Orthologs AUCell Scores", size=40, face="bold"))


print(p)

dev.off()

