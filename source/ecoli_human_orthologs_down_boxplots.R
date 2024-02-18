setwd("/data/ruoffcj/projects/drug_treatment/")
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
source("source/cole_functions.R")

cell_lines <- c("A549","K562","MCF7")

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
  
  
  genesets_name <- "ecoli_human_orthologs_down"
  
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_",genesets_name,"_aucell_scores.rds"))
  
  type1_cell_names <- colnames(data)[data$cell_group == "1"]
  type2_cell_names <- colnames(data)[data$cell_group == "2"]
  non_rac_cell_names <- colnames(data)[data$cell_group == "0"]
  
  type1_cell_values <- scores[type1_cell_names,]
  type2_cell_values <- scores[type2_cell_names,]
  non_rac_values <- scores[non_rac_cell_names,]
  
  
  boxplot_df <- as.data.frame(cbind(scores,ifelse(rownames(scores) %in% non_rac_cell_names,"Non-RAC", ifelse(rownames(scores) %in% type1_cell_names, "RAC Type 1","RAC Type 2"))))
  
  colnames(boxplot_df) <- c("value","group")
  boxplot_df$value <- as.numeric(boxplot_df$value)
  boxplot_df$group <- factor(boxplot_df$group, levels = c("RAC Type 1","RAC Type 2","Non-RAC"))
  
  my_comparisons <- list( c("Non-RAC", "RAC Type 2"))
  
  p <- ggboxplot(boxplot_df, x = "group", y = "value", fill="group",short.panel.labs = FALSE)
  
  # Use only p.format as label. Remove method name.
  
  p <- p + stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox", size=8)+
    ggtitle(curr_cell_line)+
    xlab("")+
    ylab("")+
    scale_fill_manual(values=c("red", "orange","lightblue"),name = "Cell Groups")+
    theme(legend.position="right",
          title = element_text(size=20, face = "bold"),
          axis.text = element_text(size=20),
          legend.text = element_text(size=24),
          legend.title = element_text(size=26),
          legend.key.height = unit(1.5,"cm"),
          legend.key.width = unit(1.5,"cm"))
  
  plots <- append(plots,list(p))
  
}


png(paste0("/data/ruoffcj/projects/drug_treatment/figures/ecoli_human_orthologs_down_boxplots.png"),
    width=30, height=12, units="in",res = 300)


figure <- ggarrange(plotlist = plots, ncol=3, common.legend = T, legend=c("right"))

p <- annotate_figure(figure, left = text_grob("AUCell Score", rot = 90, vjust = 1, size=35, face="bold"),
                     bottom = text_grob("", size=35, face="bold"),
                     top=text_grob("E. Coli Antimicrobial Resistance Orthologs Cluster AUCell Scores", size=40, face="bold"))


print(p)

dev.off()

