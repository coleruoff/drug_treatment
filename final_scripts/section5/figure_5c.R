args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(msigdbr)
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
library(tidyverse)
source("final_scripts/drug_treatment_functions.R")
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
# plotDirectory <- "/data/ruoffcj/projects/drug_treatment/final_figures/"

################################################################################
all_data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/all_cell_lines_data.rds"))

cell_lines <- c("A549","K562","MCF7")

supercluster_components <- readRDS(paste0(dataDirectory, "processed_data/supercluster_components.rds"))

organism_to_use <- "ecoli"

df <- list()

for(curr_cell_line in cell_lines){
  data <- all_data[[curr_cell_line]]
  
  data <- AddMetaData(data, ifelse(data$Cluster == supercluster_components[[1]][[curr_cell_line]], 1,0), col.name = "supercluster1")
  data <- AddMetaData(data, ifelse(data$Cluster == supercluster_components[[2]][[curr_cell_line]], 1,0), col.name = "supercluster2")
  
  scores <- readRDS(paste0(dataDirectory, "aucell_score_objects/",curr_cell_line,"_processed_filtered_",organism_to_use,"_human_orthologs_up_aucell_scores.rds"))
  # scores <- scale(scores)
  data <- AddMetaData(data, scores[,1], col.name = colnames(scores))
  
  # Supercluster 1
  curr_scores <- data@meta.data %>% 
    filter(supercluster1 == 1 & treatment_stage == "post") %>% 
    pull(paste0(organism_to_use,"_human_orthologs_up"))
  
  df[["scores"]] <- append(df[["scores"]],curr_scores)
  df[["sc"]] <- append(df[["sc"]],rep("Supercluster 1", length(curr_scores)))
  df[["cell_line"]] <- append(df[["cell_line"]],rep(curr_cell_line, length(curr_scores)))
  
  # Supercluster 2
  curr_scores <- data@meta.data %>% 
    filter(supercluster2 == 1 & treatment_stage == "post") %>% 
    pull(paste0(organism_to_use,"_human_orthologs_up"))
  
  df[["scores"]] <- append(df[["scores"]],curr_scores)
  df[["sc"]] <- append(df[["sc"]],rep("Supercluster 2", length(curr_scores)))
  df[["cell_line"]] <- append(df[["cell_line"]],rep(curr_cell_line, length(curr_scores)))
  
  # Non-RACs
  curr_scores <- data@meta.data %>% 
    filter(rac == "nonrac" & treatment_stage == "post") %>% 
    pull(paste0(organism_to_use,"_human_orthologs_up"))
  
  df[["scores"]] <- append(df[["scores"]],curr_scores)
  df[["sc"]] <- append(df[["sc"]],rep("Non-RAC", length(curr_scores)))
  df[["cell_line"]] <- append(df[["cell_line"]],rep(curr_cell_line, length(curr_scores)))
  
}

df <- data.frame(df)
df$rac <- ifelse(df$sc == "Non-RAC", "Non-RAC", "RAC")


my_comparisons <- list(c("Supercluster 1","Non-RAC"),c("Supercluster 2","Non-RAC"))

# p <- ggboxplot(df, x="sc",y="scores",fill="sc", facet.by = "cell_line")
p <- ggboxplot(df, x="sc",y="scores",fill="sc")


stat.test <- compare_means(
  scores ~ sc, data = df,
  method = "wilcox.test",
  alternative = "less")

stat.test <- stat.test[-1,]
stat.test$p.format[1] <- "ns"

p <- p + stat_pvalue_manual(stat.test, label = "p.format", y.position = c(.12,.1), size=5)+
  xlab("")+
  ylab("AUCell Score")+
  ylim(0,.13)+
  theme(legend.position="right",
        axis.text = element_text(size=8),
        axis.text.x = element_text(size=10),
        axis.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
        legend.key.height = unit(1.5,"mm"),
        legend.key.width = unit(1.5,"mm"))+
  NoLegend()

p
 
# p <- p + stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "wilcox", label.x = 2.2, size=8, method.args = list(alternative = "greater"))+
#   xlab("")+
#   ylab("AUCell Score")+
#   theme(legend.position="right",
#         axis.text = element_text(size=25),
#         axis.text.x = element_text(size=30),
#         axis.title = element_text(size=28),
#         legend.text = element_text(size=24),
#         legend.title = element_text(size=26),
#         legend.key.height = unit(1.5,"cm"),
#         legend.key.width = unit(1.5,"cm"))+
#   NoLegend()

tiff(paste0(plotDirectory,"figure_5c.tiff"), width=100, height = 80, units = "mm", res = 1000)

print(p)

dev.off()

# png(paste0(plotDirectory, "figure_5c.png"),
#     width=12, height=8, units="in",res = 300)
# 
# print(p)
# 
# dev.off()
