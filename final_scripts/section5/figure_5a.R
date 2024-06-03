args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(msigdbr)
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
source("final_scripts/drug_treatment_functions.R")
set.seed(42)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
plotDirectory <- "/data/ruoffcj/projects/drug_treatment/final_figures/"

################################################################################
# all_data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/all_cell_lines_data.rds"))

cell_lines <- c("A549","K562","MCF7")

supercluster_components <- readRDS(paste0(dataDirectory, "processed_data/supercluster_components.rds"))

organism_to_use <- "yeast"

df <- list()

for(curr_cell_line in cell_lines){
  data <- all_data[[curr_cell_line]]
  
  data <- AddMetaData(data, ifelse(data$Cluster == supercluster_components[[1]][[curr_cell_line]], 1,0), "supercluster1")
  data <- AddMetaData(data, ifelse(data$Cluster == supercluster_components[[2]][[curr_cell_line]], 1,0), "supercluster2")
  
  scores <- readRDS(paste0(dataDirectory, "aucell_score_objects/",curr_cell_line,"_processed_filtered_",organism_to_use,"_human_orthologs_up_aucell_scores.rds"))
  # scores <- scale(scores)
  data <- AddMetaData(data, scores)
  
  # Supercluster 1
  curr_scores <- data@meta.data %>% 
    filter(supercluster1 == 1) %>% 
    pull(paste0(organism_to_use,"_human_orthologs_up"))
  
  df[["scores"]] <- append(df[["scores"]],curr_scores)
  df[["sc"]] <- append(df[["sc"]],rep("Supercluster 1", length(curr_scores)))
  df[["cell_line"]] <- append(df[["cell_line"]],rep(curr_cell_line, length(curr_scores)))
  
  # Supercluster 2
  curr_scores <- data@meta.data %>% 
    filter(supercluster2 == 1) %>% 
    pull(paste0(organism_to_use,"_human_orthologs_up"))
  
  df[["scores"]] <- append(df[["scores"]],curr_scores)
  df[["sc"]] <- append(df[["sc"]],rep("Supercluster 2", length(curr_scores)))
  df[["cell_line"]] <- append(df[["cell_line"]],rep(curr_cell_line, length(curr_scores)))
  
  # Non-RACs
  curr_scores <- data@meta.data %>% 
    filter(rac == "nonrac") %>% 
    pull(paste0(organism_to_use,"_human_orthologs_up"))
  
  df[["scores"]] <- append(df[["scores"]],curr_scores)
  df[["sc"]] <- append(df[["sc"]],rep("Non-RAC", length(curr_scores)))
  df[["cell_line"]] <- append(df[["cell_line"]],rep(curr_cell_line, length(curr_scores)))
  
}

df <- data.frame(df)
df$rac <- ifelse(df$sc == "Non-RAC", "Non-RAC", "RAC")

# df <- df %>% 
#   filter(scores > 0)

my_comparisons <- list(c("Supercluster 1","Non-RAC"),c("Supercluster 2","Non-RAC"))

# p <- ggboxplot(df, x="sc",y="scores",fill="sc", facet.by = "cell_line")
p <- ggboxplot(df, x="sc",y="scores",fill="sc")

plot_title <- ifelse(organism_to_use == "ecoli", "E. coli Antimicrobial Resistance Human Ortholog Scores","Yeast Antifungal Resistance Human Ortholog Scores")

p <- p + stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "wilcox", label.x = 2.2, size=8, method.args = list(alternative = "greater"))+
  ggtitle(plot_title)+
  xlab("")+
  ylab("AUCell Score")+
  theme(legend.position="right",
        title = element_text(size=40, face = "bold"),
        axis.text = element_text(size=25),
        axis.text.x = element_text(size=30),
        axis.title = element_text(size=28),
        legend.text = element_text(size=24),
        legend.title = element_text(size=26),
        legend.key.height = unit(1.5,"cm"),
        legend.key.width = unit(1.5,"cm"))+
  NoLegend()


png(paste0(plotDirectory,"figure_5a.png"),
    width=20, height=12, units="in",res = 300)

print(p)

dev.off()

