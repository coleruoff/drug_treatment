library(Seurat)
library(ggpubr)
library(tidyverse)
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
# plotDirectory <- "/data/ruoffcj/projects/drug_treatment/"
args = commandArgs(trailingOnly=TRUE)
dataDirectory <- args[1]
plotDirectory <- args[2]

###############################################################################

all_data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/all_cell_lines_data.rds"))
 
cell_lines <- c("A549","K562","MCF7")

supercluster1_components <- c(9,5,8)
supercluster2_components <- c(19,11,5)
supercluster3_components <- c(14,9,13)

names(supercluster1_components) <- cell_lines
names(supercluster2_components) <- cell_lines
names(supercluster3_components) <- cell_lines

df <- list()

genesets_to_use <- c("KEGG_GLUTATHIONE_METABOLISM", "HALLMARK_FATTY_ACID_METABOLISM", "cycling_signature")

for(curr_geneset in genesets_to_use){
  for(curr_cell_line in cell_lines){
    data <- all_data[[curr_cell_line]]
    
    data$supercluster <- ifelse(data$Cluster == supercluster1_components[curr_cell_line], 1, 0)
    data$supercluster <- ifelse(data$Cluster == supercluster2_components[curr_cell_line], 2, data$supercluster)
    data$supercluster <- ifelse(data$Cluster == supercluster3_components[curr_cell_line], 3, data$supercluster)
    
    scores1 <- readRDS(paste0(dataDirectory, "aucell_score_objects/",curr_cell_line,"_processed_filtered_oren_2021_gene_signatures_aucell_scores.rds"))
    scores2 <- readRDS(paste0(dataDirectory, "aucell_score_objects/",curr_cell_line,"_processed_filtered_hallmarks_aucell_scores.rds"))
    scores3 <- readRDS(paste0(dataDirectory, "aucell_score_objects/",curr_cell_line,"_processed_filtered_consensus_cycling_signature_aucell_scores.rds"))
    scores <- cbind(scores1,scores2,scores3)
    scores <- scale(scores)
    data <- AddMetaData(data, scores)
    
    for(i in 0:3){
      df[["scores"]] <- append(df[["scores"]], data@meta.data %>% 
                                 filter(supercluster == i) %>% 
                                 pull(curr_geneset))
      
      cells_count <- data@meta.data %>% 
        count(supercluster) %>% 
        filter(supercluster == i) %>% 
        pull(n)
      
      df[["sc"]] <- append(df[["sc"]], rep(paste0("Supercluster ",i), cells_count))
      
      df[["geneset"]] <- append(df[["geneset"]], rep(curr_geneset, cells_count))
    }
  }
}


df <- data.frame(df)

df$sc <- ifelse(df$sc == "Supercluster 0", "Non-RACs", df$sc)
df$sc <- factor(df$sc, levels = c("Supercluster 1","Supercluster 2","Supercluster 3","Non-RACs"))


my_comparisons <- list(c("Supercluster 1","Non-RACs"),c("Supercluster 2","Non-RACs"),c("Supercluster 3","Non-RACs"))

p <- ggboxplot(df, x="sc",y="scores",fill="sc", facet.by = "geneset")

plot_title <- paste0("Supercluster Scores")

p <- p + stat_compare_means(comparisons = my_comparisons,label = "p.format", method = "wilcox", label.x = 2.2, size=8)+
  ggtitle(plot_title)+
  xlab("")+
  ylab("AUCell Score (z-score)")+
  theme(legend.position="right",
        title = element_text(size=20, face = "bold"),
        axis.text = element_text(size=20),
        legend.text = element_text(size=24),
        legend.title = element_text(size=26),
        legend.key.height = unit(1.5,"cm"),
        legend.key.width = unit(1.5,"cm"))+
  NoLegend()

png(paste0(plotDirectory, "final_figures/figure_3d.png"),
    width = 26,height=12, units = 'in',res = 300)

print(p)

dev.off()

