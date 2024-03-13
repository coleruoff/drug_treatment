source("source/read_in_all_cell_lines.R")
library(Seurat)
library(ggpubr)

###############################################################################
geneset_to_use <- "HALLMARK_NOTCH_SIGNALING"
colnames(scores)
cell_lines <- c("A549","K562","MCF7")

supercluster1_components <- c(9,5,8)
supercluster2_components <- c(19,11,5)
supercluster3_components <- c(14,9,13)

names(supercluster1_components) <- cell_lines
names(supercluster2_components) <- cell_lines
names(supercluster3_components) <- cell_lines

sc1_scores <- c()
sc2_scores <- c()
sc3_scores <- c()
rest_scores <- c()

for(curr_cell_line in cell_lines){
  data <- all_data[[curr_cell_line]]
  
  data <- AddMetaData(data, ifelse(data$Cluster == supercluster1_components[curr_cell_line], 1,0), "supercluster1")
  data <- AddMetaData(data, ifelse(data$Cluster == supercluster2_components[curr_cell_line], 1,0), "supercluster2")
  data <- AddMetaData(data, ifelse(data$Cluster == supercluster3_components[curr_cell_line], 1,0), "supercluster3")
  
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line,"_processed_filtered_oren_2021_gene_signatures_aucell_scores.rds"))
  scores <- scale(scores)
  data <- AddMetaData(data, scores)
  
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line,"_processed_filtered_hallmarks_aucell_scores.rds"))
  scores <- scale(scores)
  data <- AddMetaData(data, scores)
  
  df <- data@meta.data
  sc1_scores <- append(sc1_scores,df %>% 
                         filter(supercluster1 == 1) %>% 
                         pull(geneset_to_use))
  
  sc2_scores <- append(sc2_scores,df %>% 
                         filter(supercluster2 == 1) %>% 
                         pull(geneset_to_use))
  
  sc3_scores <- append(sc3_scores,df %>% 
                         filter(supercluster3 == 1) %>% 
                         pull(geneset_to_use))
  
  rest_scores <- append(rest_scores,df %>% 
                          filter(rac == "nonrac") %>%
                          pull(geneset_to_use))
  
}

df <- data.frame(rbind(cbind("Supercluster 1",sc1_scores),cbind("Supercluster 2",sc2_scores),cbind("Supercluster 3",sc3_scores),cbind("Non-RACs",rest_scores)))

colnames(df) <- c("sc","scores")
df$scores <- as.numeric(df$scores)

df <- df %>%
  filter(scores > 0)

my_comparisons <- list(c("Supercluster 1","Non-RACs"),c("Supercluster 2","Non-RACs"),c("Supercluster 3","Non-RACs"))

p <- ggboxplot(df, x="sc",y="scores",fill="sc")

plot_title <- paste0(geneset_to_use," Score")

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

p

