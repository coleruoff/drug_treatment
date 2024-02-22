library(tidyverse)
library(ggpubr)

curr_cell_line <- "MCF7"
cat(curr_cell_line,"\n")

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- RACs[[curr_cell_line]]


scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
clusters_of_interest <- RACs[[curr_cell_line]]


data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, data$Cluster, "Non-RAC"), col.name = "RAC")
#Add metadata for RAC and Cell Group
data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, "1", ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, paste0(data$Cluster, "_1"), ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), paste0(data$Cluster, "_2"), paste0(data$Cluster, "_0"))), col.name = "cell_cluster_group")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


df <- data@meta.data


df$cell_group <- ifelse(df$cell_group == 0, "Non-RAC", paste0("RAC Type ", df$cell_group))

df$cell_group <- factor(df$cell_group, levels = c("RAC Type 1","RAC Type 2","Non-RAC"))

p <- ggboxplot(df, x = "Cluster", y = "G2M.Score",fill = "cell_group")

my_comparisons <- list(c("RAC Type 1", "RAC Type 2"),c("Non-RAC", "RAC Type 1"),c("Non-RAC", "RAC Type 2"))

plot_title <- paste0(curr_cell_line, " G2M Score")

p <- p + stat_compare_means(comparisons = my_comparisons,label = "p.format", method = "wilcox", label.x = 2.2, size=8)+
  ggtitle(plot_title)+
  xlab("")+
  ylab("")+
  scale_fill_manual(values=c("red","orange", "lightblue"),name = "Cell Groups")+
  theme(legend.position="right",
        title = element_text(size=20, face = "bold"),
        axis.text = element_text(size=20),
        legend.text = element_text(size=24),
        legend.title = element_text(size=26),
        legend.key.height = unit(1.5,"cm"),
        legend.key.width = unit(1.5,"cm"))



p

# Componenents of supercluster proliferation comparison
temp1 <- df %>% 
  filter(Cluster %in% c(8,12)) %>% 
  pull(proliferation_index)


temp2 <- df %>% 
  filter(Cluster %in% c(5)) %>% 
  pull(proliferation_index)

boxplot(temp1,temp2)
wilcox.test(temp1,temp2)

############
# Cell Cycle Phase percentages

df <- data@meta.data
temp <- df %>% 
  count(rac,Cluster, Phase)

percents <- df %>% 
  count(Cluster)

temp <- merge(temp,percents,by="Cluster")

temp <- temp %>% 
  mutate("percent" = n.x/n.y)

ggplot(temp)+
  geom_col(aes(x=Cluster,y=percent,fill=Phase), position = "dodge")+
  ggtitle(paste0(curr_cell_line, " Cell Cycle Phase Percentages"))+
  facet_grid(~rac)



############
#Endocytosis

scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line,"_processed_filtered_kegg_endocytosis_aucell_scores.rds"))
data <- AddMetaData(data, scores,col.name = "endocytosis_score")

df <- data@meta.data

df$cell_group <- ifelse(df$cell_group == 0, "Non-RAC", paste0("RAC Type ", df$cell_group))

df$cell_group <- factor(df$cell_group, levels = c("RAC Type 1","RAC Type 2","Non-RAC"))

p <- ggboxplot(df, x = "Cluster", y = "endocytosis_score",fill = "cell_group")

my_comparisons <- list(c("RAC Type 1", "RAC Type 2"),c("Non-RAC", "RAC Type 1"),c("Non-RAC", "RAC Type 2"))

plot_title <- paste0(curr_cell_line, " Endocytosis Score")

p <- p + stat_compare_means(comparisons = my_comparisons,label = "p.format", method = "wilcox", label.x = 2.2, size=8)+
  ggtitle(plot_title)+
  xlab("")+
  ylab("")+
  scale_fill_manual(values=c("red","orange", "lightblue"),name = "Cell Groups")+
  theme(legend.position="right",
        title = element_text(size=20, face = "bold"),
        axis.text = element_text(size=20),
        legend.text = element_text(size=24),
        legend.title = element_text(size=26),
        legend.key.height = unit(1.5,"cm"),
        legend.key.width = unit(1.5,"cm"))

p












