library(tidyverse)
library(reshape2)
library(ggpubr)
library(rstatix)

specific_MPs <- readRDS("/data/CDSL_hannenhalli/Cole/data/specific_MPs.rds")

cell_lines <- c("A549","K562","MCF7")
curr_cell_line <- "A549"
# 
# for(curr_cell_line in cell_lines){
#   
#   cat(curr_cell_line, "\n")
#   # Read in cell line data
#   trapnell_data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/drug_treatment/processed_data/sciPlex_data/", curr_cell_line, "_processed.rds"))
#   
#   trapnell_data <- trapnell_data[,trapnell_data$treatment.stage == 'post']
#   
#   # Read in cluster group info
#   emergent_clusters <- readRDS("/data/CDSL_hannenhalli/Cole/drug_treatment/processed_data/cluster_info/emergent_clusters.rds")
#   emergent_clusters <- emergent_clusters[[curr_cell_line]]
#   
#   growing_clusters <- readRDS("/data/CDSL_hannenhalli/Cole/drug_treatment/processed_data/cluster_info/growing_non_emergent_clusters.rds")
#   growing_clusters <- unique(unlist(growing_clusters[[curr_cell_line]]))
#   
#   shrinking_clusters <- readRDS("/data/CDSL_hannenhalli/Cole/drug_treatment/processed_data/cluster_info/shrinking_clusters.rds")
#   shrinking_clusters <- unique(unlist(shrinking_clusters[[curr_cell_line]]))
#   
#   # Read in ITH scores and remove specific MPs
#   trapnell_ITH_scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/AUCell_score_objects/", curr_cell_line, "_ITH_meta_programs_AUC.rds"))
#   trapnell_ITH_scores <- trapnell_ITH_scores[!rownames(trapnell_ITH_scores) %in% specific_MPs,colnames(trapnell_ITH_scores) %in% colnames(trapnell_data)]
#   
#   # Remove "dying" cells
#   # cell_names_to_keep <- trapnell_data@meta.data %>% 
#   #   filter(Cluster != 19) %>% 
#   #   pull(cell)
#   # 
#   # 
#   # data_to_use <- trapnell_data[,cell_names_to_keep]
#   # scores <- trapnell_ITH_scores[,cell_names_to_keep]
#   
#   data_to_use <- trapnell_data
#   scores <- trapnell_ITH_scores
#   
#   
#   df <- data.frame(t(scores)) %>% 
#     cbind(.,"cluster" = data_to_use$Cluster) %>% 
#     rownames_to_column() %>% 
#     pivot_longer(!c("rowname", "cluster")) %>% 
#     mutate("emergent" = factor(ifelse(cluster %in% emergent_clusters, "Emergent","Non-emergent"))) %>% 
#     mutate("growing" = factor(ifelse(cluster %in% growing_clusters, "Growing","Non-growing"))) %>% 
#     mutate("shrinking" = factor(ifelse(cluster %in% shrinking_clusters, "Shrinking","Non-shrinking"))) %>% 
#     mutate("cell_line" = rep(curr_cell_line, nrow(.)))
#   
#   colnames(df) <- c("cell_name", "cluster", "geneset","value", "emergent","growing","shrinking", "cell_line")
#   
#   df$geneset <- sapply(df$geneset, FUN =  function(x) gsub("\\."," ", x))
#   df$geneset <- sapply(df$geneset, FUN =  function(x) gsub("  "," ", x))
#   
#   
#   if(which(curr_cell_line == cell_lines) == 1){
#     total_df <- df
#   } else {
#     total_df <- rbind(total_df, df)
#   }
# }
# 
# 
# saveRDS(total_df, "/data/CDSL_hannenhalli/Cole/drug_treatment/processed_data/all_cell_lines_MP_values_tidy.rds")


df <- readRDS("/data/CDSL_hannenhalli/Cole/drug_treatment/processed_data/all_cell_lines_MP_values_tidy.rds")
df <- readRDS("//hpcdrive.nih.gov/CDSL_hannenhalli/Cole/drug_treatment/processed_data/all_cell_lines_MP_values_tidy.rds")
geneset_names <- unique(df$geneset)

resistance_MPs <- factor(c("MP1 Cell Cycle  G2 M", "MP2 Cell Cycle  G1 S", "MP3 Cell Cylce HMG rich", "MP5 Stress", "MP6 Hypoxia",
                    "MP12 EMT I", "MP13 EMT II", "MP14 EMT III", "MP19 Epithelial Senescence","MP20 MYC"), levels = c("MP1 Cell Cycle  G2 M", "MP2 Cell Cycle  G1 S", "MP3 Cell Cylce HMG rich", "MP5 Stress", "MP6 Hypoxia",
                                                                                                                      "MP12 EMT I", "MP13 EMT II", "MP14 EMT III", "MP19 Epithelial Senescence","MP20 MYC"))

i <- geneset_names[1]

for(i in geneset_names){

  curr_mp <- gsub(" ","_",i)

  file_name <- paste0("/data/CDSL_hannenhalli/Cole/drug_treatment/results/cluster_MP_boxplots/all_cell_lines_",curr_mp,"_boxplot_growing.png")
  png(filename = file_name, height=1000, width=1000)

  plot_title <- paste0(i, ": AUCell Scores Across Cell Lines")

  bxp <- ggplot(df[df$geneset == i,])+
    geom_boxplot(aes(x=cell_line, y=value, fill=emergent))+
    ggtitle(plot_title)+
    scale_fill_manual(values = c("Non-emergent"="gray", "Emergent"="red"))+
    labs(fill = NULL)+
    xlab("")+
    ylab("AUCell Score")+
    theme(text = element_text(size = 20))


  colnames(df)

  stat.test <- df %>% 
    filter(geneset == i) %>% 
    group_by(cell_line) %>%
    wilcox_test(value ~ emergent) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")


  stat.test

  stat.test <- stat.test %>%
    add_xy_position(x = "cell_line", dodge = 0.8)

  plot(bxp + stat_pvalue_manual(stat.test,  label = "p", tip.length = 0))

  dev.off()

}

plot_title <- "AUCell Scores of Meta-Programs Associated with Resistance"

resistance_MPs <- factor(c("MP1 Cell Cycle  G2 M", "MP2 Cell Cycle  G1 S", "MP3 Cell Cylce HMG rich", "MP5 Stress", "MP6 Hypoxia",
                           "MP12 EMT I", "MP13 EMT II", "MP14 EMT III", "MP19 Epithelial Senescence","MP20 MYC"), levels = c("MP1 Cell Cycle  G2 M", "MP2 Cell Cycle  G1 S", "MP3 Cell Cylce HMG rich", "MP5 Stress", "MP6 Hypoxia",
                                                                                                                             "MP12 EMT I", "MP13 EMT II", "MP14 EMT III", "MP19 Epithelial Senescence","MP20 MYC"))

MPs_to_use <- factor(c("MP1 Cell Cycle  G2 M", "MP2 Cell Cycle  G1 S","MP12 EMT I","MP19 Epithelial Senescence"), 
                     levels = c("MP1 Cell Cycle  G2 M", "MP2 Cell Cycle  G1 S", "MP12 EMT I",  "MP19 Epithelial Senescence"))

emergent <- c(14:19)

df <- cbind(df,ifelse(df$cluster %in% emergent, df$cluster, "Non-emergent"))


colnames(df)[9] <- "type"


bxp <- ggplot(df[df$geneset %in% MPs_to_use & df$cell_line == "A549",])+
  geom_boxplot(aes(x=type, y=value, fill=emergent))+
  ggtitle(plot_title)+
  scale_fill_manual(values = c("Non-emergent"="gray", "Emergent"="red"))+
  labs(fill = NULL)+
  xlab("")+
  ylab("AUCell Score")+
  theme(text = element_text(size = 20))+
  facet_wrap(~factor(geneset, levels = resistance_MPs))



stat.test <- df %>% 
  group_by(type) %>%
  wilcox_test(value ~ type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

stat.test <- df %>% 
  wilcox_test(value ~ type) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")


stat.test <- stat.test %>%
  add_xy_position(x = "type", dodge = 0.8)

plot(bxp + stat_pvalue_manual(stat.test,  label = "p", tip.length = 0))









# 
# 
# bxp <- ggplot(df[df$geneset == i,])+
#   geom_boxplot(aes(x=cell_line, y=value, fill=emergent))+
#   ggtitle(plot_title)+
#   scale_fill_manual(values = c("Non-emergent"="gray", "Emergent"="red"))+
#   labs(fill = NULL)+
#   xlab("")+
#   ylab("AUCell Score")+
#   theme(text = element_text(size = 20))
# 
# 
# colnames(df)
# 
# stat.test <- df %>%
#   group_by(cell_line) %>%
#   t_test(value ~ emergent) %>%
#   adjust_pvalue(method = "bonferroni") %>%
#   add_significance("p.adj")
# 
# 
# stat.test
# 
# stat.test <- stat.test %>%
#   add_xy_position(x = "cell_line", dodge = 0.8)
# 
# plot(bxp + stat_pvalue_manual(stat.test,  label = "p", tip.length = 0))
# 
# # Add 10% spaces between the p-value labels and the plot border
# bxp + stat_pvalue_manual(stat.test,  label = "p", tip.length = 0) +
#   scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# 
# 
# #################################################################################
# for(i in geneset_names){
#   
#   curr_mp <- gsub(" ","_",i)
#   
#   file_name <- paste0("/data/CDSL_hannenhalli/Cole/drug_treatment/results/cluster_MP_boxplots/", curr_cell_line, "_cluster_",curr_mp,"_boxplot_emergent_clusters.png")
#   # png(filename = file_name, height=1000, width=1000)
#   
#   plot_title <- paste0(curr_cell_line, " Clusters AUCell Scores: ", i)
#   
#   #Emergent
#   plot(ggplot(df[df$geneset == i,])+
#     geom_boxplot(aes(x=cluster, y=value, fill=emergent))+
#     ggtitle(plot_title)+
#     scale_fill_manual(values = c("Non-emergent"="gray", "Emergent"="red"))+
#     labs(fill = NULL)+
#     xlab("Cluster")+
#     ylab("AUCell Score")+ 
#       theme(text = element_text(size = 20)))
#   
#   dev.off()
#   
#   file_name <- paste0("/data/CDSL_hannenhalli/Cole/drug_treatment/results/cluster_MP_boxplots/", curr_cell_line, "_cluster_",curr_mp,"_boxplot_emergent.png")
#   # png(filename = file_name, height=1000, width=1000)
#   
#   vec1 <- df %>% 
#     filter(emergent == "Emergent" & geneset==i) %>% 
#     pull(value)
#   
#   vec2 <- df %>% 
#     filter(emergent == "Non-emergent" & geneset==i) %>% 
#     pull(value)
#   
#   wilcox_res <- wilcox.test(vec1,vec2)
#   pval <- sprintf(wilcox_res$p.value, fmt="%.5f")
#   
#   plot_title <- paste0(curr_cell_line, " Clusters AUCell Scores: ", i, " (p value: ",pval,")")
#   
#   
#   plot(ggplot(df[df$geneset == i,])+
#     geom_boxplot(aes(x=cell_line, y=value, fill=emergent))+
#     ggtitle(plot_title)+
#     scale_fill_manual(values = c("Non-emergent"="gray", "Emergent"="red"))+
#     labs(fill = NULL)+
#     xlab("")+
#     ylab("AUCell Score")+ 
#     theme(text = element_text(size = 20)))
#   
#   dev.off()
#   
#   #Growing
#   file_name <- paste0("/data/CDSL_hannenhalli/Cole/drug_treatment/results/cluster_MP_boxplots/", curr_cell_line, "_cluster_",curr_mp,"_boxplot_growing_clusters.png")
#   # png(filename = file_name, height=1000, width=1000)
#   
#   plot_title <- paste0(curr_cell_line, " Clusters AUCell Scores: ", i)
#   
#   plot(ggplot(df[df$geneset == i,])+
#     geom_boxplot(aes(x=cluster, y=value, fill=growing))+
#     ggtitle(plot_title)+
#     scale_fill_manual(values = c("Non-growing"="gray", "Growing"="red"))+
#     labs(fill = NULL)+
#     xlab("Cluster")+
#     ylab("AUCell Score")+ 
#     theme(text = element_text(size = 20)))
#   
#   dev.off()
#   
#   file_name <- paste0("/data/CDSL_hannenhalli/Cole/drug_treatment/results/cluster_MP_boxplots/", curr_cell_line, "_cluster_",curr_mp,"_boxplot_growing.png")
#   # png(filename = file_name, height=1000, width=1000)
#   
#   vec1 <- df %>% 
#     filter(growing == "Growing" & geneset==i) %>% 
#     pull(value)
#   
#   vec2 <- df %>% 
#     filter(growing == "Non-growing" & geneset==i) %>% 
#     pull(value)
#   
#   wilcox_res <- wilcox.test(vec1,vec2)
#   pval <- sprintf(wilcox_res$p.value, fmt="%.5f")
#   
#   plot_title <- paste0(curr_cell_line, " Clusters AUCell Scores: ", i, " (p value: ",pval,")")
#   
#   plot(ggplot(df[df$geneset == i,])+
#     geom_boxplot(aes(x=cell_line, y=value, fill=growing))+
#     ggtitle(plot_title)+
#     scale_fill_manual(values = c("Non-growing"="gray", "Growing"="red"))+
#     labs(fill = NULL)+
#     xlab("")+
#     ylab("AUCell Score")+ 
#     theme(text = element_text(size = 20)))
#   
#   dev.off()
#   
#   #Shrinking
#   file_name <- paste0("/data/CDSL_hannenhalli/Cole/drug_treatment/results/cluster_MP_boxplots/", curr_cell_line, "_cluster_",curr_mp,"_boxplot_shrinking_clusters.png")
#   # png(filename = file_name, height=1000, width=1000)
#   
#   plot_title <- paste0(curr_cell_line, " Clusters AUCell Scores: ", i)
#   
#   plot(ggplot(df[df$geneset == i,])+
#     geom_boxplot(aes(x=cluster, y=value, fill=shrinking))+
#     ggtitle(plot_title)+
#     scale_fill_manual(values = c("Non-shrinking"="gray", "Shrinking"="red"))+
#     labs(fill = NULL)+
#     xlab("Cluster")+
#     ylab("AUCell Score")+ 
#     theme(text = element_text(size = 20)))
#   
#   dev.off()
#   
#   file_name <- paste0("/data/CDSL_hannenhalli/Cole/drug_treatment/results/cluster_MP_boxplots/", curr_cell_line, "_cluster_",curr_mp,"_boxplot_shrinking.png")
#   # png(filename = file_name, height=1000, width=1000)
#   
#   vec1 <- df %>% 
#     filter(shrinking == "Shrinking" & geneset==i) %>% 
#     pull(value)
#   
#   vec2 <- df %>% 
#     filter(shrinking == "Non-shrinking" & geneset==i) %>% 
#     pull(value)
#   
#   wilcox_res <- wilcox.test(vec1,vec2)
#   pval <- sprintf(wilcox_res$p.value, fmt="%.5f")
#   
#   plot_title <- paste0(curr_cell_line, " Clusters AUCell Scores: ", i, " (p value: ",pval,")")
#   
#   plot(ggplot(df[df$geneset == i,])+
#     geom_boxplot(aes(x=cell_line, y=value, fill=shrinking))+
#     ggtitle(plot_title)+
#     scale_fill_manual(values = c("Non-shrinking"="gray", "Shrinking"="red"))+
#     labs(fill = NULL)+
#     xlab("")+
#     ylab("AUCell Score")+ 
#     theme(text = element_text(size = 20)))
#   
#   dev.off()
# }
# 
# 
# 
# 
# 
#   
# 





