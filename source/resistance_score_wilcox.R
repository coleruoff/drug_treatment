library(tidyverse)
library(ggpubr)
library(rstatix)

curr_cell_line <- "MCF7"

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

# data <- data[,data$treatment_stage=="pre"]

scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))

active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]

all_clusters <- as.numeric(levels(data))

active_clusters <- c()

total_df <- matrix(NA, ncol=3,nrow=0)

for(i in all_clusters){
  cat(i,"\n")
  
  curr_cluster_names <- colnames(data)[data$Cluster == i & data$treatment_stage == "post"]
  
  if(length(curr_cluster_names) > 5){
    curr_cluster_scores <- scores[rownames(scores) %in% curr_cluster_names]
    rest_scores <- scores[rownames(scores) %in% colnames(data)]
    
    curr_df_segment <- cbind(rep("curr",length(curr_cluster_scores)),rep(i,length(curr_cluster_scores)), curr_cluster_scores)
    total_df <- rbind(total_df, curr_df_segment)
    
    sampled_rest <- sample(rest_scores,length(curr_cluster_scores))
    
    curr_df_segment <- cbind(rep("rest",length(sampled_rest)),rep(i,length(sampled_rest)), sampled_rest)
    total_df <- rbind(total_df, curr_df_segment)
    
    # boxplot(curr_cluster_scores,rest_scores)
    
    wilcox_res <- wilcox.test(curr_cluster_scores,rest_scores)
    
    if(wilcox_res$p.value < 0.05 & mean(curr_cluster_scores) > mean(rest_scores)){
      active_clusters <- append(active_clusters, i)
    }
  }
}


active_clusters

total_df <- as.data.frame(total_df)
colnames(total_df) <- c("group","cluster","value")
total_df$value <- as.numeric(total_df$value)
total_df$cluster <- factor(total_df$cluster, levels=all_clusters)

stat.test <- total_df %>%
  group_by(cluster) %>%
  wilcox_test(value ~ group) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")

stat.test

bxp <- ggboxplot(total_df, x = "group", y = "value", fill = "group",
                 facet.by = "cluster")

stat.test <- stat.test %>% add_xy_position(x = "group")
stat.test$custom.label <- ifelse(stat.test$p.adj <= 0.05, stat.test$p.adj, "ns")

# Visualization
bxp + 
  stat_pvalue_manual(stat.test, label = "{p}{p.adj.signif}") +
  ggtitle(curr_cell_line)



ggplot(total_df)+
  geom_boxplot(aes(x=group,y=value))+
  facet_wrap(~cluster)







cat(mean(curr_cluster_scores)/mean(rest_scores),"\n")



i <- 1
for(i in all_clusters){
  
  curr_cluster_pre_names <- colnames(data)[data$Cluster == i & data$treatment_stage =="pre"]
  curr_cluster_post_names <- colnames(data)[data$Cluster == i & data$treatment_stage =="post"]
  
  pre_scores <- scores[rownames(scores) %in% curr_cluster_pre_names]
  
  post_scores <- scores[rownames(scores) %in% curr_cluster_post_names]
  
  boxplot(pre_scores,post_scores)
  
  wilcox_res <- wilcox.test(pre_scores,post_scores)
  
  if(mean(post_scores) > mean(pre_scores)){
    active_clusters <- append(active_clusters, i)
  }
  
}

total_pre_names <- colnames(data)[data$treatment_stage =="pre"]
total_post_names <- colnames(data)[data$treatment_stage =="post"]

total_pre_scores <- scores[rownames(scores) %in% total_pre_names]
total_post_scores <- scores[rownames(scores) %in% total_post_names]

total_pre_mean <- mean(total_pre_scores)
total_post_mean <- mean(total_post_scores)

ORs <- c()
pvals <- c()

for(i in all_clusters){
  
  if(sum(data$Cluster == i & data$treatment_stage =="pre") > 5){
    curr_cluster_pre_names <- colnames(data)[data$Cluster == i & data$treatment_stage =="pre"]
    curr_cluster_post_names <- colnames(data)[data$Cluster == i & data$treatment_stage =="post"]
    
    pre_scores <- scores[rownames(scores) %in% curr_cluster_pre_names]
    
    post_scores <- scores[rownames(scores) %in% curr_cluster_post_names]
    
    
    contin_table <- matrix(c(mean(post_scores),total_post_mean,mean(pre_scores),total_pre_mean), ncol=2)
    contin_table <- contin_table*100000 
    
    res <- fisher.test(contin_table)
    
    ORs <- append(ORs, res$estimate)
    pvals <- append(pvals, res$p.value)
  }
}

df <- data.frame(cbind(paste0("Cluster ", all_clusters),ORs))
colnames(df) <- c("cluster","or")
df$cluster <- factor(df$cluster, levels = paste0("Cluster ", all_clusters))
df$or <- as.numeric(df$or)

df$color <- ifelse(df$or > 1, "RAC","Non-RAC")

plot_title <- paste0(curr_cell_line, " Resistant Active/Inactive Odds Ratio")



ggplot(df)+
  geom_col(aes(x=cluster, y=or, fill=color))+
  scale_fill_manual(name="",values=c("#56B4E9","#E69F00"))+
  geom_hline(yintercept=1, linetype="dashed", color = "red", size=1)+
  xlab("Cluster")+
  ylab("Odds Ratio")+
  ggtitle(plot_title)+
  theme(axis.text.x = element_text(angle=45, hjust=1))



dim(data)
