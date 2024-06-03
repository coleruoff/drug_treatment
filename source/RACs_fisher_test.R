

curr_cell_line <- "MCF7"

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

data <- data[,data$treatment_stage=="pre"]

scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))

active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
active_cell_names <- active_cell_names[active_cell_names %in% colnames(data)]

all_clusters <- as.numeric(levels(data))


total_num_active <- length(active_cell_names)
total_num_inactive <- ncol(data)-total_num_active

ORs <- c()
pvals <- c()

for(i in all_clusters){
  curr_cluster_names <- colnames(data)[data$Cluster == i]
  
  curr_num_active <- sum(curr_cluster_names %in% active_cell_names)
  curr_num_inactive <- length(curr_cluster_names) - curr_num_active
  
  
  contin_table <- matrix(c(curr_num_active,total_num_active,curr_num_inactive,total_num_inactive), ncol=2)
  
  res <- fisher.test(contin_table)
  
  ORs <- append(ORs, res$estimate)
  pvals <- append(pvals, res$p.value)
  
}

plot(ORs)

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

all_clusters[ORs > 1.5]
