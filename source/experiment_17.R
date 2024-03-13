source("source/read_in_all_cell_lines.R")
source("/data/ruoffcj/projects/aucell_scoring/aucell_thresholding.R")
library(AUCell)
library(Seurat)

curr_cell_line <- cell_lines[1]

cat(curr_cell_line, "\n")

data <- all_data[[curr_cell_line]]  

###################################

pre_data <- data[,data$treatment_stage == "pre"]

pre_clusters <- pre_data@meta.data %>% 
  count(Cluster) %>% 
  filter(n > .01*ncol(pre_data)) %>% 
  pull(Cluster)

pre_data <- pre_data[,pre_data$Cluster %in% pre_clusters]


dim(pre_data)

resistance_signature <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_watermelon_resistance_signature.rds")

data_to_use <- pre_data
genesets <- resistance_signature

set.seed(42) #Make sure that a random seed is set to a fixed number. This ensures reproducibility across runs.
# auc_obj <- compute_AUCell_scores(data_to_use, genesets, compute_thresholds=F, nCores = 1, assay_to_use = "RNA")
# 
# 
# #Compute separate thresholds across time-points
# computed_thresholds_df <- compute_shuffled_gene_set_AUCell_scores(data_to_use, gene_sets=genesets, nCores=1, do_sample_wise=F, q_thresh=0.95, num_controls=100, assay_to_use = "RNA")
# 
# sum(auc_obj$auc_mat > computed_thresholds_df$threshold)/ncol(pre_data)
# 
# 
# pre_data <- AddMetaData(pre_data, auc_obj$auc_mat, col.name='resistance_score')
# pre_data <- AddMetaData(pre_data, ifelse(pre_data$resistance_score>computed_thresholds_df$threshold,"active","inactive"),col.name="resistance_active")


scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))

total_active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
active_cell_names <- total_active_cell_names[total_active_cell_names %in% colnames(pre_data)]

all_clusters <- as.numeric(levels(pre_data))

total_num_active <- length(total_active_cell_names)
total_num_inactive <- ncol(data)-total_num_active

ORs <- c()
pvals <- c()

for(i in all_clusters){
  curr_cluster_names <- colnames(pre_data)[pre_data$Cluster == i]
  
  curr_num_active <- sum(curr_cluster_names %in% active_cell_names)
  curr_num_inactive <- length(curr_cluster_names) - curr_num_active
  
  contin_table <- matrix(c(curr_num_active,total_num_active,curr_num_inactive,total_num_inactive), ncol=2)
  
  res <- fisher.test(contin_table)
  
  ORs <- append(ORs, res$estimate)
  pvals <- append(pvals, res$p.value)
  
}


df <- data.frame(cbind(paste0("", all_clusters),ORs))
colnames(df) <- c("cluster","or")
df$cluster <- factor(df$cluster, levels = all_clusters)
df$or <- as.numeric(df$or)

df$color <- ifelse(df$or > 1.5, "RAC","Non-RAC")

plot_title <- curr_cell_line

p <- ggplot(df)+
  geom_col(aes(x=cluster, y=or, fill=color))+
  scale_fill_manual(name="Cluster Type",values=c("lightblue","orange"))+
  geom_hline(yintercept=1.5, linetype="dashed", color = "red", size=1)+
  xlab("")+
  ylab("")+
  ggtitle(plot_title)+
  theme(legend.position="right",
        title = element_text(size=28, face="bold"),
        axis.text = element_text(size=30),
        legend.text = element_text(size=24),
        legend.title = element_text(size=26),
        legend.key.height = unit(1.5,"cm"),
        legend.key.width = unit(1.5,"cm"))



p

pre_or <- df$or



####################

drug_classes <- as.character(unique(data$pathway_level_1))
drug_classes <- drug_classes[-which(drug_classes == "Vehicle")]
drug_classes <- drug_classes[-which(drug_classes == "Other")]

curr_class <- drug_classes[2]

drug_data <- data[,data$pathway_level_1==curr_class]

drug_clusters <- drug_data@meta.data %>% 
  count(Cluster) %>% 
  filter(n > .01*ncol(drug_data)) %>% 
  pull(Cluster)

drug_data <- drug_data[,drug_data$Cluster %in% drug_clusters]



total_active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
active_cell_names <- total_active_cell_names[total_active_cell_names %in% colnames(drug_data)]

all_clusters <- as.numeric(levels(drug_data))

total_num_active <- length(total_active_cell_names)
total_num_inactive <- ncol(data)-total_num_active

ORs <- c()
pvals <- c()

for(i in all_clusters){
  curr_cluster_names <- colnames(drug_data)[drug_data$Cluster == i]
  
  curr_num_active <- sum(curr_cluster_names %in% active_cell_names)
  curr_num_inactive <- length(curr_cluster_names) - curr_num_active
  
  contin_table <- matrix(c(curr_num_active,total_num_active,curr_num_inactive,total_num_inactive), ncol=2)
  
  res <- fisher.test(contin_table)
  
  ORs <- append(ORs, res$estimate)
  pvals <- append(pvals, res$p.value)
  
}


df <- data.frame(cbind(paste0("", all_clusters),ORs))
colnames(df) <- c("cluster","or")
df$cluster <- factor(df$cluster, levels = all_clusters)
df$or <- as.numeric(df$or)

df$color <- ifelse(df$or > 1.5, "RAC","Non-RAC")

plot_title <- curr_cell_line

p <- ggplot(df)+
  geom_col(aes(x=cluster, y=or, fill=color))+
  scale_fill_manual(name="Cluster Type",values=c("lightblue","orange"))+
  geom_hline(yintercept=1.5, linetype="dashed", color = "red", size=1)+
  xlab("")+
  ylab("")+
  ggtitle(plot_title)+
  theme(legend.position="right",
        title = element_text(size=28, face="bold"),
        axis.text = element_text(size=30),
        legend.text = element_text(size=24),
        legend.title = element_text(size=26),
        legend.key.height = unit(1.5,"cm"),
        legend.key.width = unit(1.5,"cm"))



p






pre_cluster_cell_counts <- pre_data@meta.data %>% 
  count(Cluster) %>% 
  pull(n)


pre_cluster_fractions <- pre_cluster_cell_counts/ncol(pre_data)


drug_cluster_cell_counts <- drug_data@meta.data %>% 
  count(Cluster) %>% 
  pull(n)


drug_cluster_fractions <- drug_cluster_cell_counts/ncol(drug_data)











drug_or <- df$or


temp <- cbind(drug_clusters,log(drug_cluster_fractions/pre_cluster_fractions))
colnames(temp) <- c("cluster","or_fc")
temp <- as.data.frame(temp)

temp <- temp %>% 
  filter(cluster %in% c(4,9,12,13,14,16,18,19))

temp$cluster <- factor(as.character(temp$cluster),levels = drug_clusters)

ggplot(temp)+
  geom_col(aes(x=cluster,y=or_fc))


ggplot(rbind(pre_data@meta.data,drug_data@meta.data))+
  geom_boxplot(aes(x=Cluster,y=resistance_score, fill=treatment_stage))











