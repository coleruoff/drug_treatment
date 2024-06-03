library(Seurat)
library(AUCell)

# raj_resistant_breast <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/diverse_clonal_fates_data/raj_resistant_breast_processed.rds")
# raj_resistant_melanoma <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/drug_treatment/processed_data/diverse_clonal_fates_data/raj_resistant_melanoma_processed.rds"))
# 
# raj_control_breast <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/drug_treatment/processed_data/memorySeq_data/MDA_control_tpm.rds"))
# raj_control_melanoma <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/drug_treatment/processed_data/memorySeq_data/WM989_control_tpm.rds"))

MPs <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ITH_meta_programs.rds")
specifc_mps <- c("MP39 Metal-response","MP31 Alveolar","MP29 NPC/OPC","MP28 Oligo normal","MP27 Oligo Progenitor","MP38 Glutathione","MP41 Unassigned","MP35 Hemato-related-I","MP37 Hemato-related-II","MP32 Skin-pigmentation","MP36 IG","MP16 MES (glioma)","MP15 EMT IV")
MPs <- MPs[!names(MPs) %in% specifc_mps]

files <- list.files("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/watermelon_data/")

files <- files[1:8]
total_df <- list()

for(curr_file in files){
  cat(curr_file, "\n")
  watermelon_data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/watermelon_data/",curr_file))
  
  watermelon_data$time_point <- ifelse(grepl("_10_", watermelon_data$sample_pool), 1, 0)
  
  cells_AUC <- AUCell_run(watermelon_data@assays$RNA@data, MPs)
  
  curr_df <- cbind(t(as.matrix(cells_AUC@assays@data$AUC)),watermelon_data$time_point)
  
  total_df <- append(total_df, list(curr_df))
}

total_df <- do.call("rbind",total_df)
dim(total_df)

cancer_data <- total_df
colnames(cancer_data)[ncol(cancer_data)] <- "labels"

set.seed(42)  # for reproducibility
train_idx <- sample(nrow(cancer_data), 0.7 * nrow(cancer_data))
train_data <- cancer_data[train_idx, ]
test_data <- cancer_data[-train_idx, ]

# Create the logistic regression model
model <- glm(labels ~ ., data = data.frame(train_data), family = "binomial")

saveRDS(model, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/experiment_data/experiment_29_model.rds")

model <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/experiment_data/experiment_29_model.rds")


summary(model)

cancer_data["labels",]
# Make predictions on the testing set
predictions <- predict(model, newdata = data.frame(test_data), type = "response")

# Convert probabilities to binary predictions
threshold <- 0.5
binary_predictions <- ifelse(predictions > threshold, 1, 0)

# Evaluate the model
accuracy <- mean(binary_predictions == test_data[,"labels"])



temp <- summary(model)
sort(temp$coefficients[,"Estimate"],decreasing = T)



################################################################
setwd("/data/ruoffcj/projects/drug_treatment/")
library(tidyverse)
library(Seurat)
library(ggpubr)

curr_cell_line <- "K562"
use_pre <- F

cell_lines <- c("A549","K562","MCF7")

plots <- list()
for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line,"\n")
  
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  if(use_pre){
    pre_clusters <- data@meta.data %>%
      count(treatment_stage,Cluster) %>%
      filter(treatment_stage == "pre" & n > 10) %>%
      pull(Cluster)
    
    data <- data[,data$treatment_stage=='pre' & data$Cluster %in% pre_clusters]
  } else{
    data <- data[,data$treatment_stage=='post']
  }
  
  
  trapnell_scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line, "_processed_filtered_ITH_meta_programs_aucell_scores.rds"))
  
  trapnell_scores <- trapnell_scores[,!colnames(trapnell_scores) %in% specifc_mps]
  
  predictions <- predict(model, newdata = data.frame(trapnell_scores), type = "response")
  
  binary_predictions <- ifelse(predictions > threshold, 1, 0)
  
  
  active_cell_names <- names(binary_predictions)[binary_predictions==1]
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
  
  
  df <- data.frame(cbind(paste0("", all_clusters),ORs))
  colnames(df) <- c("cluster","or")
  df$cluster <- factor(df$cluster, levels = all_clusters)
  df$or <- as.numeric(df$or)
  
  df$color <- ifelse(df$or > 1.5, "RAC","Non-RAC")
  
  # plot_title <- paste0("Resistant Active/Inactive Odds Ratio (", curr_cell_line, " Pre-Treatment)")
  plot_title <- curr_cell_line
  # 
  # if(use_pre){
  #   png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/active_inactive_or_figures/",curr_cell_line,"_barplots_pre.png"),
  #       width=1000, height = 500)
  # } else {
  #   png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/active_inactive_or_figures/",curr_cell_line,"_barplots.png"),
  #       width=1000, height = 500)
  # }
  
  
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
  
  plots <- append(plots,list(p))
}


# png(paste0("/data/ruoffcj/projects/drug_treatment/final_figures/figure_1c.png"),
#      width=20, height=20, units= "in", res=300)

figure <- ggarrange(plotlist = plots, nrow=3, common.legend = T,legend=c("right"))

p <- annotate_figure(figure, left = text_grob("Odds Ratio", rot = 90, vjust = 1, size=35, face="bold"),
                     bottom = text_grob("Clusters", size=35, face="bold"),
                     top=text_grob("Resistant Active/Inactive Odds Ratio", size=40, face="bold"))

print(p)










