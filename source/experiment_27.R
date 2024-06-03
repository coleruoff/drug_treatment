library(Seurat)

files <- list.files("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/watermelon_data/")

files <- files[1:8]
files <- files[-which(files %in% c("BT474_BREAST_processed.rds", "EFM192A_BREAST_processed.rds","MMACSF_SKIN_processed.rds"))]

all_de_list <- list()

for(curr_file in files){
  cat(curr_file, "\n")
  watermelon_data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/watermelon_data/",curr_file))
  
  watermelon_data$cycling <- ifelse(grepl("Non",watermelon_data$sample_pool),"noncycling",
                                    ifelse(grepl("Cyc",watermelon_data$sample_pool), "cycling","pre-treatment"))
  
  watermelon_data$time_point <- ifelse(grepl("_10_", watermelon_data$sample_pool), 1, 0)
  
  Idents(watermelon_data) <- watermelon_data$cycling
  
  de_res <- FindMarkers(watermelon_data, ident.1 = "cycling",ident.2 = "noncycling", test.use = "MAST")
  
  all_de_list <- append(all_de_list, list(de_res))
}

names(all_de_list) <- mapply(files, FUN = function(x) gsub("_processed.rds","",x))


# saveRDS(all_de_list, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/oren_cell_lines_cycling_vs_non_de.rds")

all_de_list <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/oren_cell_lines_cycling_vs_non_de.rds")


cycling_de_genes <- list()
for(i in all_de_list){
  temp <- i %>% 
    filter(p_val_adj < 0.05 & avg_log2FC > 0) %>% 
    arrange(desc(avg_log2FC)) %>% 
    rownames_to_column("gene") %>% 
    pull(gene)
  
  cycling_de_genes <- append(cycling_de_genes,list(temp))
    
}


consensus_cycling_signature <- list("cycling_signature"=find_consensus_geneset(cycling_de_genes, 4))

data <- all_data[["A549"]]

cells_AUC <- AUCell_run(data@assays$RNA$data, consensus_cycling_signature)


saveRDS(consensus_cycling_signature,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/consensus_cycling_signature.rds")

########################################################################################################################

source("source/read_in_all_cell_lines.R")
library(Seurat)
library(ggpubr)

###############################################################################
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
  
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line,"_processed_filtered_consensus_cycling_signature_aucell_scores.rds"))
  scores <- scale(scores)
  data <- AddMetaData(data, scores)
  
  
  df <- data@meta.data
  sc1_scores <- append(sc1_scores,df %>% 
                         filter(supercluster1 == 1) %>% 
                         pull(cycling_signature))
  
  sc2_scores <- append(sc2_scores,df %>% 
                         filter(supercluster2 == 1) %>% 
                         pull(cycling_signature))
  
  
  sc3_scores <- append(sc3_scores,df %>% 
                         filter(supercluster3 == 1) %>% 
                         pull(cycling_signature))
  
  
  rest_scores <- append(rest_scores,df %>% 
                          filter(rac == "nonrac") %>% 
                          pull(cycling_signature))
  
}

df <- data.frame(rbind(cbind("Supercluster 1",sc1_scores),cbind("Supercluster 2",sc2_scores),cbind("Supercluster 3",sc3_scores),cbind("Non-RAC",rest_scores)))

colnames(df) <- c("sc","scores")
df$scores <- as.numeric(df$scores)

# df <- df %>% 
#   filter(scores > 0)

my_comparisons <- list(c("Supercluster 1","Non-RAC"),c("Supercluster 2","Non-RAC"),c("Supercluster 3","Non-RAC"))

p <- ggboxplot(df, x="sc",y="scores",fill="sc")

plot_title <- paste0("Consensus Cycling Signature Score")

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


##################################################################
# PROLIFERATION INDEX

sc1_scores <- c()
sc2_scores <- c()
sc3_scores <- c()
rest_scores <- c()

for(curr_cell_line in cell_lines){
  data <- all_data[[curr_cell_line]]
  
  data <- AddMetaData(data, ifelse(data$Cluster == supercluster1_components[curr_cell_line], 1,0), "supercluster1")
  data <- AddMetaData(data, ifelse(data$Cluster == supercluster2_components[curr_cell_line], 1,0), "supercluster2")
  data <- AddMetaData(data, ifelse(data$Cluster == supercluster3_components[curr_cell_line], 1,0), "supercluster3")
  
  
  
  df <- data@meta.data
  sc1_scores <- append(sc1_scores,df %>% 
                         filter(supercluster1 == 1) %>% 
                         pull(proliferation_index))
  
  sc2_scores <- append(sc2_scores,df %>% 
                         filter(supercluster2 == 1) %>% 
                         pull(proliferation_index))
  
  
  sc3_scores <- append(sc3_scores,df %>% 
                         filter(supercluster3 == 1) %>% 
                         pull(proliferation_index))
  
  
  rest_scores <- append(rest_scores,df %>% 
                          filter(rac == "nonrac") %>% 
                          pull(proliferation_index))
  
}

df <- data.frame(rbind(cbind("Supercluster 1",sc1_scores),cbind("Supercluster 2",sc2_scores),cbind("Supercluster 3",sc3_scores),cbind("Non-RAC",rest_scores)))

colnames(df) <- c("sc","scores")
df$scores <- as.numeric(df$scores)

# df <- df %>% 
#   filter(scores > 0)

my_comparisons <- list(c("Supercluster 1","Non-RAC"),c("Supercluster 2","Non-RAC"),c("Supercluster 3","Non-RAC"))

p <- ggboxplot(df, x="sc",y="scores",fill="sc")

plot_title <- paste0("Proliferation Score")

p <- p + stat_compare_means(comparisons = my_comparisons,label = "p.format", method = "wilcox", label.x = 2.2, size=8)+
  ggtitle(plot_title)+
  xlab("")+
  ylab("")+
  theme(legend.position="right",
        title = element_text(size=20, face = "bold"),
        axis.text = element_text(size=20),
        legend.text = element_text(size=24),
        legend.title = element_text(size=26),
        legend.key.height = unit(1.5,"cm"),
        legend.key.width = unit(1.5,"cm"))+
  NoLegend()

p


