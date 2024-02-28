cell_lines <- c("A549","K562","MCF7")

curr_cell_line <- cell_lines[1]
cat(curr_cell_line,"\n")

RACs <- list(c(4,9,12,13,14,16,18),c(4,5,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")

emergent <- list(c(14:19),c(9,11),c(13,15,18))
names(emergent) <- cell_lines

final_df <- data.frame(matrix(NA,ncol=2,nrow=0))
for(curr_cell_line in cell_lines){
  
  #Read in cell line data
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
  threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
  active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
  clusters_of_interest <- RACs[[curr_cell_line]]
  
  #Add metadata for RAC and Cell Group
  data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, "1", ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, paste0(data$Cluster, "_1"), ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), paste0(data$Cluster, "_2"), paste0(data$Cluster, "_0"))), col.name = "cell_cluster_group")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & data$Cluster %in% emergent[[curr_cell_line]], "emergent_rac", ifelse(data$rac == "rac" & (!data$Cluster %in% emergent[[curr_cell_line]]), "non_emergent_rac", "non_rac")), col.name = "emergent_rac")
  
  
  df <- as.data.frame(cbind(data@assays$RNA$data["ABCB1",],data$emergent_rac,curr_cell_line))
  
  
  colnames(df) <- c("ABCB1","cell_group","cell_line")
  
  df$ABCB1 <- as.numeric(df$ABCB1)
  
  final_df <- rbind(final_df, df)
}



final_df <- final_df %>%
  filter(ABCB1 > 0)

ggplot(final_df)+
  geom_boxplot(aes(x=cell_group, y=ABCB1))+
  facet_wrap(~cell_line)
