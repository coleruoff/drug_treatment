library(readxl)
library(tidyverse)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"

cell_lines <- c("A549","K562","MCF7")

crispr_ko_data <- read_xlsx("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/gottesman_crispr_data/CRISPR_KO_summary_negative.xlsx", sheet = 5)

supercluster_components <- readRDS(paste0(dataDirectory, "processed_data/supercluster_components.rds"))

supercluster_top_tfs <- readRDS(paste0(dataDirectory, "genesets/rac_supercluster_top_tf_list.rds"))

i <- 2
for(i in 1:2){
  
  # Get current supercluster top TFs
  curr_tfs <- supercluster_top_tfs[[i]]
  
  curr_tf_activity <- list()
  for(curr_cell_line in cell_lines){
    
    # Get TF activity for all current TFs in each component cell line for current supercluster
    tf_heatmap <- readRDS(paste0(dataDirectory, "processed_data/", curr_cell_line, "_tf_heatmap.rds"))
    
    curr_tf_activity <- append(curr_tf_activity, list(tf_heatmap[,supercluster_components[[i]][[curr_cell_line]]]))
    
    
  }
  
  common_tfs <- intersect(intersect(names(curr_tf_activity[[1]]),names(curr_tf_activity[[2]])),names(curr_tf_activity[[3]]))
  
  for(i in 1:3){
    curr_tf_activity[[i]] <- curr_tf_activity[[i]][names(curr_tf_activity[[i]]) %in% common_tfs]
  }
  
  temp <- data.frame(curr_tf_activity)
  
  # Calculate mean TF activity across cell lines
  mean_sc_tf_activity <- rowMeans(temp)
  
  names(sort(mean_sc_tf_activity, decreasing = T)[1:500])
  
  # Get current TFs logFC value from CRISPR experiment
  curr_tf_crispr <- crispr_ko_data %>% 
    filter(id %in% names(mean_sc_tf_activity)) %>% 
    dplyr::select(id, `neg|score`) %>% 
    data.frame()
  rownames(curr_tf_crispr) <- curr_tf_crispr$id
  
  # Only keep TF activity that are in CRISPR screen
  mean_sc_tf_activity <- mean_sc_tf_activity[names(mean_sc_tf_activity) %in% curr_tf_crispr$id]
  
  curr_tf_crispr <- curr_tf_crispr[names(mean_sc_tf_activity),]
  
  all(curr_tf_crispr$id == names(mean_sc_tf_activity))
  
  
  
  plot(log(mean_sc_tf_activity),abs(curr_tf_crispr$neg.score))
  abline(v=.5, col="blue")
  abline(h=.4, col="blue")
  
  cat(cor(abs(curr_tf_crispr$neg.score), log(mean_sc_tf_activity)))
  
  
  
}


temp <- data.frame(cbind(log(mean_sc_tf_activity),abs(curr_tf_crispr$neg.score)))


h <- .2
v <- 0

temp %>% 
  filter(X1 > v) %>% 
  nrow()


temp %>% 
  filter(X2 < h) %>% 
  nrow()


p1 <- 53/1117

p2 <- 452/1117


p1*p2

temp %>% 
  filter(X2 < .2 & X1 > .5) %>% 
  nrow()

(24/1117)/(p1*p2)






colnames(temp) <- cell_lines

pairs(cbind(temp,mean_sc_tf_activity))



summary(mean_sc_tf_activity)
temp <- data.frame(curr_tf_activity)

boxplot(1:100,101:1000)
wilcox.test(1:100,101:1000)



all(supercluster_top_tfs[[1]] %in% crispr_ko_data$id)

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, 1, 2), c("white", "red", "black"))
Heatmap(tf_heatmap, col = col_fun)

