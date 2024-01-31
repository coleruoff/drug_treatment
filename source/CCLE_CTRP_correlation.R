library(tidyverse)
library(AUCell)
set.seed(42)
################################################################################

genesets_to_use <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_signatures.rds")

################################################################################
# Create dataframe of cell line scores for each gene set
################################################################################
CCLE_data <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/CCLE_data/CCLE_Object_normalized.rds")

dim(CCLE_data)

cells_AUC <- AUCell_run(CCLE_data@assays$RNA@data, genesets_to_use) 

cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, nCores=1, assign=TRUE)

colnames(cells_AUC@assays@data$AUC)

#Add trimmed cell line names to CCLE metadata
CCLE_data$trimmed_cell_line <- sapply(CCLE_data$Cell_line, FUN = function(x) strsplit(x,"_")[[1]][1])

#Score each cell line for all signatures
trimmed_CCLE_cell_lines_names <- unique(CCLE_data$trimmed_cell_line)

CCLE_df <- data.frame(matrix(NA,ncol=3, nrow=0))

#Score each cell line with AUCell and calculate average raw AUCell score for each cell line
for(curr_cell_line in trimmed_CCLE_cell_lines_names[1:100]){
  cat(curr_cell_line,"\n")
  
  curr_cell_names <- colnames(CCLE_data)[CCLE_data$trimmed_cell_line == curr_cell_line]
  
  mean_raw_scores <- rowMeans(cells_AUC@assays@data$AUC[,curr_cell_names])
  
  curr_df <- as.data.frame(cbind(curr_cell_line,mean_raw_scores)) %>% 
    rownames_to_column("geneset")
  
  
  CCLE_df <- rbind(CCLE_df,curr_df)
}

colnames(CCLE_df) <- c("geneset","cell_line","mean_score")

################################################################################
# Create viability dataframe
################################################################################

CTRP_data <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/CTRP_all_data.rds")

#Get viability data for each cell line
CTRP_cell_line_viability <- CTRP_data %>% 
  dplyr::select(c("ccl_name","area_under_curve","master_cpd_id"))

#Calculate mean viability across drugs for cell lines
CTRP_df <- CTRP_cell_line_viability %>% 
  group_by(ccl_name) %>% 
  mutate(avg_viability = mean(as.numeric(area_under_curve))) %>% 
  dplyr::select(ccl_name,avg_viability) %>% 
  distinct()

CTRP_df <- CTRP_cell_line_viability %>% 
  group_by(ccl_name) %>% 
  summarize(avg_viability = mean(as.numeric(area_under_curve))) %>% 
  dplyr::select(ccl_name,avg_viability) 


colnames(CTRP_df)[1] <- "cell_line"




final_df <- merge(CCLE_df, CTRP_df, by="cell_line")
final_df$mean_score <- as.numeric(final_df$mean_score)

ggplot(final_df, aes(x=avg_viability, y=mean_score))+
  geom_point()+
  geom_smooth(method = lm)+
  facet_wrap(~geneset, scales="free_y")
