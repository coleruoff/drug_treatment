library(readxl)
library(tidyverse)
library(AUCell)
set.seed(42)
################################################################################
ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
  row_names = rownames(X)
  num_genes = nrow(X)
  gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})
  
  # Ranks for genes
  R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')
  
  # Calculate enrichment score (es) for each sample (column)
  es = apply(R, 2, function(R_col) {
    gene_ranks = order(R_col, decreasing = TRUE)
    
    # Calc es for each gene set
    es_sample = sapply(gene_sets, function(gene_set_idx) {
      # pos: match (within the gene set)
      # neg: non-match (outside the gene set)
      indicator_pos = gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos
      
      rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha
      
      step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)
      
      step_cdf_diff = step_cdf_pos - step_cdf_neg
      
      # Normalize by gene number
      if (scale) step_cdf_diff = step_cdf_diff / num_genes
      
      # Use ssGSEA or not
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    unlist(es_sample)
  })
  
  if (length(gene_sets) == 1) es = matrix(es, nrow = 1)
  
  # Normalize by absolute diff between max and min
  if (norm) es = es / diff(range(es))
  
  # Prepare output
  rownames(es) = names(gene_sets)
  colnames(es) = colnames(X)
  return(es)
}

genesets_to_use <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_signatures.rds")

CCLE_data <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/CCLE_data/CCLE_Object_normalized.rds")
GDSC_data <- read_xlsx("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/GDSC/GDSC2_fitted_dose_response_27Oct23.xlsx")

################################################################################

#Add trimmed cell line names to CCLE metadata
CCLE_data$trimmed_cell_line <- sapply(CCLE_data$Cell_line, FUN = function(x) strsplit(x,"_")[[1]][1])

shared_cell_lines <- unique(CCLE_data$trimmed_cell_line)[unique(CCLE_data$trimmed_cell_line) %in% unique(GDSC_data$CELL_LINE_NAME)]

CCLE_data_filtered <- CCLE_data[,CCLE_data$trimmed_cell_line %in% shared_cell_lines]

dim(CCLE_data_filtered)

CCLE_scores <- ssgsea(as.matrix(CCLE_data_filtered@assays$RNA@data), genesets_to_use)

colnames(CCLE_scores)


################################################################################
# Create dataframe of cell line scores for each gene set
################################################################################

#Score each cell line for all signatures
trimmed_CCLE_cell_lines_names <- unique(CCLE_data_filtered$trimmed_cell_line)

CCLE_df <- data.frame(matrix(NA,ncol=3, nrow=0))

#Score each cell line with AUCell and calculate average raw AUCell score for each cell line

curr_cell_line <- trimmed_CCLE_cell_lines_names[1]
for(curr_cell_line in trimmed_CCLE_cell_lines_names){
  cat(curr_cell_line,"\n")
  
  curr_sample_names <- colnames(CCLE_data)[CCLE_data$trimmed_cell_line == curr_cell_line]
  
  mean_scores <- rowMeans(CCLE_scores[,curr_sample_names])
  
  curr_df <- as.data.frame(cbind(curr_cell_line,mean_scores)) %>% 
    rownames_to_column("geneset")
  
  
  CCLE_df <- rbind(CCLE_df,curr_df)
}

colnames(CCLE_df) <- c("geneset","cell_line","mean_score")

################################################################################
# Create viability dataframe
################################################################################

GDSC_data <- GDSC_data %>% 
  filter(CELL_LINE_NAME %in% shared_cell_lines)


colnames(GDSC_data)

GDSC_df <- GDSC_data %>% 
  dplyr::select(c("CELL_LINE_NAME","AUC")) %>% 
    group_by(CELL_LINE_NAME) %>% 
    mutate(avg_auc = mean(as.numeric(AUC))) %>% 
    dplyr::select(CELL_LINE_NAME,avg_auc) %>% 
    distinct()


colnames(GDSC_df)[1] <- "cell_line"



##############################################

final_df <- merge(CCLE_df, GDSC_df, by="cell_line")
final_df$mean_score <- as.numeric(final_df$mean_score)





ggplot(final_df, aes(x=avg_auc, y=mean_score))+
  geom_point()+
  geom_smooth(method = lm)+
  facet_wrap(~geneset, scales="free_y")




for(i in unique(final_df$geneset)){
  
  temp <- final_df %>% 
    filter(geneset == i) %>% 
    select(avg_auc,mean_score)
  
  cat(cor(temp[,1],temp[,2]), "\n")
  
}

