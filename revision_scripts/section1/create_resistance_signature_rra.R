source("revision_scripts/drug_treatment_functions.R")

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"

################################################################################

all_cell_line_de_res <- readRDS(paste0(dataDirectory, "de_results/oren_all_cell_line_de_res.rds"))

# Remove PC9 because we use it for validation
all_cell_line_de_res <- all_cell_line_de_res[-1]

#No DE genes(?)
all_cell_line_de_res <- all_cell_line_de_res[-7]

# Remove breast cancer cell line due to few cells and DE genes
# all_cell_line_de_res <- all_cell_line_de_res[-6]

# Remove intestine because its not significantly overlapped with others
# all_cell_line_de_res <- all_cell_line_de_res[-3]

all_rank_lists <- list()

for(i in 1:length(all_cell_line_de_res)){
  
  
  # if(i == 6){
  #   curr_mat <- all_cell_line_de_res[[i]] %>%
  #     filter(cluster == 10 & avg_log2FC > 0 & p_val < 0.05) %>%
  #     arrange(desc(avg_log2FC)) %>%
  #     select(gene)
  # } else {
  #   curr_mat <- all_cell_line_de_res[[i]] %>%
  #     filter(cluster == 10 & avg_log2FC > 0 & p_val_adj < 0.05) %>%
  #     arrange(desc(avg_log2FC)) %>%
  #     select(gene)
  # }
  
  curr_mat <- all_cell_line_de_res[[i]] %>%
    filter(cluster == 10 & avg_log2FC > 0 & p_val_adj < 0.05) %>%
    arrange(desc(avg_log2FC)) %>%
    select(gene)
  
  
  # curr_mat <- all_cell_line_de_res[[i]] %>%
  #   filter(cluster == 10 & avg_log2FC > 0) %>%
  #   arrange((p_val_adj)) %>%
  #   select(gene)
  
  
  all_rank_lists <- append(all_rank_lists,curr_mat)
  
}

res <- aggregateRanks(all_rank_lists, method = "RRA")

consensus_genes <- subset(res, Score < 0.05)

consensus_genes <- consensus_genes$Name

resistance_signature <- list("resistance_signature" = consensus_genes)

saveRDS(resistance_signature, paste0(dataDirectory, "genesets/resistance_signature.rds"))







