args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])


crispr_ko_data <- read_xlsx("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/gottesman_crispr_data/CRISPR_KO_summary_negative.xlsx", sheet = 5)

crispr_ko_ranks <- crispr_ko_data %>% 
  arrange(`neg|fdr`) %>% 
  mutate(ranks = 1:nrow(.)) %>% 
  mutate(adj_p_value = `neg|fdr`) %>% 
  dplyr::select(id,ranks,adj_p_value)


final_df <- readRDS(paste0(dataDirectory, "supercluster_tf_top_bottom_crispr_ranks_data.rds"))

top_genes <- list()
list_names <- c()

for(sc_to_use in paste0("Supercluster ", 1:2)){
  curr_ranks <- final_df %>% 
    filter(geneset == sc_to_use & top == "Top") %>% 
    pull(rank)
  
  curr_genes <- crispr_ko_ranks %>% 
    filter(ranks %in% curr_ranks & adj_p_value < 0.05) %>%
    arrange(ranks) %>%
    pull(id) 
  
  top_genes <- append(top_genes, list(curr_genes))
}


names(top_genes) <- paste0("Supercluster ", 1:2)

max_len <- max(lengths(top_genes))

top_genes <- lapply(top_genes, FUN = function(x) append(x,rep("",max_len-length(x))))

top_genes <- data.frame(top_genes)

write.xlsx(top_genes, file=paste0(dataDirectory, "supplementary_data/supplementary_table3.xlsx"), sheetName="supercluster_top_tfs", row.names=FALSE)
