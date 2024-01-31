

cell_lines <- c("A549","K562","MCF7")

all_active_vs_inactive_signatures <- list()
for(curr_cell_line in cell_lines){
  de_result <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line,"_all_resistant_cells_de.rds"))
  
  
  active_vs_inactive_signature <- de_result %>% 
    filter(cluster=="resistant" & p_val_adj < 0.05, avg_log2FC > 0) %>% 
    arrange(desc(avg_log2FC)) %>% 
    pull(gene)
  
  
  all_active_vs_inactive_signatures <- append(all_active_vs_inactive_signatures, list(active_vs_inactive_signature))
  
}


names(all_active_vs_inactive_signatures) <- paste0(cell_lines, "_active_vs_inactive_signature")




saveRDS(all_active_vs_inactive_signatures, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/active_vs_inactive_signatures.rds"))
