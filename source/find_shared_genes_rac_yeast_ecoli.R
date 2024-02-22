
cell_lines <- c("A549","K562","MCF7")
global_rac_type_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_rac_type_signatures.rds")
global_rac_type1_signatures <- global_rac_type_signatures[grepl("type1",names(global_rac_type_signatures))]
names(global_rac_type1_signatures) <- cell_lines


upregulated_yeast_human_orthologs_symbol <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/orthologs/upregulated_yeast_human_orthologs_symbol.rds")
length(upregulated_yeast_human_orthologs_symbol)

upregulated_ecoli_human_orthologs_symbol <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/orthologs/upregulated_ecoli_human_orthologs_symbol.rds")
length(upregulated_ecoli_human_orthologs_symbol)


curr_cell_line <- "MCF7"

yeast_shared_genes <- list()
ecoli_shared_genes <- list()
for(curr_cell_line in cell_lines){
  curr_yeast_shared_genes <- upregulated_yeast_human_orthologs_symbol[upregulated_yeast_human_orthologs_symbol %in% global_rac_type1_signatures[[curr_cell_line]]]
  
  yeast_shared_genes <- append(yeast_shared_genes,list(curr_yeast_shared_genes))
  
  
  curr_ecoli_shared_genes <- upregulated_ecoli_human_orthologs_symbol[upregulated_ecoli_human_orthologs_symbol %in% global_rac_type1_signatures[[curr_cell_line]]]
  
  ecoli_shared_genes <- append(ecoli_shared_genes,list(curr_ecoli_shared_genes))
}



sort(global_rac_type1_signatures[["MCF7"]])
