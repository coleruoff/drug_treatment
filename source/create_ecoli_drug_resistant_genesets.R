library(readxl)
library(tidyverse)

AMR_de <- read_xlsx("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/ecoli_data/AMR_de_genes.xlsx")


colnames(AMR_de) <- AMR_de[1,]
AMR_de <- AMR_de[-1,]


experiments <- c("TETvs","MMCvs","IPMvs","CAZvs","KANvs","CIPvs","PMEvs","ERYvs","CHLvs")


ecoli_AMR_genesets_up <- list()
ecoli_AMR_genesets_down <- list()
for(curr_experiment in experiments){
  curr_experiment_data <- AMR_de[,grepl(curr_experiment, colnames(AMR_de)) | grepl("Gene",colnames(AMR_de))]
  
  curr_experiment_data <- curr_experiment_data[,-c(2,3)]
  
  colnames(curr_experiment_data) <- c("gene_id","log2FC","pval","padj","significant","gene_name")
  
  curr_experiment_genes_up <- curr_experiment_data %>% 
    filter(as.numeric(padj) < 0.05 & as.numeric(log2FC) > 0) %>% 
    pull(gene_id)
  
  ecoli_AMR_genesets_up <- append(ecoli_AMR_genesets_up, list(curr_experiment_genes_up))
  
  curr_experiment_genes_down <- curr_experiment_data %>% 
    filter(as.numeric(padj) < 0.05 & as.numeric(log2FC) < 0) %>% 
    pull(gene_id)
  
  ecoli_AMR_genesets_down <- append(ecoli_AMR_genesets_down, list(curr_experiment_genes_down))
  
  
}



names(ecoli_AMR_genesets_up) <- gsub("vs","",experiments)

names(ecoli_AMR_genesets_down) <- gsub("vs","",experiments)



saveRDS(ecoli_AMR_genesets_up, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/ecoli_data/ecoli_AMR_genesets_up.rds")
saveRDS(ecoli_AMR_genesets_down, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/ecoli_data/ecoli_AMR_genesets_down.rds")


