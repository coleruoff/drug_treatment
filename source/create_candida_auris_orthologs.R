library(readxl)
source("source/human_to_organism_conversion_function.R")

yeast_de <- read_xlsx("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/yeast_data/candida_auris_de.XLSX", col_names = T)
colnames(yeast_de) <- yeast_de[1,]
yeast_de <- yeast_de[-1,]

genes <- yeast_de$gene_id

upregulated_genes <- yeast_de %>% 
  filter(as.numeric(padj) < 0.05 & up_down == "Up") %>% 
  pull(gene_id)

length(upregulated_genes)

candida_auris_upregulated_orthologs <- convert_yeast_to_human_genes(upregulated_genes)

list(candida_auris_upregulated_orthologs)

saveRDS(candida_auris_upregulated_orthologs, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/yeast_upregulated_orthologs.rds")
