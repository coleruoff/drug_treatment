yeast_shared_genes <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/yeast_RAC_shared_genes.rds")

ecoli_shared_genes <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ecoli_RAC_shared_genes.rds")


write.csv(yeast_shared_genes, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/yeast_rac_shared_genes.csv")



capture.output(yeast_shared_genes, file = "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/yeast_rac_shared_genes.csv")

capture.output(ecoli_shared_genes, file = "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/ecoli_rac_shared_genes.csv")
