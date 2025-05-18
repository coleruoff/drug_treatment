args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"

################################################################################

AMR_de <- read_xlsx(paste0(dataDirectory, "raw_data/ecoli_data/AMR_de_genes.xlsx"))

colnames(AMR_de) <- AMR_de[1,]
AMR_de <- AMR_de[-1,]

experiments <- c("TETvs","MMCvs","IPMvs","CAZvs","KANvs","CIPvs","PMEvs","ERYvs","CHLvs")

ecoli_AMR_genesets_up <- list()
for(curr_experiment in experiments){
  curr_experiment_data <- AMR_de[,grepl(curr_experiment, colnames(AMR_de)) | grepl("Gene",colnames(AMR_de))]
  
  curr_experiment_data <- curr_experiment_data[,-c(2,3)]
  
  colnames(curr_experiment_data) <- c("gene_id","log2FC","pval","padj","significant","gene_name")
  
  curr_experiment_genes_up <- curr_experiment_data %>% 
    filter(as.numeric(padj) < 0.05 & as.numeric(log2FC) > 0) %>% 
    arrange(desc(as.numeric(log2FC))) %>% 
    pull(gene_id)
  
  ecoli_AMR_genesets_up <- append(ecoli_AMR_genesets_up, list(curr_experiment_genes_up))
  
}

names(ecoli_AMR_genesets_up) <- gsub("vs","",experiments)

################################################################################

ecoli_AMR_consensus_geneset_up <- find_consensus_geneset(ecoli_AMR_genesets_up,5)

################################################################################

#Read in conversion table from OMA browser
human_ecoli_conversion <- read.table(paste0(dataDirectory, "raw_data/ecoli_data/human_ecoli_source_id.txt"), fill = T, sep='\t')[,1:3]
colnames(human_ecoli_conversion) <- c("HUMAN", "ECOLI", "MAPPING")
human_ecoli_conversion$HUMAN <- gsub("\\..*", "", human_ecoli_conversion$HUMAN)

human_ecoli_conversion <- human_ecoli_conversion %>% 
  filter(ECOLI %in% ecoli_AMR_consensus_geneset_up)

# Gene human orthologs of ecoli genes
ecoli_human_orthologs_ensembl <- human_ecoli_conversion$HUMAN

# Add HUMAN ensembl to symbol conversion
human_ensembl_symbol_conversion <- AnnotationDbi::select(org.Hs.eg.db, 
                                                         keys = ecoli_human_orthologs_ensembl,
                                                         columns = c("ENSEMBL", "SYMBOL"),
                                                         keytype = "ENSEMBL")

ecoli_human_orthologs_symbol <- human_ensembl_symbol_conversion$SYMBOL

ecoli_human_orthologs_up <- list("ecoli_human_orthologs_up"=ecoli_human_orthologs_symbol)

saveRDS(ecoli_human_orthologs_up, paste0(dataDirectory, "genesets/ecoli_human_orthologs_up.rds"))















