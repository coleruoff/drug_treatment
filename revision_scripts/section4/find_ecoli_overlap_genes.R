args = commandArgs(trailingOnly=TRUE)
# dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])

source("revision_scripts/drug_treatment_functions.R")
set.seed(42)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"

################################################################################
# Read in conversion table
ecoli_conversion_table <- read.table(paste0(dataDirectory, "raw_data/ecoli_data/human_ecoli_source_id.txt"), fill = T, sep='\t')[,1:3]
colnames(ecoli_conversion_table) <- c("HUMAN", "ECOLI", "MAPPING")
ecoli_conversion_table$HUMAN <- gsub("\\..*", "", ecoli_conversion_table$HUMAN)


# Read in overlapping genes
ecoli_overlap_signatures <- readRDS(paste0(dataDirectory, "genesets/ecoli_sc_overlap.rds"))

ecoli_overlap_signatures <- ecoli_overlap_signatures[[1]]

# Get SYMBOL to ENTREZID conversion for human genes
human_entrezid_symbol_conversion <- AnnotationDbi::select(org.Hs.eg.db, 
                                                         keys = ecoli_overlap_signatures[[1]],
                                                         columns = c("SYMBOL", "ENSEMBL"),
                                                         keytype = "SYMBOL")
# Select ENTREZID human genes
ecoli_overlap_ensembl<- human_entrezid_symbol_conversion$ENSEMBL


# Get ENTREZID ecoli overlap genes
ecoli_overlap_ecoli_genes_entrez <- ecoli_conversion_table %>% 
  filter(HUMAN %in% ecoli_overlap_ensembl) %>% 
  pull(ECOLI)


ecoli_overlap_ecoli_genes_entrez
