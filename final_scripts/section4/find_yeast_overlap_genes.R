args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"

################################################################################

yeast_overlap_signatures <- readRDS(paste0(dataDirectory, "genesets/yeast_sc_overlap.rds"))

#Read in conversion table from OMA browser
yeast_conversion_table <- read.table(paste0(dataDirectory, "raw_data/yeast_data/human_candida_auris_source_id.txt"), fill = T, sep='\t')[,1:3]
colnames(yeast_conversion_table) <- c("HUMAN","YEAST", "MAPPING")

# Trim ensembl gene names
yeast_conversion_table$HUMAN <- gsub("\\..*", "", yeast_conversion_table$HUMAN)

genes_of_interest <- yeast_overlap_signatures[[2]]


human_ensembl_symbol_conversion <- AnnotationDbi::select(org.Hs.eg.db, 
                                                         keys = genes_of_interest,
                                                         columns = c("ENSEMBL", "SYMBOL"),
                                                         keytype = "SYMBOL")


genes_of_interest_ensembl <- human_ensembl_symbol_conversion$ENSEMBL


genes_of_interest_yeast <- yeast_conversion_table %>% 
  filter(HUMAN %in% genes_of_interest_ensembl) %>% 
  pull(YEAST)

genes_of_interest_yeast



