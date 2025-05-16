args = commandArgs(trailingOnly=TRUE)
# dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])

source("revision_scripts/drug_treatment_functions.R")
set.seed(42)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"

################################################################################

yeast_de <- read_xlsx(paste0(dataDirectory, "raw_data/yeast_data/candida_auris_de.XLSX"), col_names = T)
colnames(yeast_de) <- yeast_de[1,]
yeast_de <- yeast_de[-1,]

upregulated_genes <- yeast_de %>%
  filter(as.numeric(padj) < 0.05 & log2FoldChange > 0) %>%
  pull(gene_id)

#Read in conversion table from OMA browser
human_yeast_conversion <- read.table(paste0(dataDirectory, "raw_data/yeast_data/human_candida_auris_source_id.txt"), fill = T, sep='\t')[,1:3]
colnames(human_yeast_conversion) <- c("HUMAN","YEAST", "MAPPING")

# Trim ensembl gene names
human_yeast_conversion$HUMAN <- gsub("\\..*", "", human_yeast_conversion$HUMAN)

#Select genes that map to one yeast gene and those that are present in yeast upregulated list
upregulated_yeast_human_orthologs_table <- human_yeast_conversion %>% 
  filter(YEAST %in% upregulated_genes)

# Gene human orthologs of yeast genes
yeast_human_orthologs_ensembl <- upregulated_yeast_human_orthologs_table$HUMAN

# Add HUMAN ensembl to symbol conversion
human_ensembl_symbol_conversion <- AnnotationDbi::select(org.Hs.eg.db, 
                                                         keys = yeast_human_orthologs_ensembl,
                                                         columns = c("ENSEMBL", "SYMBOL"),
                                                         keytype = "ENSEMBL")

yeast_human_orthologs_symbol <- human_ensembl_symbol_conversion$SYMBOL

yeast_human_orthologs_up <- list("yeast_human_orthologs_up"=yeast_human_orthologs_symbol)

saveRDS(yeast_human_orthologs_up, paste0(dataDirectory, "genesets/yeast_human_orthologs_up.rds"))
