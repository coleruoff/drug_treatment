args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])

library(tidyverse)

files_to_delete <- c(paste0(dataDirectory,"genesets/common_resistance_signature.rds"),
                     paste0(dataDirectory,"genesets/raj_watermelon_resistance_signature.rds"),
                     paste0(dataDirectory, "genesets/all_resistant_genesets.rds"),
                     paste0(dataDirectory, "genesets/drug_resistance_signatures.rds"))

for(to_be_deleted in files_to_delete){
  if (file.exists(to_be_deleted)) {
    #Delete file if it exists
    file.remove(to_be_deleted)
    cat(to_be_deleted)
  }
}

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

################################################################################

# Goyal 2023 breast

breast_de_res <- readRDS(paste0(dataDirectory,"de_results/raj_breast_vs_control_deseq_result.rds"))

breast_resistant_genes <- data.frame(breast_de_res) %>% 
  filter(padj < 0.05 & log2FoldChange > 0) %>% 
  rownames_to_column() %>% 
  pull(rowname)

# Goyal 2023 melanoma

melanoma_de_res <- readRDS(paste0(dataDirectory,"de_results/raj_melanoma_vs_control_deseq_result.rds"))

melanoma_resistant_genes <- data.frame(melanoma_de_res) %>% 
  filter(padj < 0.05 & log2FoldChange > 0) %>% 
  rownames_to_column() %>% 
  pull(rowname)

# Raj wm9 2017
wm9_nodrug_vs_resistant_de <- read.csv(paste0(dataDirectory, "raw_data/raj_2017_data/wm9_resNoDrugVsResistant.csv"), row.names = 1)

wm9_resistant_genes <- wm9_nodrug_vs_resistant_de %>% 
  filter(padj < 0.05 & log2FoldChange > 0) %>% 
  arrange(desc(log2FoldChange)) %>% 
  pull(GeneSymbol)

wm9_nodrug_vs_drug_de <- read.csv(paste0(dataDirectory, "raw_data/raj_2017_data/wm9_resNoDrugVsDrug.csv"), row.names = 1)

wm9_early_drug_genes <- wm9_nodrug_vs_drug_de %>% 
  filter(padj < 0.05 & log2FoldChange > 2) %>% 
  arrange(desc(log2FoldChange)) %>% 
  pull(GeneSymbol)

wm9_resistant_genes_refined <- wm9_resistant_genes[!wm9_resistant_genes %in% wm9_early_drug_genes]

# Raj wm983b 2017
wm983b_nodrug_vs_resistant_de <- read.csv(paste0(dataDirectory,"raw_data/raj_2017_data/wm983b_resNoDrugVsResistant.csv"), row.names = 1)

wm983b_resistant_genes <- wm983b_nodrug_vs_resistant_de %>% 
  filter(padj < 0.05 & log2FoldChange > 0) %>% 
  arrange(desc(log2FoldChange)) %>% 
  pull(GeneSymbol)

wm983b_nodrug_vs_drug_de <- read.csv(paste0(dataDirectory,"raw_data/raj_2017_data/wm983b_resNoDrugVsDrug.csv"), row.names = 1)

wm983b_early_drug_genes <- wm983b_nodrug_vs_drug_de %>% 
  filter(padj < 0.05 & log2FoldChange > 2) %>% 
  arrange(desc(log2FoldChange)) %>% 
  pull(GeneSymbol)

wm983b_resistant_genes_refined <- wm983b_resistant_genes[!wm983b_resistant_genes %in% wm983b_early_drug_genes]

# Get intersection of all upregulated genes
initial_resistance_genes <- intersect(intersect(wm9_resistant_genes_refined,wm983b_resistant_genes_refined),intersect(breast_resistant_genes,melanoma_resistant_genes))

initial_resistance_signature <- list("initial_resistance_signature" = initial_resistance_genes)

saveRDS(initial_resistance_signature, paste0(dataDirectory,"genesets/initial_resistance_signature.rds"))

################################################################################
# Add Oren 2021 day 14 genes to resistance signature

watermelon_de <- readRDS(paste0(dataDirectory,"de_results/watermelon_time_point_de.rds"))

watermelon_day14_genes <- watermelon_de %>% 
  filter(p_val_adj < 0.05 & cluster == 14 & avg_log2FC > 0) %>% 
  arrange(desc(avg_log2FC)) %>% 
  pull(gene)

all_resistant_genes <- intersect(initial_resistance_genes, watermelon_day14_genes)

all_resistant_genes <- list("raj_watermelon_resistance_signature" = all_resistant_genes)

saveRDS(all_resistant_genes, paste0(dataDirectory,"genesets/raj_watermelon_resistance_signature.rds"))

################################################################################

all_resistant_genesets <- list(wm9_resistant_genes_refined,wm983b_resistant_genes_refined,breast_resistant_genes,melanoma_resistant_genes,watermelon_day14_genes)
names(all_resistant_genesets) <- c("shaffer_2017_w9","shaffer_2017_w983b","goyal_2023_breast","goyal_2023_melanoma","oren_2021_pc9")

saveRDS(all_resistant_genesets, paste0(dataDirectory, "genesets/all_resistant_genesets.rds"))

####
drug_resistance_signatures <- list(common_resistance_genes, all_resistant_genes$raj_watermelon_resistance_signature)
names(drug_resistance_signatures) <- c("shaffer_goyal_signature","shaffer_goyal_oren_signature")

saveRDS(drug_resistance_signatures, paste0(dataDirectory, "genesets/drug_resistance_signatures.rds"))






