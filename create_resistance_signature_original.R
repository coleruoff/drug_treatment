library(tidyverse)

# Raj 2023 breast

breast_de_res <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/raj_breast_vs_control_deseq_result.rds"))
# breast_de_res <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/raj_breast_resistant_vs_control_de_res.rds")

breast_resistant_genes <- data.frame(breast_de_res) %>% 
  filter(padj < 0.05 & log2FoldChange > 0) %>% 
  rownames_to_column() %>% 
  pull(rowname)

length(breast_resistant_genes)

# Raj 2023 melanoma

melanoma_de_res <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/raj_melanoma_vs_control_deseq_result.rds"))
# melanoma_de_res <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/raj_melanoma_resistant_vs_control_de_res.rds")

melanoma_resistant_genes <- data.frame(melanoma_de_res) %>% 
  filter(padj < 0.05 & log2FoldChange > 0) %>% 
  rownames_to_column() %>% 
  pull(rowname)

length(melanoma_resistant_genes)

# Raj wm9 2017
wm9_nodrug_vs_resistant_de <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/raj_2017_data/wm9_resNoDrugVsResistant.csv", row.names = 1)

wm9_resistant_markers <- wm9_nodrug_vs_resistant_de %>% 
  filter(padj < 0.05 & log2FoldChange > 0) %>% 
  arrange(desc(log2FoldChange)) %>% 
  pull(GeneSymbol)

wm9_nodrug_vs_drug_de <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/raj_2017_data/wm9_resNoDrugVsDrug.csv", row.names = 1)

wm9_early_drug_markers <- wm9_nodrug_vs_drug_de %>% 
  filter(padj < 0.05 & log2FoldChange > 2) %>% 
  arrange(desc(log2FoldChange)) %>% 
  pull(GeneSymbol)

wm9_resistant_markers_refined <- wm9_resistant_markers[!wm9_resistant_markers %in% wm9_early_drug_markers]

length(wm9_resistant_markers_refined)

# Raj wm983b 2017
wm983b_nodrug_vs_resistant_de <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/raj_2017_data/wm983b_resNoDrugVsResistant.csv", row.names = 1)

wm983b_resistant_markers <- wm983b_nodrug_vs_resistant_de %>% 
  filter(padj < 0.05 & log2FoldChange > 0) %>% 
  arrange(desc(log2FoldChange)) %>% 
  pull(GeneSymbol)

wm983b_nodrug_vs_drug_de <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/raj_2017_data/wm983b_resNoDrugVsDrug.csv", row.names = 1)

wm983b_early_drug_markers <- wm983b_nodrug_vs_drug_de %>% 
  filter(padj < 0.05 & log2FoldChange > 2) %>% 
  arrange(desc(log2FoldChange)) %>% 
  pull(GeneSymbol)

wm983b_resistant_markers_refined <- wm983b_resistant_markers[!wm983b_resistant_markers %in% wm983b_early_drug_markers]

length(wm983b_resistant_markers_refined)

common_resistance_genes <- intersect(intersect(wm9_resistant_markers_refined,wm983b_resistant_markers_refined),intersect(breast_resistant_genes,melanoma_resistant_genes))

length(common_resistance_genes)

common_resistance_signature <- list("common_resistance_signature" = common_resistance_genes)

saveRDS(common_resistance_signature, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/common_resistance_signature.rds")
saveRDS(common_resistance_signature, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/genesets/raj_watermelon_resistance_signature.rds")



sort(common_resistance_genes)

watermelon_de <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/watermelon_time_point_de.rds")


watermelon_day14_genes <- watermelon_de %>% 
  filter(p_val_adj < 0.05 & cluster == 14 & avg_log2FC > 0) %>% 
  arrange(desc(avg_log2FC)) %>% 
  pull(gene)

sort(all_resistant_genes$raj_watermelon_resistance_signature)
all_resistant_genes <- intersect(common_resistance_genes, watermelon_day14_genes)

all_resistant_genes <- list("raj_watermelon_resistance_signature" = all_resistant_genes)

saveRDS(all_resistant_genes, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_watermelon_resistance_signature.rds")
saveRDS(all_resistant_genes, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/genesets/raj_watermelon_resistance_signature.rds")

