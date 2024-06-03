library(tidyverse)

# Raj 2023 breast
raj_resistant_cancer_type <- "breast"

breast_de_res <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/drug_treatment/raj_", raj_resistant_cancer_type,"_vs_control_deseq_result.rds"))

breast_resistant_genes <- data.frame(breast_de_res) %>% 
  filter(padj < 0.05 & log2FoldChange > 0) %>% 
  rownames_to_column() %>% 
  pull(rowname)

# Raj 2023 melanoma
raj_resistant_cancer_type <- "melanoma"

melanoma_de_res <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/drug_treatment/raj_", raj_resistant_cancer_type,"_vs_control_deseq_result.rds"))

melanoma_resistant_genes <- data.frame(melanoma_de_res) %>% 
  filter(padj < 0.05 & log2FoldChange > 0) %>% 
  rownames_to_column() %>% 
  pull(rowname)

# Raj wm9 2017
wm9_nodrug_vs_resistant_de <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/raj_2017_data/wm9_resNoDrugVsResistant.csv", row.names = 1)

wm9_resistant_genes <- wm9_nodrug_vs_resistant_de %>% 
  filter(padj < 0.05 & log2FoldChange > 0) %>% 
  arrange(desc(log2FoldChange)) %>% 
  pull(GeneSymbol)

wm9_nodrug_vs_drug_de <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/raj_2017_data/wm9_resNoDrugVsDrug.csv", row.names = 1)

wm9_early_drug_genes <- wm9_nodrug_vs_drug_de %>% 
  filter(padj < 0.05 & log2FoldChange > 2) %>% 
  arrange(desc(log2FoldChange)) %>% 
  pull(GeneSymbol)

wm9_resistant_genes <- wm9_resistant_genes[!wm9_resistant_genes %in% wm9_early_drug_genes]

# Raj wm983b 2017
wm983b_nodrug_vs_resistant_de <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/raj_2017_data/wm983b_resNoDrugVsResistant.csv", row.names = 1)

wm983b_resistant_genes <- wm983b_nodrug_vs_resistant_de %>% 
  filter(padj < 0.05 & log2FoldChange > 0) %>% 
  arrange(desc(log2FoldChange)) %>% 
  pull(GeneSymbol)

wm983b_nodrug_vs_drug_de <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/raj_2017_data/wm983b_resNoDrugVsDrug.csv", row.names = 1)

wm983b_early_drug_genes <- wm983b_nodrug_vs_drug_de %>% 
  filter(padj < 0.05 & log2FoldChange > 2) %>% 
  arrange(desc(log2FoldChange)) %>% 
  pull(GeneSymbol)

wm983b_resistant_gene <- wm983b_resistant_genes[!wm983b_resistant_genes %in% wm983b_early_drug_genes]

# Watermelon day 14
watermelon_de <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/watermelon_time_point_de.rds")

watermelon_resistant_genes <- watermelon_de %>% 
  filter(p_val_adj < 0.05 & cluster == 14 & avg_log2FC > 0) %>% 
  arrange(desc(avg_log2FC)) %>% 
  pull(gene)


#################################################################################
raj_breast_ranks <- data.frame(cbind(breast_resistant_genes,1:length(breast_resistant_genes)))
colnames(raj_breast_ranks) <- c("gene","rank")

raj_melanoma_ranks <- data.frame(cbind(melanoma_resistant_genes,1:length(melanoma_resistant_genes)))
colnames(raj_melanoma_ranks) <- c("gene","rank")

wm9_ranks <- data.frame(cbind(wm9_resistant_genes,1:length(wm9_resistant_genes)))
colnames(wm9_ranks) <- c("gene","rank")

wm983b_ranks <- data.frame(cbind(wm983b_resistant_genes,1:length(wm983b_resistant_genes)))
colnames(wm983b_ranks) <- c("gene","rank")

watermelon_ranks <- data.frame(cbind(watermelon_resistant_genes,1:length(watermelon_resistant_genes)))
colnames(watermelon_ranks) <- c("gene","rank")


temp <- merge(raj_breast_ranks,raj_melanoma_ranks, by="gene",all=T)

temp <- merge(temp,wm9_ranks, by="gene",all=T)

temp <- merge(temp,wm983b_ranks, by="gene",all=T)

ranks_df <- merge(temp,watermelon_ranks, by="gene",all=T)

colnames(ranks_df) <- c("gene","breast","melanoma","wm9","wm983b","watermelon")

ranks_df[is.na(ranks_df)] <- -1

ranks_df <- ranks_df %>% 
  column_to_rownames("gene")

ranks_df <- data.matrix(ranks_df)

row_total_ranks <- rowSums(ranks_df)
names(row_total_ranks) <- rownames(ranks_df)


row_total_ranks <- sort(row_total_ranks, decreasing = T)

ranked_resistance_signature <- list("ranked_resistance_signature" = names(row_total_ranks)[1:200])

saveRDS(ranked_resistance_signature, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ranked_resistance_signature.rds")

















common_resistance_genes <- intersect(intersect(wm9_resistant_markers_refined,wm983b_resistant_markers_refined),intersect(breast_resistant_genes,melanoma_resistant_genes))

common_resistance_signature <- list("common_resistance_signature" = common_resistance_genes)

# saveRDS(common_resistance_signature, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/common_resistance_signature.rds")





write.table(ranked_resistance_signature$ranked_resistance_signature, "test.txt", row.names = F, quote = F, col.names = F)

watermelon_de <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/watermelon_time_point_de.rds")


watermelon_day14_genes <- watermelon_de %>% 
  filter(p_val_adj < 0.05 & cluster == 14 & avg_log2FC > 0) %>% 
  arrange(desc(avg_log2FC)) %>% 
  pull(gene)


all_resistant_genes <- intersect(common_resistance_genes, watermelon_day14_genes)

all_resistant_genes <- list("raj_watermelon_resistance_signature" = all_resistant_genes)

saveRDS(all_resistant_genes, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_watermelon_resistance_signature.rds")


source("source/cole_functions.R")

mat <- calc_jaccard_matrix(watermelon_time_point_de,watermelon_time_point_de)

Heatmap(mat)

length(intersect(melanoma_resistant_genes,wm9_resistant_markers_refined))/length(union(melanoma_resistant_genes,wm9_resistant_markers_refined))

