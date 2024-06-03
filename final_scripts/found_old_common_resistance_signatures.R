
watermelon_de_new <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/de_results/watermelon_time_point_de.rds")
watermelon_de_old <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/watermelon_time_point_de.rds")


watermelon_day14_genes <- watermelon_de %>% 
  filter(p_val_adj < 0.05 & cluster == 14 & avg_log2FC > 0) %>% 
  arrange(desc(avg_log2FC)) %>% 
  pull(gene)

sort(all_resistant_genes$raj_watermelon_resistance_signature)
all_resistant_genes <- intersect(common_resistance_genes, watermelon_day14_genes)

all_resistant_genes <- list("raj_watermelon_resistance_signature" = all_resistant_genes)




###############################################################################
sum(real_common_geneset %in% common_resistance_genes)
sum(real_common_geneset %in% all_resistant_genes$raj_watermelon_resistance_signature)



all_resistant_genesets_og <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/all_resistant_genesets.rds")
real_common_geneset <- intersect(intersect(all_resistant_genesets_og[[1]],all_resistant_genesets_og[[2]]),intersect(all_resistant_genesets_og[[3]],all_resistant_genesets_og[[4]]))
length(real_common_geneset)



temp <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_common_resistance_signature.rds")

all_resistant_genesets <- list(wm9_resistant_genes_refined,wm983b_resistant_genes_refined,breast_resistant_genes,melanoma_resistant_genes,watermelon_day14_genes)
names(all_resistant_genesets) <- c("shaffer_2017_w9","shaffer_2017_w983b","goyal_2023_breast","goyal_2023_melanoma","oren_2021_pc9")

common_resistance_genes <- intersect(intersect(wm9_resistant_markers_refined,wm983b_resistant_markers_refined),intersect(breast_resistant_genes,all_resistant_genesets_og[[4]]))

length(common_resistance_genes)

length(intersect(real_common_geneset,common_resistance_genes))

curr_resistance_signature <- all_resistant_genes$raj_watermelon_resistance_signature

og_resistance_signature <- intersect(all_resistant_genesets_og[[5]],real_common_geneset)


sum(og_resistance_signature %in% all_resistant_genes)
sum(og_resistance_signature %in% all_resistant_genes)


################################################################################

breast_de_res_1 <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/raj_breast_resistant_vs_control_de_res.rds")

breast_de_res_2 <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/de_results/raj_breast_vs_control_deseq_result.rds"))


dim(breast_de_res_1)
dim(breast_de_res_2)

all.equal(as.data.frame(breast_de_res_1), as.data.frame(breast_de_res_2))


melanoma_de_res_1 <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/raj_melanoma_resistant_vs_control_de_res.rds")

melanoma_de_res_2 <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/de_results/raj_melanoma_vs_control_deseq_result.rds"))


potential_genes <- data.frame(melanoma_de_res_1) %>% 
  filter(padj < 0.05 & log2FoldChange > 0) %>% 
  rownames_to_column() %>% 
  pull(rowname) 



intersect(potential_genes, all_resistant_genesets_og[[4]])

dim(melanoma_de_res_1)
dim(melanoma_de_res_2)

as.data.frame(melanoma_de_res_2) %>% 
  arrange(desc(log2FoldChange)) %>% 
  rownames_to_column() %>% 
  filter(rowname %in% all_resistant_genesets_og[[4]]) %>% 
  pull(log2FoldChange) %>% 
  min()


temp <- as.data.frame(melanoma_de_res_1) %>% 
  rownames_to_column() %>% 
  mutate("in_old" = rowname %in% all_resistant_genesets_og[[4]]) 

as.data.frame(melanoma_de_res_1) %>% 
  rownames_to_column() %>% 
  filter(log2FoldChange > 0 & padj < 0.05)
  





as.data.frame(melanoma_de_res_2) %>% 
  filter(log2FoldChange > 0 & padj < 0.05) %>% 
  nrow()




melanoma_resistant_genes <- data.frame(melanoma_de_res_1) %>% 
  filter(padj < 0.05 & log2FoldChange > 0) %>% 
  rownames_to_column() %>% 
  pull(rowname)


melanoma_resistant_genes <- data.frame(melanoma_de_res_2) %>% 
  filter(padj < 0.05 & log2FoldChange > 0) %>% 
  rownames_to_column() %>% 
  pull(rowname)


all.equal(as.data.frame(melanoma_de_res_1), as.data.frame(melanoma_de_res_2))





melanoma_de_res <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/de_results/raj_melanoma_vs_control_deseq_result.rds"))
melanoma_de_res <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/raj_melanoma_resistant_vs_control_de_res.rds")

melanoma_resistant_genes <- data.frame(melanoma_de_res) %>% 
  filter(padj < 0.05 & log2FoldChange > 3) %>% 
  rownames_to_column() %>% 
  pull(rowname)

length(melanoma_resistant_genes)
