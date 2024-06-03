library(tidyverse)


# THIS SCRIPT IS NOT USED IN THE PROJECT. USED FOR EXPLORATION


nodrug_vs_resistant_de <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/resNoDrugVsResistant.csv", row.names = 1)

nodrug_vs_drug_de <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/resNoDrugVsDrug.csv", row.names = 1)

drug_vs_resistant_de <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/resDrugVsResistant.csv", row.names = 1)

resistant_markers <- nodrug_vs_resistant_de %>% 
  filter(padj < 0.05 & log2FoldChange > 2) %>% 
  arrange(desc(log2FoldChange)) %>% 
  pull(GeneSymbol)

early_drug_markers <- nodrug_vs_drug_de %>% 
  filter(padj < 0.05 & log2FoldChange > 2) %>% 
  arrange(desc(log2FoldChange)) %>% 
  pull(GeneSymbol)

# I dont know which condition is the base condition
idk_de <- drug_vs_resistant_de %>% 
  filter(padj < 0.05 & log2FoldChange > 1.4) %>% 
  arrange(desc(log2FoldChange)) %>% 
  pull(GeneSymbol)


resistant_markers

resistant_markers_refined <- resistant_markers[!resistant_markers %in% early_drug_markers]

early_drug_markers_refined <- early_drug_markers[!early_drug_markers %in% resistant_markers]


raj_resistant_2017 <- list("raj_resistant_2017" = resistant_markers_refined)

raj_early_response_2017 <- list("raj_early_response_2017" = early_drug_markers_refined)



raj_melanoma <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_melanoma_resistant_vs_control_de_signature_500.rds")


sum(raj_melanoma %in% resistant_markers)/length(raj_melanoma)


length(intersect(resistant_markers, early_drug_markers))/(length(union(resistant_markers, early_drug_markers)))
