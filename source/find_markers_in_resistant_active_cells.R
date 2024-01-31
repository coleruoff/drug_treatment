source("source/cole_functions.R")

cell_lines <- c("A549","K562","MCF7")

#################################################################################
# Read in Trapnell data
#################################################################################
curr_cell_line <- cell_lines[3]

data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))

num_clusters <- nlevels(data$Cluster)

#################################################################################
# Read in AUCell score results and thresholds
#################################################################################

resistant_cancer_type <- "common"

scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line, "_processed_filtered_raj_",resistant_cancer_type,"_resistance_signature_aucell_scores.rds"))

threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line,"_processed_filtered_raj_",resistant_cancer_type,"_resistance_signature_aucell_thesholds.rds"))

all_active_cells <- rownames(scores)[scores > threshold$threshold]



data <- AddMetaData(data, metadata=ifelse(rownames(data@meta.data) %in% all_active_cells, 1, 0), col.name = "active_resistance")


Idents(data) <- data$active_resistance


resistance_de_markers <- FindMarkers(data, ident.1 = 1)



resistant_active_markers <- resistance_de_markers %>% 
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>% 
  rownames()

#################################################################################
# 
#################################################################################
genes_to_use <- resistant_active_markers

go_enrich <- enrichGO(gene = genes_to_use,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)


dotplot(go_enrich,showCategory=30)+
  ggtitle(paste0("Functional Enrichment of Markers of Resistant Active ", curr_cell_line, " Trapnell Cells"))


#################################################################################
#
#################################################################################
raj_common_resistance_signature <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_common_resistance_signature.rds")

genes_to_use <- intersect(raj_common_resistance_signature,resistant_active_markers)

go_enrich <- enrichGO(gene = genes_to_use,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)


dotplot(go_enrich,showCategory=30)+
  ggtitle(paste0("Functional Enrichment of Markers of Resistant Active ", curr_cell_line, " Trapnell Cells and Raj Common Resistance Signature Intersection"))



