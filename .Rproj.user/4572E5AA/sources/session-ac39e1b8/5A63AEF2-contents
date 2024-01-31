source("source/cole_functions.R")
library(egg)
library(ggpubr)
require(grid)
set.seed(42)

################################################################################
# GSEA Heatmap
################################################################################

################################################################################
#GENESETS 1

curr_cell_line <- "A549"

emergent <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/trapnell_cluster_groups/emergent_clusters.rds")
emergent <- emergent[[curr_cell_line]]

non_emergent <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/trapnell_cluster_groups/nonemergent_clusters.rds")
non_emergent <- non_emergent[[curr_cell_line]]

trapnell_fc_results <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/fold_change_results/", curr_cell_line, "_nonemergent_vs_all_nonemergent_fc.rds"))

genesets_with_ranks <- list()
for(curr_fc_res in trapnell_fc_results){
  
  ranks <- curr_fc_res$avg_log2FC
  
  names(ranks) <- rownames(curr_fc_res)
  
  ranks <- sort(ranks, decreasing = T)
  
  genesets_with_ranks <- append(genesets_with_ranks, list(ranks))
  
}

names(genesets_with_ranks) <- paste0(curr_cell_line,"_cluster_", non_emergent)

names(trapnell_fc_results)

################################################################################
#GENESETS 2

raj_de_genesets <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_breast_clusters_vs_control_de_genesets_500.rds")

raj_cluster_markers <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_breast_cluster_de_markers.rds")

genesets2 <- raj_cluster_markers

#Remove genes from genesets2 that are not in the ranked gene list
for(i in 1:length(genesets2)){
  genesets2[[i]] <- genesets2[[i]][genesets2[[i]] %in% names(genesets_with_ranks[[1]])]
}

################################################################################

result <- create_GSEA_matrix(genesets_with_ranks, genesets2)

gsea_results <- result[[2]]

result <- result[[1]]


plot_title <- paste0("Raj Resistant Cluster DE Genes Enrichment in ", curr_cell_line, " Non-Emergent Clusters")
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# png(filename = paste0("output/figures/trapnell_raj_figures/",curr_cell_line,"_nonemergent_raj_GSEA_heatmap.png"),
#     height = 1000, width = 1600, res = 100)



plot_pretty_heatmap(result, plot_title,legend_title = "NES", col_fun = col_fun)

dev.off()

gsea_results[[13]]$leadingEdge

unique(unlist(gsea_results[[13]]$leadingEdge))


find_consensus_geneset(gsea_results[[13]]$leadingEdge, 3)


Heatmap(calc_jaccard_matrix(raj_de_genesets,raj_de_genesets))


genes_to_use <- gsea_results[[13]]$leadingEdge[[3]]



go_enrich <- enrichGO(gene = genes_to_use,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)

dotplot(go_enrich)




