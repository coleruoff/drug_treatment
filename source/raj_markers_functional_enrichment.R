source("source/cole_functions.R")

raj_markers_genesets <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_melanoma_resistance_refined_cluster_markers.rds")


raj_markers_genesets_melanoma <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_melanoma_resistance_refined_cluster_markers.rds")

raj_markers_genesets_breast <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_breast_resistance_refined_cluster_markers.rds")


Heatmap(calc_jaccard_matrix(raj_markers_genesets_breast,raj_markers_genesets_melanoma), cluster_rows = F,cluster_columns = F)


dotplots <- list()
for(curr_raj_geneset in 1:length(raj_markers_genesets)){
  
  cat(curr_raj_geneset, "\n")
  genes_to_use <- raj_markers_genesets[[curr_raj_geneset]]
  
  go_enrich <- enrichGO(gene = genes_to_use,
                        OrgDb = "org.Hs.eg.db",
                        keyType = 'SYMBOL',
                        readable = T,
                        ont = "ALL",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.10)
  
  curr_cluster_name <- names(raj_markers_genesets)[curr_raj_geneset]
  functional_enrichment_plot_title <- str_to_title(paste(strsplit(curr_cluster_name, "_")[[1]][2], strsplit(curr_cluster_name, "_")[[1]][3],sep = " "))
  
  
  dotplots <- append(dotplots,list(dotplot(go_enrich)+
                                     ggtitle(functional_enrichment_plot_title)))
}


dotplots[[10]]
plt <- egg::ggarrange(plots = dotplots, nrow=2)

main_title <- "Functional Enrichment of Raj Resistance Cluster Signatures"
plt <- annotate_figure(plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))

plot(plt)




genes_to_use <- raj_markers_genesets[[1]]

go_enrich <- enrichGO(gene = genes_to_use,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)

dotplot(go_enrich)