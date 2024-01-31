source("source/cole_functions.R")
library(clusterProfiler)
library(ggpubr)
library(tidyverse)

curr_cell_line <- "MCF7"

################################################################################
# Functional enrichment of all genes upregulated in all resistance active cells

all_de_result <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line,"_all_resistant_cells_de.rds"))

genes_to_use <- all_de_result %>% 
  filter(cluster == "resistant" & p_val_adj < .05 & avg_log2FC > .5) %>% 
  arrange(desc(avg_log2FC)) %>% 
  pull(gene)

go_enrich <- enrichGO(gene = genes_to_use,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)

dotplot(go_enrich, showCategory=20) + ggtitle(paste0(curr_cell_line, " Functional Enrichment of Genes Upregulated in All Resistance Active Cells"))

broad_go_terms <- dotplot(go_enrich, showCategory=30)


broad_go_terms <- broad_go_terms$data$ID


#################################################################################
# Create active component signatures

curr_cell_line <- "A549"

de_result <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line,"_resistant_clusters_de.rds"))

# Intracluster DE Upregulated genes
intracluster_active_de_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/",curr_cell_line, "intracluster_active_de_signatures.rds"))

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- RACs[[curr_cell_line]]

active_signatures <- list()

for(curr_cluster in clusters_of_interest){
  
  cat(curr_cluster, "\n")
  
  inactive_signature <- de_result %>% 
    filter(cluster == curr_cluster & p_val_adj < .05 & avg_log2FC > .5) %>% 
    arrange(desc(avg_log2FC)) %>% 
    pull(gene)
  
  active_signature <- de_result %>% 
    filter(cluster == paste0(curr_cluster,"_resistant") & p_val_adj < .05 & avg_log2FC > .5) %>% 
    arrange(desc(avg_log2FC)) %>% 
    pull(gene)
  
  #Find shared genes between active and inactive componenents
  # shared_genes <- intersect(active_signature, inactive_signature)
  # 
  # #Remove shared genes
  # active_signature <- active_signature[!active_signature %in% shared_genes]
  # 
  # #Add genes from DE between active and inactive within this cluster
  # active_signature <- unique(append(active_signature, intracluster_active_de_signatures[[paste0(curr_cluster,"_active_signature")]]))
  
  active_signatures <- append(active_signatures, list(active_signature))
  
}

names(active_signatures) <- paste0(clusters_of_interest, "_active_signature")

saveRDS(active_signatures, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/",curr_cell_line,"_active_signatures.rds"))
#################################################################################

dotplots <- list()

for(i in 1:length(active_signatures)){
  cat(i, "\n")
  
  genes_to_use <- active_signatures[[i]]
  
  go_enrich <- enrichGO(gene = genes_to_use,
                        OrgDb = "org.Hs.eg.db",
                        keyType = 'SYMBOL',
                        readable = T,
                        ont = "ALL",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.10)
  
  # go_enrich@result <- go_enrich@result[!go_enrich@result$ID %in% broad_go_terms,]
  
  if(nrow(go_enrich) > 0){
    p <- dotplot(go_enrich) + ggtitle(names(active_signatures)[i])
    
    
    dotplots <- append(dotplots, list(p))
  }
  
  
}


plt <- egg::ggarrange(plots = dotplots, ncol = 4)

main_title <- paste0(curr_cell_line, " Functional Enrichment of Resistant Active Clusters Signatures")
plt <- annotate_figure(plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))

plot(plt)


Heatmap(calc_jaccard_matrix(active_signatures,active_signatures))
