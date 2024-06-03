library(msigdbr)
library(tidyverse)
library(clusterProfiler)
library(ggpubr)


GSEA_dotplots <- function(ranked_genelists, t2g_df){
  
  all_GSEA_results <- list()
  
  for(i in names(ranked_genelists)){
    geneList <- ranked_genelists[[i]]
    
    em <- GSEA(geneList, TERM2GENE = t2g_df)
    
    all_GSEA_results <- append(all_GSEA_results, em)
  }
  
  names(all_GSEA_results) <- names(ranked_genelists)
  
  ################################################################################
  #Plot GSEA results
  
  dotplots <- list()
  for(i in names(ranked_genelists)){
    
    if(nrow(all_GSEA_results[[i]]) > 0){
      
      p <- dotplot(all_GSEA_results[[i]], showCategory=30) + ggtitle(i)
      
      dotplots <- append(dotplots, list(p))
    }
    
  }
  
  
  return(dotplots)
  
}

################################################################################
# Read in ranked genelist

curr_cell_line <- "A549"

trapnell_fc_results <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/fold_change_results/", curr_cell_line, "_all_clusters_fc.rds"))

names_trim1 <- gsub(paste0(curr_cell_line,"_cluster"), "", names(trapnell_fc_results))
names_trim2 <- gsub("_fc_results", "",names_trim1)

all_clusters <- names(trapnell_fc_results)

genesets_with_ranks <- list()
for(curr_fc_res in trapnell_fc_results){
  
  ranks <- curr_fc_res$avg_log2FC
  
  names(ranks) <- rownames(curr_fc_res)
  
  ranks <- sort(ranks, decreasing = T)
  
  genesets_with_ranks <- append(genesets_with_ranks, list(ranks))
  
}

names(genesets_with_ranks) <- paste0(curr_cell_line,"_cluster_", all_clusters)

names(genesets_with_ranks) <- names(trapnell_fc_results)


pre_treatment_cluster_ranked_markers <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_pre_treatment_cluster_markers.rds"))

post_treatment_cluster_ranked_markers <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_post_treatment_cluster_markers.rds"))

# jaccard_mat <- calc_jaccard_matrix(pre_treatment_cluster_ranked_markers, post_treatment_cluster_ranked_markers)

################################################################################
#Set up term to gene df

msigdbr_show_species()

m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame


hallmark_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

head(hallmark_t2g)


################################################################################
# Run GSEA

dotplots <- GSEA_dotplots(genesets_with_ranks, hallmark_t2g)

plt <- egg::ggarrange(plots = dotplots, nrow=5)

main_title <- paste0(curr_cell_line, " Functional Enrichment of Trapnell Pre-Treatment Clusters")
plt <- annotate_figure(plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))

dev.off()
plot(plt)

names(genesets_with_ranks)[1]




dotplots[[3]]



