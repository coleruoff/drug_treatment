library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
source("source/cole_functions.R")

genesets_characterization <- function(genesets_to_use, universe_to_use, num_pathways=10){
  m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, human_gene_symbol)
  # hallmark_uni <- unique(c(unlist(all_signatures), m_t2g$human_gene_symbol))
  
  mp_t2g <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ith_meta_programs_t2g.rds")
  
  # mp_uni <- unique(c(unlist(all_signatures), mp_t2g$human_gene_symbol))
  
  hallmarks_dotplots <- list()
  mps_dotplots <- list()
  go_dotplots <- list()
  go_results <- list()
  
  for(i in 1:length(genesets_to_use)){
    
    curr_cluster <- names(genesets_to_use)[i]
    cat(curr_cluster, "\n")
    
    curr_geneset <- genesets_to_use[[i]]
    
    
    # Hallmarks Enrichment
    hallmark_enrichment_results <- enricher(curr_geneset, TERM2GENE=m_t2g,
                                            universe = universe_to_use)
    
    # MPs Enrichment
    mp_enrichment_results <- enricher(curr_geneset, TERM2GENE=mp_t2g,
                                      universe = universe_to_use)
    
    # GO Functional Enrichment
    ego <- enrichGO(gene          = curr_geneset,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = .01,
                    qvalueCutoff  = .05,
                    readable      = TRUE,
                    keyType = "SYMBOL",
                    universe = universe_to_use)
    
    go_results <- append(go_results,list(ego))
    
    if(nrow(hallmark_enrichment_results) > 0){
      p <- dotplot(hallmark_enrichment_results,
                   showCategory = num_pathways) + ggtitle(paste0("", curr_cluster))
      hallmarks_dotplots <- append(hallmarks_dotplots, list(p))
    }
    if(!is.null(mp_enrichment_results) && nrow(mp_enrichment_results) > 0){
      p <- dotplot(mp_enrichment_results,
                   showCategory = num_pathways)+ggtitle(paste0("", curr_cluster))
      mps_dotplots <- append(mps_dotplots, list(p))
    }
    if(nrow(ego) > 0){
      p <- dotplot(ego,
                   showCategory = num_pathways) + ggtitle(paste0("", curr_cluster))
      go_dotplots <- append(go_dotplots, list(p))
    }
  }
  
  
  ret_list <- list(hallmarks_dotplots, mps_dotplots, go_dotplots,go_results)
  names(ret_list) <- c("hallmarks_dotplots", "mps_dotplots", "go_dotplots","go_results")
  
  return(ret_list) 
}

cell_lines <- c("A549","K562","MCF7")

# cell_line_universes <- list()
# for(curr_cell_line in cell_lines){
#   
#   cat(curr_cell_line, "\n")
#   
#   #find genes expressed in > 1% of cells
#   data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/",curr_cell_line,"_processed_filtered.rds"))
#   gene_universe <- names(apply(data@assays$RNA@data, 1, FUN = function(x) sum(x>0) >.01*ncol(data)))
#   
#   
#   cell_line_universes <- append(cell_line_universes, list(gene_universe))
# }
# 
# names(cell_line_universes) <- cell_lines
# 
# gene_universe <- intersect(cell_line_universes[[1]], cell_line_universes[[2]])
# gene_universe <- intersect(gene_universe, cell_line_universes[[3]])
# 
# saveRDS(gene_universe, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cell_line_gene_universe_intersection.rds")
# saveRDS(cell_line_universes, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cell_line_universes.rds")

gene_universe_intersection <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cell_line_gene_universe_intersection.rds")
cell_line_universes <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cell_line_universes.rds")

################################################################################
# Plotting for RAC type 1 superclusters signatures
type1_supercluster_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/type1_supercluster_signatures.rds")

names1 <- gsub("type1_", "", names(type1_supercluster_signatures))
names2 <- gsub("_", " ", names1)
names3 <- gsub("supercluster", "Supercluster ", names2)
names(type1_supercluster_signatures) <- gsub("signature", "Signature ", names3)

all_dotplots <- genesets_characterization(type1_supercluster_signatures, universe_to_use = gene_universe_intersection, num_pathways=20)

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/supercluster_enrichment/type1_superclusters_hallmarks.png"), width = 1000,height = 800)
hallmarks_plt <- egg::ggarrange(plots = all_dotplots$hallmarks_dotplots, ncol = 2)
main_title <- paste0("Cancer Hallmarks Enrichment of RAC Type 1 Superclusters Signatures")
hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(hallmarks_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/supercluster_enrichment/type1_superclusters_mps.png"), width = 1000,height = 800)
mps_plt <- egg::ggarrange(plots = all_dotplots$mps_dotplots, ncol = 2)
main_title <- paste0("ITH Meta-programs Enrichment of RAC Type 1 Superclusters Signatures")
mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(mps_plt)
dev.off()

png(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/signature_enrichment_dotplots/supercluster_enrichment/type1_superclusters_functional_enrichment.png"), width = 1000,height = 800)
go_plt <- egg::ggarrange(plots = all_dotplots$go_dotplots, ncol = 2)
main_title <- paste0("GO Functional Enrichment of RAC Type 1 Superclusters Signatures")
go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
plot(go_plt)
dev.off()

