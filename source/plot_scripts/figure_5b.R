setwd("/data/ruoffcj/projects/drug_treatment/")
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
                   showCategory = num_pathways, font.size=20) + 
        ggtitle(paste0("", curr_cluster))+
        theme(plot.title = element_text(size=20))
      
      hallmarks_dotplots <- append(hallmarks_dotplots, list(p))
    }
    if(!is.null(mp_enrichment_results) && nrow(mp_enrichment_results) > 0){
      p <- dotplot(mp_enrichment_results,
                   showCategory = num_pathways, font.size=20)+
        ggtitle(paste0("", curr_cluster))+
        theme(plot.title = element_text(size=20))
      mps_dotplots <- append(mps_dotplots, list(p))
    }
    if(nrow(ego) > 0){
      p <- dotplot(ego,
                   showCategory = num_pathways, font.size=20) + 
        ggtitle(paste0("", curr_cluster))+
        theme(plot.title = element_text(size=20))
      
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
# Plotting for shared genes between RACs signatures and yeast orthologs


type1_supercluster_upregulated_genesets <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/type1_superclusters_upregulated_genesets.rds")
yeast_human_orthologs_up <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/yeast_human_orthologs_up.rds")

shared_supercluster1_genes <- intersect(type1_supercluster_upregulated_genesets[[1]],yeast_human_orthologs_up$yeast_human_orthologs_up)
shared_supercluster2_genes <- intersect(type1_supercluster_upregulated_genesets[[2]],yeast_human_orthologs_up$yeast_human_orthologs_up)

shared_signatures <- list("Supercluster 1" = shared_supercluster1_genes,"Supercluster 2" = shared_supercluster2_genes)


all_dotplots <- list()
for(i in 1:length(shared_signatures)){
  
  curr_signature <- shared_signatures[i]
  
  curr_dotplots <- genesets_characterization(curr_signature, universe_to_use = gene_universe_intersection)
  
  all_dotplots[["hallmarks_dotplots"]] <- append(all_dotplots[["hallmarks_dotplots"]],curr_dotplots[["hallmarks_dotplots"]])
  all_dotplots[["mps_dotplots"]] <- append(all_dotplots[["mps_dotplots"]],curr_dotplots[["mps_dotplots"]])
  all_dotplots[["go_dotplots"]] <- append(all_dotplots[["go_dotplots"]],curr_dotplots[["go_dotplots"]])
  all_dotplots[["go_results"]] <- append(all_dotplots[["go_results"]],curr_dotplots[["go_results"]])
  
}

# saveRDS(shared_genes, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/yeast_RAC_shared_genes.rds")

hallmarks_plt <- ggarrange(plotlist = all_dotplots$hallmarks_dotplots, ncol = 2, common.legend = T, legend=c("right"))
main_title <- paste0("\nCancer Hallmarks")
hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 30))


mps_plt <- ggarrange(plotlist = all_dotplots$mps_dotplots, ncol = 2, common.legend = T,legend=c("right"))
main_title <- paste0("ITH Meta-programs")
mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 30))


go_plt <- ggarrange(plotlist = all_dotplots$go_dotplots, ncol = 2, common.legend = T,legend=c("right"))
main_title <- paste0("GO Pathways")
go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 30))


plots <- list(hallmarks_plt, mps_plt, go_plt)


png(paste0("/data/ruoffcj/projects/drug_treatment/final_figures/figure_5b.png"),
    width=26, height=24, units= "in", res = 300)



figure <- ggarrange(plotlist = plots, nrow=3, common.legend = T,legend=c("right"))

p <- annotate_figure(figure, left = text_grob("", rot = 90, vjust = 1, size=35, face="bold"),
                     bottom = text_grob("", size=35, face="bold"),
                     top=text_grob("Enrichment of Shared Genes Between RAC Type 1 Superclusters Signatures \nand Yeast Antifungal Resistance Orthologs", size=35, face="bold"))

print(p)

dev.off()
