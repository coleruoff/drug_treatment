################################################################################
# Enrichment of yeast human orthologs in HUMAN
################################################################################

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

################################################################################

yeast_human_orthologs_up <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/yeast_human_orthologs_up.rds")

names(yeast_human_orthologs_up) <- NULL

all_dotplots <- genesets_characterization(yeast_human_orthologs_up)


hallmarks_plt <- ggarrange(plotlist = all_dotplots$hallmarks_dotplots, ncol = 1, common.legend = T, legend=c("right"))
main_title <- paste0("\nCancer Hallmarks")
hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 30))


mps_plt <- ggarrange(plotlist = all_dotplots$mps_dotplots, ncol = 1, common.legend = T,legend=c("right"))
main_title <- paste0("ITH Meta-programs")
mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 30))


go_plt <- ggarrange(plotlist = all_dotplots$go_dotplots, ncol = 1, common.legend = T,legend=c("right"))
main_title <- paste0("GO Pathways")
go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 30))


plots <- list(hallmarks_plt, mps_plt, go_plt)


png(paste0("/data/ruoffcj/projects/drug_treatment/figures/yeast_orthologs_enrichment_plots.png"),
    width=26, height=14, units= "in", res = 300)


figure <- ggarrange(plotlist = plots, ncol=3, common.legend = T,legend=c("right"))

p <- annotate_figure(figure, left = text_grob("", rot = 90, vjust = 1, size=35, face="bold"),
                     bottom = text_grob("", size=35, face="bold"),
                     top=text_grob("Enrichment of Yeast Antifungal Resistance Orthologs\n", size=35, face="bold"))

print(p)

dev.off()

