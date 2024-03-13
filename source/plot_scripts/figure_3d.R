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
  
  specifc_mps <- c("MP39 Metal-response","MP31 Alveolar","MP29 NPC/OPC","MP28 Oligo normal","MP27 Oligo Progenitor","MP38 Glutathione","MP41 Unassigned","MP35 Hemato-related-I","MP37 Hemato-related-II","MP32 Skin-pigmentation","MP36 IG","MP16 MES (glioma)","MP15 EMT IV")
  
  mp_t2g <- mp_t2g %>% 
    filter(!gs_name %in% specifc_mps)
  
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
      p <- barplot(hallmark_enrichment_results,
                   showCategory = num_pathways, font.size=20) + 
        ggtitle(paste0("", curr_cluster))+
        theme(plot.title = element_text(size=20))
      
      hallmarks_dotplots <- append(hallmarks_dotplots, list(p))
    }
    if(!is.null(mp_enrichment_results) && nrow(mp_enrichment_results) > 0){
      p <- barplot(mp_enrichment_results,
                   showCategory = num_pathways, font.size=20)+
        ggtitle(paste0("", curr_cluster))+
        theme(plot.title = element_text(size=20))
      
      mps_dotplots <- append(mps_dotplots, list(p))
    }
    if(nrow(ego) > 0){
      p <- barplot(ego,
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

gene_universe_intersection <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cell_line_gene_universe_intersection.rds")
cell_line_universes <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cell_line_universes.rds")

################################################################################
# Plotting for RAC superclusters signatures
supercluster_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_signatures.rds")

names1 <- gsub("type1_", "", names(supercluster_signatures))
names2 <- gsub("_", " ", names1)
names3 <- gsub("supercluster", "Supercluster ", names2)
names(supercluster_signatures) <- gsub("signature", "Signature ", names3)

all_dotplots <- genesets_characterization(supercluster_signatures, universe_to_use = gene_universe_intersection)


hallmarks_plt <- ggarrange(plotlist = all_dotplots$hallmarks_dotplots, ncol = 1, common.legend = T, legend=c("right"))
main_title <- paste0("\nCancer Hallmarks")
hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 30))


mps_plt <- ggarrange(plotlist = all_dotplots$mps_dotplots, ncol = 1, common.legend = T,legend=c("right"))
main_title <- paste0("\nITH Meta-programs")
mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 30))


go_plt <- ggarrange(plotlist = all_dotplots$go_dotplots, ncol = 1, common.legend = T,legend=c("right"))
main_title <- paste0("\nGO Pathways")
go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 30))


plots <- list(hallmarks_plt,mps_plt, go_plt)


# png(paste0("/data/ruoffcj/projects/drug_treatment/final_figures/figure_3d.png"),
#     width=30, height=24, units= "in", res = 300)



figure <- ggarrange(plotlist = plots, ncol=3, nrow=1, common.legend = T,legend=c("right"))

p <- annotate_figure(figure, left = text_grob("", rot = 90, vjust = 1, size=35, face="bold"),
                     bottom = text_grob("", size=35, face="bold"),
                     top=text_grob("RAC Superclusters Signatures Enrichment", size=40, face="bold"))

print(p)

dev.off()



