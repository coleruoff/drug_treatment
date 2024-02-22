library(msigdbr)
library(org.Hs.eg.db)
library(ggpubr)

genesets_characterization <- function(genesets_to_use){
  m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, human_gene_symbol)
  # hallmark_uni <- unique(c(unlist(all_signatures), m_t2g$human_gene_symbol))
  
  mp_t2g <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ith_meta_programs_t2g.rds")
  
  # mp_uni <- unique(c(unlist(all_signatures), mp_t2g$human_gene_symbol))
  
  hallmarks_dotplots <- list()
  mps_dotplots <- list()
  go_dotplots <- list()
  
  for(i in 1:length(genesets_to_use)){
    
    curr_cluster <- names(genesets_to_use)[i]
    cat(curr_cluster, "\n")
    
    curr_geneset <- genesets_to_use[[i]]
    
    
    # Hallmarks Enrichment
    hallmark_enrichment_results <- enricher(curr_geneset, TERM2GENE=m_t2g) #,
    #universe = hallmark_uni)
    
    # MPs Enrichment
    mp_enrichment_results <- enricher(curr_geneset, TERM2GENE=mp_t2g)#,
    # universe = mp_uni)
    
    # GO Functional Enrichment
    ego <- enrichGO(gene          = curr_geneset,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = .01,
                    qvalueCutoff  = .05,
                    readable      = TRUE,
                    keyType = "SYMBOL")
    
    if(nrow(hallmark_enrichment_results) > 0){
      p <- dotplot(hallmark_enrichment_results) + ggtitle(paste0("", curr_cluster))
      hallmarks_dotplots <- append(hallmarks_dotplots, list(p))
    }
    if(!is.null(mp_enrichment_results) && nrow(mp_enrichment_results) > 0){
      p <- dotplot(mp_enrichment_results)+ggtitle(paste0("", curr_cluster))
      mps_dotplots <- append(mps_dotplots, list(p))
    }
    if(nrow(ego) > 0){
      p <- dotplot(ego) + ggtitle(paste0("", curr_cluster))
      go_dotplots <- append(go_dotplots, list(p))
    }
  }
  
  
  ret_list <- list(hallmarks_dotplots, mps_dotplots, go_dotplots)
  names(ret_list) <- c("hallmarks_dotplots", "mps_dotplots", "go_dotplots")
  
  return(ret_list) 
}

################################################################################

active_supercluster_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/active_supercluster_signatures.rds")

rac_supercluster_signature <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_signature.rds")

rac_supercluster_consensus_signature <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_consensus_signature.rds")

################################################################################
# RAC supercluster
################################################################################

all_dotplots <- genesets_characterization(rac_supercluster_signature)

png("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/supercluster_figures/rac_supercluster_hallmarks.png",
    width=800,height=800)
hallmarks_plt <- egg::ggarrange(plots = all_dotplots$hallmarks_dotplots, ncol = 2)
main_title <- paste0("Hallmarks Enrichment of Supercluster Signature")
hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 28))
plot(hallmarks_plt)
dev.off()

# png("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/supercluster_figures/rac_supercluster_mps.png",
#     width=1000,height=3000)
# mps_plt <- egg::ggarrange(plots = all_dotplots$mps_dotplots, ncol = 2)
# main_title <- paste0("ITH Meta-programs Enrichment of Supercluster Signature")
# mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 14))
# plot(mps_plt)
# dev.off()

png("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/supercluster_figures/rac_supercluster_go.png",
    width=800,height=800)
go_plt <- egg::ggarrange(plots = all_dotplots$go_dotplots, ncol = 2)
main_title <- paste0("GO Functional Enrichment of Supercluster Signature")
go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 28))
plot(go_plt)
dev.off()

################################################################################
# RAC consensus supercluster
################################################################################

all_dotplots <- genesets_characterization(rac_supercluster_consensus_signature)

png("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/supercluster_figures/rac_consensus_supercluster_hallmarks.png",
    width=800,height=800)
hallmarks_plt <- egg::ggarrange(plots = all_dotplots$hallmarks_dotplots, ncol = 2)
main_title <- paste0("Hallmarks Enrichment of Supercluster Consensus Signature")
hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 24))
plot(hallmarks_plt)
dev.off()

png("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/supercluster_figures/rac_consensus_supercluster_mps.png",
    width=800,height=800)
mps_plt <- egg::ggarrange(plots = all_dotplots$mps_dotplots, ncol = 2)
main_title <- paste0("ITH Meta-programs Enrichment of Supercluster Consensus Signature")
mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 24))
plot(mps_plt)
dev.off()

png("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/supercluster_figures/rac_consensus_supercluster_go.png",
    width=800,height=800)
go_plt <- egg::ggarrange(plots = all_dotplots$go_dotplots, ncol = 2)
main_title <- paste0("GO Functional Enrichment of Supercluster Consensus Signature")
go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 24))
plot(go_plt)
dev.off()


################################################################################

################################################################################
# Active superclusters
################################################################################

all_dotplots <- genesets_characterization(active_supercluster_signatures)

png("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/supercluster_figures/active_superclusters_hallmarks.png",
    width=1000,height=1000)
hallmarks_plt <- egg::ggarrange(plots = all_dotplots$hallmarks_dotplots, ncol = 2)
main_title <- paste0("Hallmarks Enrichment of Resistance Active Supercluster Signatures")
hallmarks_plt <- annotate_figure(hallmarks_plt, top = text_grob(main_title, color = "black", face = "bold", size = 28))
plot(hallmarks_plt)
dev.off()

png("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/supercluster_figures/active_superclusters_mps.png",
    width=1000,height=1000)
mps_plt <- egg::ggarrange(plots = all_dotplots$mps_dotplots, ncol = 2)
main_title <- paste0("ITH Meta-programs Enrichment of Resistance Active Supercluster Signatures")
mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 24))
plot(mps_plt)
dev.off()

png("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/supercluster_figures/active_superclusters_go.png",
    width=1000,height=1400)
go_plt <- egg::ggarrange(plots = all_dotplots$go_dotplots, ncol = 2)
main_title <- paste0("GO Functional Enrichment of Resistance Active Supercluster Signatures")
go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 28))
plot(go_plt)
dev.off()

