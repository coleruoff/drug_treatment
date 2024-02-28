setwd("/data/ruoffcj/projects/drug_treatment/")
library(tidyverse)
library(clusterProfiler)
library(ggpubr)
library(msigdbr)
library(org.Hs.eg.db)

yeast_human_orthologs_up <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/yeast_human_orthologs_up.rds")

supercluster_signature <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_signature.rds")
supercluster_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_signatures.rds")


m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, human_gene_symbol)
# hallmark_uni <- unique(c(unlist(all_signatures), m_t2g$human_gene_symbol))

mp_t2g <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ith_meta_programs_t2g.rds")

gene_universe_intersection <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cell_line_gene_universe_intersection.rds")

# yeast

i <- 1
for(i in 1:2){
  curr_geneset <- intersect(supercluster_signatures[[i]],yeast_human_orthologs_up$yeast_human_orthologs_up)
  
  # e coli
  # curr_geneset <- intersect(supercluster_signature$supercluster_signature,ecoli_human_orthologs_up$ecoli_human_orthologs_up)
  
  
  universe_to_use <- gene_universe_intersection
  
  # Hallmarks Enrichment
  # hallmark_enrichment_results <- enricher(curr_geneset, TERM2GENE=m_t2g,
  #                                         universe = universe_to_use)
  # 
  # 
  # hallmark_plot <- barplot(hallmark_enrichment_results)
  
  
  # MPs Enrichment
  mp_enrichment_results <- enricher(curr_geneset, TERM2GENE=mp_t2g,
                                    universe = universe_to_use)
  
  
  
  
  mp_plot <- barplot(mp_enrichment_results,font.size = 20)
  
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
  
  
  
  go_plot <- barplot(ego,font.size = 20) 
  plot_list <- list(mp_plot,go_plot)
}
curr_geneset <- intersect(supercluster_signature$supercluster_signature,yeast_human_orthologs_up$yeast_human_orthologs_up)

# e coli
# curr_geneset <- intersect(supercluster_signature$supercluster_signature,ecoli_human_orthologs_up$ecoli_human_orthologs_up)


universe_to_use <- gene_universe_intersection

# Hallmarks Enrichment
# hallmark_enrichment_results <- enricher(curr_geneset, TERM2GENE=m_t2g,
#                                         universe = universe_to_use)
# 
# 
# hallmark_plot <- barplot(hallmark_enrichment_results)


# MPs Enrichment
mp_enrichment_results <- enricher(curr_geneset, TERM2GENE=mp_t2g,
                                  universe = universe_to_use)




mp_plot <- barplot(mp_enrichment_results,font.size = 20)

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



go_plot <- barplot(ego,font.size = 20) 
plot_list <- list(mp_plot,go_plot)

figure <- ggarrange(plotlist = plot_list, ncol=2, nrow=1, common.legend = T,legend=c("right"))

p <- annotate_figure(figure, left = text_grob("", rot = 90, vjust = 1, size=35, face="bold"),
                     bottom = text_grob("", size=35, face="bold"),
                     top=text_grob("Enrichment of Shared Genes Between RAC Supercluster\nand Yeast Antifungal Resistance Orthologs", size=40, face="bold"))


png(paste0("/data/ruoffcj/projects/drug_treatment/final_figures/figure_5b.png"),
    width=20, height=12, units= "in", res = 300)

print(p)

dev.off()











