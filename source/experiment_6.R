library(forcats)
library(org.Hs.eg.db)
library(clusterProfiler)
library(msigdbr)
library(ggpubr)


m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, human_gene_symbol)
# hallmark_uni <- unique(c(unlist(all_signatures), m_t2g$human_gene_symbol))

mp_t2g <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ith_meta_programs_t2g.rds")

specifc_mps <- c("MP39 Metal-response","MP31 Alveolar","MP29 NPC/OPC","MP28 Oligo normal","MP27 Oligo Progenitor","MP38 Glutathione","MP41 Unassigned","MP35 Hemato-related-I","MP37 Hemato-related-II","MP32 Skin-pigmentation","MP36 IG","MP16 MES (glioma)","MP15 EMT IV")

mp_t2g <- mp_t2g %>% 
  filter(!gs_name %in% specifc_mps)


gene_universe_intersection <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cell_line_gene_universe_intersection.rds")


supercluster_signature <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_signature.rds")
supercluster_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_signatures.rds")

curr_geneset <- supercluster_signature$supercluster_signature
curr_geneset <- supercluster_signatures$supercluster2_signature
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
                     top=text_grob("RAC Supercluster Signature Enrichment", size=40, face="bold"))


png(paste0("/data/ruoffcj/projects/drug_treatment/final_figures/figure_3b.png"),
    width=20, height=12, units= "in", res = 300)

print(p)

dev.off()


write.table(sort(curr_geneset),"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/test.txt",
            quote = F, row.names = F,col.names = F)

