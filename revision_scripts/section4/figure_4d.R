args = commandArgs(trailingOnly=TRUE)
# dataDirectory <- paste0(args[1],"final_data/")
# plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("revision_scripts/drug_treatment_functions.R")
set.seed(42)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_figures/"

################################################################################

supercluster_signatures <- readRDS(paste0(dataDirectory, "genesets/supercluster_up_signatures.rds"))

ecoli_human_orthologs_up <- readRDS(paste0(dataDirectory, "genesets/ecoli_human_orthologs_up.rds"))

overlap_signatures <- list()

for(i in c(2)){
  curr_geneset <- intersect(supercluster_signatures[[i]], ecoli_human_orthologs_up$ecoli_human_orthologs_up)
  
  overlap_signatures <- append(overlap_signatures,list(curr_geneset))
}


names(overlap_signatures) <- c("Supercluster 2")

final_plot <- plot_enrichment_ora(overlap_signatures)

jpeg(paste0(plotDirectory,"figure_4d.jpg"), width=250, height = 90, units = "mm", res = 1000)
print(final_plot)
dev.off()







# 
# 
#   all_results <- list()
#   
#   # Hallmarks Enrichment
#   hallmark_enrichment_results <- enricher(curr_geneset, TERM2GENE=m_t2g)
#                                           # universe = universe_to_use)
#   
#   # MPs Enrichment
#   mp_enrichment_results <- enricher(curr_geneset, TERM2GENE=mp_t2g)
#                                     # universe = universe_to_use)
#   
#   # GO Functional Enrichment
#   go_enrichment_results <- enrichGO(gene          = curr_geneset,
#                                     OrgDb         = org.Hs.eg.db,
#                                     ont           = "BP",
#                                     pAdjustMethod = "BH",
#                                     pvalueCutoff  = .01,
#                                     qvalueCutoff  = .05,
#                                     readable      = TRUE,
#                                     keyType = "SYMBOL")
#                                     # universe = universe_to_use)
#   
#   
#   plots <- list()
#   
#   if(!is.null(hallmark_enrichment_results)){
#     p1 <- dotplot(hallmark_enrichment_results, font.size=3)+
#       theme_classic()+
#       theme(axis.text.x = element_text(size = 5),
#             legend.text = element_text(size=3),
#             legend.title = element_text(size=5),
#             legend.key.size = unit(3, "mm"),
#             axis.text = element_text(size=6),
#             axis.title = element_text(size=6),
#             axis.line = element_line(linewidth=.2),
#             axis.ticks = element_line(linewidth = .2))+ 
#       scale_size_area(max_size = 3)
#     
#     plots <- append(plots, list(p1))
#   }
#   
#   if(!is.null(mp_enrichment_results)){
#     p2 <- dotplot(mp_enrichment_results, font.size=5)+
#       theme_classic()+
#       theme(axis.text.x = element_text(size = 5),
#             legend.text = element_text(size=3),
#             legend.title = element_text(size=5),
#             legend.key.size = unit(3, "mm"),
#             axis.text = element_text(size=6),
#             axis.title = element_text(size=6),
#             axis.line = element_line(linewidth=.2),
#             axis.ticks = element_line(linewidth = .2)) + 
#       scale_size_area(max_size = 3)
#     
#     plots <- append(plots, list(p2))
#   }
#   
#   if(!is.null(go_enrichment_results)){
#     p3 <- dotplot(go_enrichment_results, font.size=5)+
#       theme_classic()+
#       theme(axis.text.x = element_text(size = 5),
#             legend.text = element_text(size=3),
#             legend.title = element_text(size=5),
#             legend.key.size = unit(3, "mm"),
#             axis.text = element_text(size=6),
#             axis.title = element_text(size=6),
#             axis.line = element_line(linewidth=.2),
#             axis.ticks = element_line(linewidth = .2))+ 
#       scale_size_area(max_size = 3)
#     
#     plots <- append(plots, list(p3))
#   }
#   
#   
#   if(length(plots) > 0){
#     figure <- ggarrange(plotlist = plots, ncol=length(plots))
#     
#     figure <- annotate_figure(figure, top = text_grob(paste0("Supercluster ", i), face = "bold", size = 8))
#     figure
#     
#     all_plots <- append(all_plots, list(figure))
#   }
# }
# 
# saveRDS(overlap_signatures, paste0(dataDirectory, "genesets/ecoli_sc_overlap.rds"))
# 
# p <- ggarrange(plotlist = all_plots, nrow=length(all_plots))
# 
# p
# 
# jpeg(paste0(plotDirectory,"figure_4d.jpg"), width=250, height = 90, units = "mm", res = 1000)
# print(p)
# dev.off()
