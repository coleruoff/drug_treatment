args = commandArgs(trailingOnly=TRUE)
# dataDirectory <- paste0(args[1],"final_data/")
# plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("revision_scripts/drug_treatment_functions.R")

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_figures/"

################################################################################


supercluster_signatures <- readRDS(paste0(dataDirectory, "genesets/supercluster_down_signatures.rds"))

names(supercluster_signatures) <- c("Supercluster 1","Supercluster 2","Supercluster 3")

final_plot <- plot_enrichment_ora(supercluster_signatures)

jpeg(paste0(plotDirectory,"figure_S2b.jpg"), width=250, height = 200, units = "mm", res = 1000)
print(final_plot)
dev.off()




# 
# all_results <- list()
# 
# 
# for(curr_sc in 1:length(supercluster_signatures)){
#   
#   curr_geneset <- supercluster_signatures[[curr_sc]]
#   
#   # Hallmarks Enrichment
#   hallmark_enrichment_results <- enricher(curr_geneset, TERM2GENE=m_t2g)
#   # universe = universe_to_use)
#   
#   # MPs Enrichment
#   mp_enrichment_results <- enricher(curr_geneset, TERM2GENE=mp_t2g)
#   # universe = universe_to_use)
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
#   # universe = universe_to_use)
#   
#   
#   all_results[[paste0("sc",curr_sc)]] <- list(hallmark_enrichment_results,
#                                               mp_enrichment_results,
#                                               go_enrichment_results)
#   
# }
# 
# titles <- list("Cancer Hallmarks","ITH Meta-programs","GO Pathways")
# all_plots <- list()
# 
# for(i in 1:3){
#   
#   
#   plots <- list()
#   for(j in 1:3){
#     df <- all_results[[j]][[i]]
#     
#     if(nrow(df) > 0){
#       p <- dotplot(df, font.size=3, showCategory = 5)+
#         ggtitle(paste0("Supercluster ", j))+
#         theme_classic()+
#         theme(plot.title = element_text(size=6),
#               axis.text.x = element_text(size = 5),
#               legend.text = element_text(size=3),
#               legend.title = element_text(size=5),
#               legend.key.size = unit(3, "mm"),
#               axis.text = element_text(size=6),
#               axis.title = element_text(size=6),
#               axis.line = element_line(linewidth=.2),
#               axis.ticks = element_line(linewidth = .2))+ 
#         scale_size_area(max_size = 3)
#       
#       plots <- append(plots, list(p))
#     }
#     
#   }
#   
#   plot <- ggarrange(plotlist = plots, nrow=3, common.legend = T, legend = "right")
#   
#   plot <- annotate_figure(plot, top = text_grob(titles[i], size = 8))
#   
#   all_plots <- append(all_plots, list(plot))
# }
# 
# final_plot <- ggarrange(plotlist = all_plots, ncol=3,common.legend = T, legend = "right")
# 
# jpeg(paste0(plotDirectory,"figure_S2b.jpg"), width=250, height = 200, units = "mm", res = 1000)
# print(final_plot)
# dev.off()

################################################################################

