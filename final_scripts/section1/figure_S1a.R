args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")

set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_figures/"

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
#################################################################################

resistance_signature <- readRDS(paste0(dataDirectory, "genesets/resistance_signature.rds"))

names(resistance_signature) <- "Resistance Signature"

final_plot <- plot_enrichment_ora(resistance_signature)

jpeg(paste0(plotDirectory,"figure_S1a.jpg"), width=250, height=80, units = "mm", res = 1000)
print(final_plot)
dev.off()




# 
# 
# curr_geneset <- resistance_signature$resistance_signature
# 
# hallmark_enrichment_results <- enricher(curr_geneset, TERM2GENE=m_t2g)
# 
# # MPs Enrichment
# mp_enrichment_results <- enricher(curr_geneset, TERM2GENE=mp_t2g)
# 
# 
# # GO Functional Enrichment
# go_enrichment_results <- enrichGO(gene          = curr_geneset,
#                                   OrgDb         = org.Hs.eg.db,
#                                   ont           = "ALL",
#                                   pAdjustMethod = "BH",
#                                   pvalueCutoff  = .01,
#                                   qvalueCutoff  = .05,
#                                   readable      = TRUE,
#                                   keyType = "SYMBOL")
# 
# 
# # hallmark_enrichment_results@result$Description[1] <- "OXIDATIVE\nPHOSPHORYLATION"
# 
# p1 <- dotplot(hallmark_enrichment_results, font.size=3)+
#   theme_classic()+
#   theme(axis.text.x = element_text(size = 5),
#         legend.text = element_text(size=3),
#         legend.title = element_text(size=5),
#         legend.key.size = unit(3, "mm"),
#         axis.text = element_text(size=6),
#         axis.title = element_text(size=6),
#         axis.line = element_line(linewidth=.2),
#         axis.ticks = element_line(linewidth = .2))+ 
#   scale_size_area(max_size = 3)
# 
# p2 <- dotplot(mp_enrichment_results, font.size=3)+
#   theme_classic()+
#   theme(axis.text.x = element_text(size = 5),
#         legend.text = element_text(size=3),
#         legend.title = element_text(size=5),
#         legend.key.size = unit(3, "mm"),
#         axis.text = element_text(size=6),
#         axis.title = element_text(size=6),
#         axis.line = element_line(linewidth=.2),
#         axis.ticks = element_line(linewidth = .2))+ 
#   scale_size_area(max_size = 3)
# 
# p3 <- dotplot(go_enrichment_results, font.size=3)+
#   theme_classic()+
#   theme(axis.text.x = element_text(size = 5),
#         legend.text = element_text(size=3),
#         legend.title = element_text(size=5),
#         legend.key.size = unit(3, "mm"),
#         axis.text = element_text(size=6),
#         axis.title = element_text(size=6),
#         axis.line = element_line(linewidth=.2),
#         axis.ticks = element_line(linewidth = .2))+ 
#   scale_size_area(max_size = 3)
# 
# 
# 
# p <- p1+p2+p3
# p
# 
# 
# jpeg(paste0(plotDirectory,"figure_S1a.jpg"), width=250, height=80, units = "mm", res = 1000)
# print(p)
# dev.off()


