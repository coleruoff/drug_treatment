args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_figures/"
################################################################################
# Plotting for supercluster enrichment

supercluster_signatures <- readRDS(paste0(dataDirectory, "genesets/supercluster_up_signatures.rds"))

names(supercluster_signatures) <- c("Supercluster 1","Supercluster 2","Supercluster 3")

final_plot <- plot_enrichment_ora(supercluster_signatures)


jpeg(paste0(plotDirectory,"figure_2b.jpg"), width=300, height = 180, units = "mm", res = 1000)
print(final_plot)
dev.off()

# supercluster_signatures <- readRDS(paste0(dataDirectory, "genesets/supercluster_down_signatures.rds"))

# names(supercluster_signatures) <- c("Supercluster 1","Supercluster 2","Supercluster 3","Supercluster 4")
# names(supercluster_signatures) <- c("Supercluster 1","Supercluster 2","Supercluster 3")

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
#     p <- barplot(df, font.size=3, showCategory = 5)+
#       ggtitle(paste0("Supercluster ", j))+
#       theme_classic()+
#       scale_fill_gradient(low="firebrick4", high="royalblue4")+
#       theme(plot.title = element_text(size=10),
#             axis.text.x = element_text(size = 8),
#             legend.text = element_text(size=6),
#             legend.title = element_text(size=8),
#             legend.key.size = unit(4, "mm"),
#             axis.text = element_text(size=8),
#             axis.title = element_text(size=8),
#             axis.line = element_line(linewidth=.2),
#             axis.ticks = element_line(linewidth = .2))+ 
#       scale_size_area(max_size = 3)
#       
#     plots <- append(plots, list(p))
#   }
#   
#   plot <- ggarrange(plotlist = plots, nrow=3, common.legend = T, legend = "right")
#   
#   plot <- annotate_figure(plot, top = text_grob(titles[i], size = 12))
#   
#   all_plots <- append(all_plots, list(plot))
# }
# 
# final_plot <- ggarrange(plotlist = all_plots, ncol=3, common.legend = T, legend = "right")

# jpeg(paste0(plotDirectory,"figure_2b.jpg"), width=300, height = 180, units = "mm", res = 1000)
# print(final_plot)
# dev.off()

################################################################################






# 
# all_plots <- list()
# 
# all_min_max <- list()
# for(i in 1:3){
#   
#   # For each type of enrichment, create heatmap
#   curr_results <- list(all_results[["sc1"]][[i]],all_results[["sc2"]][[i]],all_results[["sc3"]][[i]])
#   # curr_results <- list(all_results[["sc1"]][[i]],all_results[["sc2"]][[i]],all_results[["sc3"]][[i]],all_results[["sc4"]][[i]])
#   names(curr_results) <- c("Supercluster 1","Supercluster 2","Supercluster 3")
#   
#   # fun_results <- create_enrichment_heatmap(curr_results, titles[i])
#   fun_results <- create_enrichment_dotplot(curr_results, titles[i])
#   
#   heatmap <- fun_results[[1]]
#   
#   all_min_max <- append(all_min_max, list(fun_results[[2]]))
#   
#   heatmap@row_names_param$gp$fontsize <- 4
#   heatmap@column_names_param$gp$fontsize <- 4
#   heatmap@column_title_param$gp$fontsize <- 6
#   heatmap@matrix_legend_param$title_gp$fontsize <- 6
#   heatmap@matrix_legend_param$labels_gp$fontsize <- 6
#   
#   
#   # Add heatmap to curr cell line list of heatmaps
#   all_plots <- append(all_plots, heatmap)
# }
# 
# total_min <- 0
# total_max <- 0
# for(i in length(all_min_max)){
#   
#   if(all_min_max[[i]][1] < total_min){
#     total_min <- all_min_max[[i]][1] 
#   }
#   
#   if(all_min_max[[i]][2] > total_max){
#     total_max <- all_min_max[[i]][2] 
#   }
#   
# }
# 
# 
# col_fun = colorRamp2(c(0,total_max), c("white", "red1"))
# 
# all_plots[[1]]@matrix_color_mapping@col_fun <- col_fun
# all_plots[[2]]@matrix_color_mapping@col_fun <- col_fun
# all_plots[[3]]@matrix_color_mapping@col_fun <- col_fun
# 
# 
# all_plots[[3]]@heatmap_param$show_heatmap_legend <- T
# 
# ht_grob1 = grid.grabExpr(draw(all_plots[[1]]))#, column_title = "Cancer Hallmarks", column_title_gp = gpar(fontsize=6)))
# ht_grob2 = grid.grabExpr(draw(all_plots[[2]]))#, column_title = "ITH Meta-programs", column_title_gp = gpar(fontsize=6)))
# ht_grob3 = grid.grabExpr(draw(all_plots[[3]]))#, column_title = "GO Pathways", column_title_gp = gpar(fontsize=6)))
# 
# 
# p <- plot_grid(ht_grob1,ht_grob2, ht_grob3, nrow=1, rel_widths = c(.9,.9,1))
# 
# 
# p

jpeg(paste0(plotDirectory,"figure_2b.jpg"), width=200, height = 100, units = "mm", res = 1000)
print(p)
dev.off()

