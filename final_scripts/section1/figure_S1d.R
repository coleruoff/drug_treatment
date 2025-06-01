args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")

set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/"

################################################################################
# Plotting for Global RAC enrichment

all_signatures <- readRDS(paste0(dataDirectory, "genesets/global_rac_signatures.rds"))

final_plot <- plot_enrichment_ora(all_signatures)

jpeg(paste0(plotDirectory,"figure_S1d.jpg"), width=350, height = 200, units = "mm", res = 600)
print(final_plot)
dev.off()


# 
# 
# all_results <- list()
# 
# for(curr_cell_line in cell_lines){
#   
#   # Get current cell line RAC ranks
#   ranks <- all_ranks[[curr_cell_line]]
#   
#   ###########
#   # Hallmarks
#   if(curr_cell_line == "K562"){
#     y_label <- "Cancer Hallmarks"
#   }
#   
#   hallmarks_gsea <- GSEA(ranks, TERM2GENE = m_t2g)
#   
#   # hallmark_plot <- create_barplot(hallmarks_gsea, y_label, curr_cell_line)
#   
#   ###############
#   # Meta-programs
#   if(curr_cell_line == "K562"){
#     y_label <- "ITH Meta-Programs"
#   }
#   
#   mp_gsea <- GSEA(ranks, TERM2GENE = mp_t2g)
#   
#   # mp_plot <- create_barplot(mp_gsea, y_label, curr_cell_line)
#   
#   #############
#   # GO Pathways
#   if(curr_cell_line == "K562"){
#     y_label <- "GO Pathways"
#   }
#   
#   go_gsea <- gseGO(geneList     = ranks,
#                    OrgDb        = org.Hs.eg.db,
#                    ont          = "BP",
#                    minGSSize    = 100,
#                    maxGSSize    = 500,
#                    pvalueCutoff = 0.01,
#                    verbose      = FALSE,
#                    keyType = "SYMBOL")
#   
#   
#   
#   all_results[[curr_cell_line]] <- list(hallmarks_gsea,mp_gsea,go_gsea)
#   
# }
# 
# 
# titles <- list("Cancer Hallmarks","ITH Meta-programs","GO Pathways")
# all_plots <- list()
# 
# cell_lines
# 
# for(i in 1:3){
#   
#   
#   plots <- list()
#   for(j in 1:3){
#     df <- all_results[[j]][[i]]
#     
#     p <- dotplot(df, font.size=3, showCategory = 5)+
#       ggtitle(paste0(cell_lines[j]))+
#       theme_classic()+
#       theme(plot.title = element_text(size=6),
#             axis.text.x = element_text(size = 5),
#             legend.text = element_text(size=3),
#             legend.title = element_text(size=5),
#             legend.key.size = unit(3, "mm"),
#             axis.text = element_text(size=6),
#             axis.title = element_text(size=6),
#             axis.line = element_line(linewidth=.2),
#             axis.ticks = element_line(linewidth = .2))+ 
#       scale_size_area(max_size = 3)
#     
#     plots <- append(plots, list(p))
#   }
#   
#   plot <- ggarrange(plotlist = plots, nrow=3, common.legend = T, legend = "right")
#   
#   plot <- annotate_figure(plot, top = text_grob(titles[i], size = 8))
#   
#   all_plots <- append(all_plots, list(plot))
# }
# 
# final_plot <- ggarrange(plotlist = all_plots, ncol=3, common.legend = T, legend = "right")
# 
# jpeg(paste0(plotDirectory,"figure_S1d.jpg"), width=250, height = 200, units = "mm", res = 1000)
# print(final_plot)
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 










# all_min_max <- list()
# i <- 1
# for(i in 1:3){
#   
#   # For each type of enrichment, create heatmap
#   curr_results <- list(all_results[["A549"]][[i]],all_results[["K562"]][[i]],all_results[["MCF7"]][[i]])
#   names(curr_results) <- cell_lines
#   
#   fun_results <- create_enrichment_heatmap(curr_results, titles[i])
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
#   # Add heatmap to curr cell line list of heatmaps
#   all_plots <- append(all_plots, heatmap)
# }
# 
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
# col_fun = colorRamp2(c(0,total_max), c("white", "red1"))
# 
# all_plots[[1]]@matrix_color_mapping@col_fun <- col_fun
# all_plots[[2]]@matrix_color_mapping@col_fun <- col_fun
# all_plots[[3]]@matrix_color_mapping@col_fun <- col_fun
# 
# 
# all_plots[[3]]@heatmap_param$show_heatmap_legend <- T
# 
# ht_grob1 = grid.grabExpr(draw(all_plots[[1]]))
# ht_grob2 = grid.grabExpr(draw(all_plots[[2]]))
# ht_grob3 = grid.grabExpr(draw(all_plots[[3]]))
# 
# 
# 
# p <- plot_grid(ht_grob1,ht_grob2, ht_grob3, nrow=1, rel_widths = c(2, 2, 2))
# 
# p
# 
# 
# jpeg(paste0(plotDirectory,"figure_S1d.jpg"), width=200, height = 100, units = "mm", res = 1000)
# print(p)
# dev.off()



