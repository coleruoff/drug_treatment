args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
source("final_scripts/drug_treatment_functions.R")
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
# plotDirectory <- "/data/ruoffcj/projects/drug_treatment/final_figures/"

################################################################################

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, human_gene_symbol)

mp_t2g <- readRDS(paste0(dataDirectory,"genesets/ith_meta_programs_t2g.rds"))

specifc_mps <- c("MP39 Metal-response","MP31 Alveolar","MP29 NPC/OPC","MP28 Oligo normal","MP27 Oligo Progenitor","MP38 Glutathione","MP41 Unassigned","MP35 Hemato-related-I","MP37 Hemato-related-II","MP32 Skin-pigmentation","MP36 IG","MP16 MES (glioma)","MP15 EMT IV")

mp_t2g <- mp_t2g %>% 
  filter(!gs_name %in% specifc_mps)

gene_universe_intersection <- readRDS(paste0(dataDirectory, "cell_line_gene_universe_intersection.rds"))

ecoli_human_orthologs_table <- readRDS(paste0(dataDirectory, "ortholog_mapping/ecoli_human_orthologs_table.rds"))

gene_universe_intersection <- gene_universe_intersection[gene_universe_intersection %in% ecoli_human_orthologs_table$HUMAN_SYMBOL]

################################################################################
# Plotting for RAC superclusters signatures

ecoli_human_orthologs_up <- readRDS(paste0(dataDirectory, "genesets/ecoli_human_orthologs_up.rds"))
supercluster_signatures <- readRDS(paste0(dataDirectory, "genesets/rac_supercluster_signatures.rds"))

curr_geneset <- intersect(supercluster_signatures[[1]],ecoli_human_orthologs_up$ecoli_human_orthologs_up)

universe_to_use <- gene_universe_intersection
all_results <- list()

# Hallmarks Enrichment
hallmark_enrichment_results <- enricher(curr_geneset, TERM2GENE=m_t2g,
                                        universe = universe_to_use)

# MPs Enrichment
mp_enrichment_results <- enricher(curr_geneset, TERM2GENE=mp_t2g,
                                  universe = universe_to_use)

# GO Functional Enrichment
go_enrichment_results <- enrichGO(gene          = curr_geneset,
                                  OrgDb         = org.Hs.eg.db,
                                  ont           = "BP",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = .01,
                                  qvalueCutoff  = .05,
                                  readable      = TRUE,
                                  keyType = "SYMBOL",
                                  universe = universe_to_use)


p1 <- dotplot(hallmark_enrichment_results)
p2 <- dotplot(mp_enrichment_results)
p3 <- dotplot(go_enrichment_results)

p3+theme(axis.text.y = element_text(size=8))

plots <- list(p3)


figure <- ggarrange(plotlist = plots, ncol=1, common.legend = T,legend=c("right"))

p <- annotate_figure(figure,
                     top=text_grob("Shared Genes Between E. coli Resistance Orthologs and Supercluster 1 Signature Enrichment", size=14, face="bold"))


png(paste0(plotDirectory, "figure_5d.png"),
    width = 10,height=10, units = 'in',res = 300)

print(p)

dev.off()


# all_results[[paste0("sc",1)]] <- list(hallmark_enrichment_results,
#                                             mp_enrichment_results,
#                                             go_enrichment_results)
#   
# 
# 
# titles <- list("Cancer Hallmarks","ITH Meta-programs","GO Pathways")
# all_plots <- list()
# 
# for(i in c(3)){
#   
#   # For each type of enrichment, create heatmap
#   curr_results <- list(all_results[["sc1"]][[i]])
#   names(curr_results) <- c("Supercluster 1")
#   
#   # curr_results <- curr_results[!sapply(curr_results, is.null)]
#   heatmap <- create_enrichment_heatmap(curr_results, titles[i])
#   
#   # Add heatmap to curr cell line list of heatmaps
#   all_plots <- append(all_plots, heatmap)
# }
# 
# 
# 
# curr_title <- paste0("Shared Genes Between E. Coli Resistance Orthologs and Supercluster 1 Signature Enrichment")
# 
# ht_grob1 = grid.grabExpr(draw(all_plots[[1]], padding = unit(c(0, 80, 0, 0), "mm")))
# # ht_grob2 = grid.grabExpr(draw(all_plots[[2]], padding = unit(c(0, 80, 0, 0), "mm")))
# # ht_grob3 = grid.grabExpr(draw(all_plots[[3]], padding = unit(c(0, 80, 0, 0), "mm")))
# 
# 
# 
# png(paste0(plotDirectory, "final_figures/figure_5d.png"),
#     width=36, height=15, units= "in", res = 300)
# 
# grid.newpage()
# 
# top.vp <- viewport(layout=grid.layout(2, 3,
#                                       widths=unit(c(.2,2,.2), c("null", "null", "null")),
#                                       heights=unit(c(.5,5), c("null", "null", "null"))))
# 
# title1 <- viewport(layout.pos.col = 2, layout.pos.row = 1, name = "title1")
# plot1 <- viewport(layout.pos.col = 2, layout.pos.row = 2, name = "plot1")
# # plot2 <- viewport(layout.pos.col = 3, layout.pos.row = 2, name = "plot2")
# # plot3 <- viewport(layout.pos.col = 3, layout.pos.row = 2, name = "plot3")
# 
# splot <- vpTree(top.vp, vpList(title1, plot1, plot2, plot3))
# 
# pushViewport(splot)
# 
# seekViewport("plot1")
# grid.draw(ht_grob1)
# 
# # seekViewport("plot2")
# # grid.draw(ht_grob2)
# 
# # seekViewport("plot3")
# # grid.draw(ht_grob3)
# 
# seekViewport("title1")
# grid.text(curr_title, gp = gpar(fontsize = 40))
# 
# dev.off()
