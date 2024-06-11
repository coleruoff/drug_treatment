args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
library(tidyverse)
source("final_scripts/drug_treatment_functions.R")
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
# plotDirectory <- "/data/ruoffcj/projects/drug_treatment/final_figures/"

################################################################################

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, human_gene_symbol)

mp_t2g <- readRDS(paste0(dataDirectory, "genesets/ith_meta_programs_t2g.rds"))

specifc_mps <- c("MP39 Metal-response","MP31 Alveolar","MP29 NPC/OPC","MP28 Oligo normal","MP27 Oligo Progenitor","MP38 Glutathione","MP41 Unassigned","MP35 Hemato-related-I","MP37 Hemato-related-II","MP32 Skin-pigmentation","MP36 IG","MP16 MES (glioma)","MP15 EMT IV")

mp_t2g <- mp_t2g %>% 
  filter(!gs_name %in% specifc_mps)

gene_universe_intersection <- readRDS(paste0(dataDirectory, "cell_line_gene_universe_intersection.rds"))

#Filter universe to only genes that have a yeast ortholog

#Read in conversion table from OMA browser
human_yeast_conversion <- read.table(paste0(dataDirectory, "raw_data/yeast_data/human_candida_auris_source_id.txt"), fill = T, sep='\t')[,1:3]
colnames(human_yeast_conversion) <- c("HUMAN","YEAST", "MAPPING")

# Trim ensembl gene names
human_yeast_conversion$HUMAN <- gsub("\\..*", "", human_yeast_conversion$HUMAN)


# Add HUMAN ensembl to symbol conversion
human_ensembl_symbol_conversion <- AnnotationDbi::select(org.Hs.eg.db, 
                                                         keys = gene_universe_intersection,
                                                         columns = c("ENSEMBL", "SYMBOL"),
                                                         keytype = "SYMBOL")


gene_universe_intersection <- human_ensembl_symbol_conversion$SYMBOL[human_ensembl_symbol_conversion$ENSEMBL %in% human_yeast_conversion$HUMAN]

################################################################################
# Creating supercluster signatures that consist of all genes that appear in any
# of the 3 componenet clusters upregulated genes

cell_lines <- c("A549","K562","MCF7")

num_clusters <- c(19,12,18)
names(num_clusters) <- cell_lines

all_signatures <- list()
for(curr_cell_line in cell_lines){
  
  de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_de.rds"))
  
  curr_signatures <- list()
  for(curr_cluster in 1:num_clusters[curr_cell_line]){
    
    # de_res <- readRDS(paste0(dataDirectory, "de_results/", curr_cell_line, "_cluster_", curr_cluster, "_de.rds"))
    
    curr_cluster_signature <- de_res %>% 
      filter(cluster == curr_cluster & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
      arrange(desc(avg_log2FC)) %>% 
      pull(gene) 
    
    
    curr_signatures[[curr_cluster]] <- curr_cluster_signature
    
  }
  
  all_signatures[[curr_cell_line]] <- curr_signatures
}
################################################################################

supercluster_components <- readRDS(paste0(dataDirectory, "processed_data/supercluster_components.rds"))

supercluster1_signatures <- c(list(all_signatures[["A549"]][[supercluster_components[[1]][["A549"]]]]),
                              list(all_signatures[["K562"]][[supercluster_components[[1]][["K562"]]]]),
                              list(all_signatures[["MCF7"]][[supercluster_components[[1]][["MCF7"]]]]))

supercluster1_signature <- list("supercluster1_signature" = find_consensus_geneset(supercluster1_signatures,1))

##############

supercluster2_signatures <- c(list(all_signatures[["A549"]][[supercluster_components[[2]][["A549"]]]]),
                              list(all_signatures[["K562"]][[supercluster_components[[2]][["K562"]]]]),
                              list(all_signatures[["MCF7"]][[supercluster_components[[2]][["MCF7"]]]]))

supercluster2_signature <- list("supercluster2_signature" = find_consensus_geneset(supercluster2_signatures,1))

supercluster_signatures <- c(supercluster1_signature,supercluster2_signature)

################################################################################
# Plotting for RAC superclusters signatures

yeast_human_orthologs_up <- readRDS(paste0(dataDirectory, "genesets/yeast_human_orthologs_up.rds"))
# supercluster_signatures <- readRDS(paste0(dataDirectory, "genesets/rac_supercluster_signatures.rds"))

names2 <- gsub("_", " ", names(supercluster_signatures))
names3 <- gsub("supercluster", "Supercluster ", names2)
names(supercluster_signatures) <- gsub("signature", "Signature ", names3)

curr_geneset <- intersect(supercluster_signatures[[2]], yeast_human_orthologs_up$yeast_human_orthologs_up)
saveRDS(curr_geneset, paste0(dataDirectory, "genesets/yeast_sc2_overlap.rds"))

universe_to_use <- gene_universe_intersection

names(supercluster_signatures) <- c("Supercluster 1","Supercluster 2")

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

# plots <- list(p1,p2,p3)
plots <- list(p1,p2,p3)


figure <- ggarrange(plotlist = plots, ncol=3, common.legend = T,legend=c("right"))

p <- annotate_figure(figure,
                     top=text_grob("Shared Genes Between Yeast Resistance Orthologs and Supercluster 2 Signature Enrichment", size=28, face="bold"))


png(paste0(plotDirectory, "figure_5b.png"),
    width = 24,height=10, units = 'in',res = 300)

print(p)

dev.off()

# 
# all_results[[paste0("sc",1)]] <- list(hallmark_enrichment_results,
#                                             mp_enrichment_results,
#                                             go_enrichment_results)
#   
# 
# 
# titles <- list("Cancer Hallmarks","ITH Meta-programs","GO Pathways")
# all_plots <- list()
# 
# for(i in 1:3){
#   
#   # For each type of enrichment, create heatmap
#   curr_results <- list(all_results[["sc1"]][[i]])
#   names(curr_results) <- c("Supercluster 1")
#   
#   heatmap <- create_enrichment_heatmap(curr_results, titles[i])
#   
#   # Add heatmap to curr cell line list of heatmaps
#   all_plots <- append(all_plots, heatmap)
# }
# 
# 
# 
# curr_title <- paste0("Shared Genes Between Yeast Resistance Orthologs and Supercluster 1 Signature Enrichment")
# 
# ht_grob1 = grid.grabExpr(draw(all_plots[[1]], padding = unit(c(0, 80, 0, 0), "mm")))
# ht_grob2 = grid.grabExpr(draw(all_plots[[2]], padding = unit(c(0, 80, 0, 0), "mm")))
# ht_grob3 = grid.grabExpr(draw(all_plots[[3]], padding = unit(c(0, 80, 0, 0), "mm")))
# 
# 
# 
# png(paste0(plotDirectory, "figure_5b.png"),
#     width=36, height=15, units= "in", res = 300)
# 
# grid.newpage()
# 
# top.vp <- viewport(layout=grid.layout(2, 3,
#                                       widths=unit(c(1,1,1), c("null", "null", "null")),
#                                       heights=unit(c(.5,5), c("null", "null", "null"))))
# 
# title1 <- viewport(layout.pos.col = 2, layout.pos.row = 1, name = "title1")
# plot1 <- viewport(layout.pos.col = 1, layout.pos.row = 2, name = "plot1")
# plot2 <- viewport(layout.pos.col = 2, layout.pos.row = 2, name = "plot2")
# plot3 <- viewport(layout.pos.col = 3, layout.pos.row = 2, name = "plot3")
# 
# splot <- vpTree(top.vp, vpList(title1, plot1, plot2, plot3))
# 
# pushViewport(splot)
# 
# seekViewport("plot1")
# grid.draw(ht_grob1)
# 
# seekViewport("plot2")
# grid.draw(ht_grob2)
# 
# seekViewport("plot3")
# grid.draw(ht_grob3)
# 
# seekViewport("title1")
# grid.text(curr_title, gp = gpar(fontsize = 40))
# 
# dev.off()


