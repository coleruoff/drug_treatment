args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(tidyverse)
library(ggpubr)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
source("final_scripts/drug_treatment_functions.R")
set.seed(42)

# source("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_scripts/drug_treatment_functions.R")

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/"

# Hallmarks term2gene list
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, human_gene_symbol)

new_geneset_names <- sapply(m_t2g$gs_name, FUN = function(x) gsub("HALLMARK_", "", x))
new_geneset_names <- sapply(new_geneset_names, FUN = function(x) gsub("_", " ", x))
m_t2g$gs_name <- new_geneset_names

# ITH Meta-programs term2gene list
mp_t2g <- readRDS(paste0(dataDirectory, "genesets/ith_meta_programs_t2g.rds"))

#Remove 'specific' MPs
specifc_mps <- c("MP39 Metal-response","MP31 Alveolar","MP29 NPC/OPC","MP28 Oligo normal","MP27 Oligo Progenitor","MP38 Glutathione","MP41 Unassigned","MP35 Hemato-related-I","MP37 Hemato-related-II","MP32 Skin-pigmentation","MP36 IG","MP16 MES (glioma)","MP15 EMT IV")
mp_t2g <- mp_t2g %>%
  filter(!gs_name %in% specifc_mps)

# Read in gene universes for enrichment analyses
gene_universe_intersection <- readRDS(paste0(dataDirectory, "cell_line_gene_universe_intersection.rds"))
cell_line_universes <- readRDS(paste0(dataDirectory, "cell_line_universes.rds"))

################################################################################

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

supercluster_components <- readRDS(paste0(dataDirectory, "processed_data/supercluster_components.rds"))

# Create signatures for genes that are distinct to each cell line

total_sc_distinct <- list()
for(curr_sc in 1:2){
  
  a <- supercluster_components[[curr_sc]][1]
  b <- supercluster_components[[curr_sc]][2]
  c <- supercluster_components[[curr_sc]][3]
  
  curr_supercluster_signatures <- c(list(all_signatures[["A549"]][[a]]),list(all_signatures[["K562"]][[b]]),list(all_signatures[["MCF7"]][[c]]))
  supercluster_2_signature <- list("supercluster_signature" = find_consensus_geneset(curr_supercluster_signatures,2))
  
  sc_distinct <- list()
  sc_distinct[["A549"]] <- curr_supercluster_signatures[[1]][!curr_supercluster_signatures[[1]] %in% supercluster_2_signature$supercluster_signature]
  sc_distinct[["K562"]] <- curr_supercluster_signatures[[2]][!curr_supercluster_signatures[[2]] %in% supercluster_2_signature$supercluster_signature]
  sc_distinct[["MCF7"]] <- curr_supercluster_signatures[[3]][!curr_supercluster_signatures[[3]] %in% supercluster_2_signature$supercluster_signature]
  
  total_sc_distinct[[curr_sc]] <- sc_distinct
  
}


# Create enrichment object for each supercluster
#     (hallmarks, ITH MPs, and GO BPs)
all_results <- list()

for(curr_sc in 1:2){
  cat("Supercluster ", curr_sc,"\n")
  
  curr_results <- list()
  for(curr_cell_line in cell_lines){
    
    universe_to_use <- cell_line_universes[[curr_cell_line]]
    
    curr_geneset <- total_sc_distinct[[curr_sc]][[curr_cell_line]]
    
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
    
    
    
    curr_results[[paste0(curr_cell_line, "_", supercluster_components[[curr_sc]][[curr_cell_line]])]] <- list(hallmark_enrichment_results,
                                                                                                              mp_enrichment_results,
                                                                                                              go_enrichment_results)
  }
  
  all_results[[curr_sc]] <- curr_results
}



# Create heatmaps for each supercluster
titles <- list("Cancer Hallmarks","ITH Meta-programs","GO Pathways")
all_plots <- list()
all_min_max <- list()

curr_sc <- 1

for(curr_sc in 1:2){ #Superclusters
  
  
  i <- 2
  curr_sc_plots <- list()
  for(i in 1:3){ #Enrichment Type
    curr_results <- list()
    
    for(j in 1:3){ #Supercluster Components
      curr_results <- append(curr_results, all_results[[curr_sc]][[j]][[i]])
    }
    
    names(curr_results) <- cell_lines#names(all_results[[curr_sc]])
    
    
    
    nrow(curr_results$MCF7)
    
    fun_results <- create_enrichment_heatmap(curr_results, titles[i])
    
    heatmap <- fun_results[[1]]
    
    all_min_max <- append(all_min_max, list(fun_results[[2]]))
    
    heatmap@row_names_param$gp$fontsize <- 5
    heatmap@column_names_param$gp$fontsize <- 8
    heatmap@column_title_param$gp$fontsize <- 10
    heatmap@matrix_legend_param$title_gp$fontsize <- 10
    heatmap@matrix_legend_param$labels_gp$fontsize <- 10
    
    
    # Add heatmap to curr cell line list of heatmaps
    curr_sc_plots <- append(curr_sc_plots, heatmap)
    
  }
  
  all_plots[[curr_sc]] <- curr_sc_plots
}

total_min <- 0
total_max <- 0
for(i in length(all_min_max)){
  
  if(all_min_max[[i]][1] < total_min){
    total_min <- all_min_max[[i]][1] 
  }
  
  if(all_min_max[[i]][2] > total_max){
    total_max <- all_min_max[[i]][2] 
  }
  
}

col_fun <- colorRamp2(c(0,total_max), c("white", "red1"))

# png(paste0(plotDirectory, "figure_S2c.png"),
#     width=36, height=30, units= "in", res = 300)


# grid.newpage()
# 
# top.vp <- viewport(layout=grid.layout(4, 3,
#                                       widths=unit(c(1, 1, 1), c("null", "null", "null")),
#                                       heights=unit(c(.5,5,.5,5), c("null", "null", "null"))))
# 
# title1 <- viewport(layout.pos.col = 2, layout.pos.row = 1, name = "title1")
# plot1 <- viewport(layout.pos.col = 1, layout.pos.row = 2, name = "plot1")
# plot2 <- viewport(layout.pos.col = 2, layout.pos.row = 2, name = "plot2")
# plot3 <- viewport(layout.pos.col = 3, layout.pos.row = 2, name = "plot3")
# title2 <- viewport(layout.pos.col = 2, layout.pos.row = 3, name = "title2")
# plot4 <- viewport(layout.pos.col = 1, layout.pos.row = 4, name = "plot4")
# plot5 <- viewport(layout.pos.col = 2, layout.pos.row = 4, name = "plot5")
# plot6 <- viewport(layout.pos.col = 3, layout.pos.row = 4, name = "plot6")
# 
# splot <- vpTree(top.vp, vpList(title1, plot1, plot2, plot3, title2, plot4, plot5, plot6))
# 
# pushViewport(splot)

# Plot heatmaps for each supercluster

plot_rows <- list()
curr_sc <- 1
for(curr_sc in 1:2){
  
  all_plots[[curr_sc]][[1]]@matrix_color_mapping@col_fun <- col_fun

  all_plots[[curr_sc]][[2]]@matrix_color_mapping@col_fun <- col_fun

  all_plots[[curr_sc]][[3]]@matrix_color_mapping@col_fun <- col_fun
  
  curr_title <- paste0("Supercluster ", curr_sc)
  
  
  all_plots[[curr_sc]][[1]]@heatmap_param$show_heatmap_legend <- T
  all_plots[[curr_sc]][[2]]@heatmap_param$show_heatmap_legend <- T
  all_plots[[curr_sc]][[3]]@heatmap_param$show_heatmap_legend <- T
  
  ht_grob1 = grid.grabExpr(draw(all_plots[[curr_sc]][[1]]))
  ht_grob2 = grid.grabExpr(draw(all_plots[[curr_sc]][[2]]))
  ht_grob3 = grid.grabExpr(draw(all_plots[[curr_sc]][[3]]))
  
  
  p <- plot_grid(ht_grob1,ht_grob2, ht_grob3, nrow=1, rel_widths = c(2, 2, 2))
  
  plot_rows <- append(plot_rows, list(p))
  
  # if(curr_sc == 1){
  #   seekViewport("plot1")
  #   grid.draw(ht_grob1)
  #   
  #   seekViewport("plot2")
  #   grid.draw(ht_grob2)
  #   
  #   seekViewport("plot3")
  #   grid.draw(ht_grob3)
  #   
  #   seekViewport("title1")
  #   grid.text(curr_title, gp = gpar(fontsize = 40, fontface = "bold"))
  # } else {
  #   seekViewport("plot4")
  #   grid.draw(ht_grob1)
  #   
  #   seekViewport("plot5")
  #   grid.draw(ht_grob2)
  #   
  #   seekViewport("plot6")
  #   grid.draw(ht_grob3)
  #   
  #   seekViewport("title2")
  #   grid.text(curr_title, gp = gpar(fontsize = 40, fontface = "bold"))
  # }
}



p <- plot_grid(plot_rows[[1]],plot_rows[[2]],nrow=2)

tiff(paste0(plotDirectory,"figure_S2c.tiff"), width=190, height = 200, units = "mm", res = 1000)

print(p)

dev.off()







