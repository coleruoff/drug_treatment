args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
source("final_scripts/drug_treatment_functions.R")
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
library(tidyverse)
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
# plotDirectory <- "/data/ruoffcj/projects/drug_treatment/final_figures/"

################################################################################

# Hallmarks term2gene list
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, human_gene_symbol)

# ITH Meta-programs term2gene list
mp_t2g <- readRDS(paste0(dataDirectory, "genesets/ith_meta_programs_t2g.rds"))

#Remove 'specific' MPs
specifc_mps <- c("MP39 Metal-response","MP31 Alveolar","MP29 NPC/OPC","MP28 Oligo normal","MP27 Oligo Progenitor","MP38 Glutathione","MP41 Unassigned","MP35 Hemato-related-I","MP37 Hemato-related-II","MP32 Skin-pigmentation","MP36 IG","MP16 MES (glioma)","MP15 EMT IV")
mp_t2g <- mp_t2g %>% 
  filter(!gs_name %in% specifc_mps)

cell_lines <- c("A549","K562","MCF7")

# Read in gene universes for enrichment analyses
gene_universe_intersection <- readRDS(paste0(dataDirectory, "cell_line_gene_universe_intersection.rds"))
cell_line_universes <- readRDS(paste0(dataDirectory, "cell_line_universes.rds"))

################################################################################
# Plotting for Global RAC enrichment

all_ranks <- readRDS(paste0(dataDirectory, "genesets/global_rac_ranks.rds"))

all_results <- list()

for(curr_cell_line in cell_lines){
  
  # Get current cell line RAC ranks
  ranks <- all_ranks[[curr_cell_line]]
  
  ###########
  # Hallmarks
  if(curr_cell_line == "K562"){
    y_label <- "Cancer Hallmarks"
  }
  
  hallmarks_gsea <- GSEA(ranks, TERM2GENE = m_t2g)
  
  # hallmark_plot <- create_barplot(hallmarks_gsea, y_label, curr_cell_line)
  
  ###############
  # Meta-programs
  if(curr_cell_line == "K562"){
    y_label <- "ITH Meta-Programs"
  }
  
  mp_gsea <- GSEA(ranks, TERM2GENE = mp_t2g)
  
  # mp_plot <- create_barplot(mp_gsea, y_label, curr_cell_line)
  
  #############
  # GO Pathways
  if(curr_cell_line == "K562"){
    y_label <- "GO Pathways"
  }
  
  go_gsea <- gseGO(geneList     = ranks,
                   OrgDb        = org.Hs.eg.db,
                   ont          = "BP",
                   minGSSize    = 100,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.01,
                   verbose      = FALSE,
                   keyType = "SYMBOL")
  
  
  
  all_results[[curr_cell_line]] <- list(hallmarks_gsea,mp_gsea,go_gsea)
  
}


titles <- list("Cancer Hallmarks","ITH Meta-programs","GO Pathways")
all_plots <- list()
i <- 1
for(i in 1:3){
  
  # For each type of enrichment, create heatmap
  curr_results <- list(all_results[["A549"]][[i]],all_results[["K562"]][[i]],all_results[["MCF7"]][[i]])
  names(curr_results) <- cell_lines
  
  heatmap <- create_enrichment_heatmap(curr_results, titles[i])
  
  # Add heatmap to curr cell line list of heatmaps
  all_plots <- append(all_plots, heatmap)
}




ht_grob1 = grid.grabExpr(draw(all_plots[[1]], padding = unit(c(0, 70, 0, 0), "mm")))
ht_grob2 = grid.grabExpr(draw(all_plots[[2]], padding = unit(c(0, 60, 0, 0), "mm")))
ht_grob3 = grid.grabExpr(draw(all_plots[[3]], padding = unit(c(0, 70, 0, 0), "mm")))


png(paste0(plotDirectory, "figure_S2b.png"),
    width=30, height=10, units= "in", res = 300)

grid.newpage()

top.vp <- viewport(layout=grid.layout(2, 3,
                                      widths=unit(c(1, 1, 1), c("null", "null", "null")),
                                      heights=unit(c(.5,5), c("null", "null", "null"))))

title1 <- viewport(layout.pos.col = 2, layout.pos.row = 1, name = "title1")
plot1 <- viewport(layout.pos.col = 1, layout.pos.row = 2, name = "plot1")
plot2 <- viewport(layout.pos.col = 2, layout.pos.row = 2, name = "plot2")
plot3 <- viewport(layout.pos.col = 3, layout.pos.row = 2, name = "plot3")

splot <- vpTree(top.vp, vpList(title1, plot1, plot2, plot3))

pushViewport(splot)

seekViewport("plot1")
grid.draw(ht_grob1)

seekViewport("plot2")
grid.draw(ht_grob2)

seekViewport("plot3")
grid.draw(ht_grob3)


dev.off()



