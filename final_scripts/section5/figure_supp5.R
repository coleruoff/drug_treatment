args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"data/")
plotDirectory <- paste0(args[1],"figures/")

setwd(args[1])
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
library(tidyverse)
source("source/drug_treatment_functions.R")
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

################################################################################
# Plotting for ortholog signatures

orthologs <- list()

yeast_human_orthologs_up <- readRDS(paste0("genesets/yeast_human_orthologs_up.rds"))

ecoli_human_orthologs_up <- readRDS(paste0("genesets/ecoli_human_orthologs_up.rds"))


orthologs[["yeast_upregulated_orthologs"]] <- yeast_human_orthologs_up$yeast_human_orthologs_up
orthologs[["ecoli_upregulated_orthologs"]] <- ecoli_human_orthologs_up$ecoli_human_orthologs_up


all_results <- list()

for(i in 1:2){
  
  curr_geneset <- orthologs[[i]]
  
  # Hallmarks Enrichment
  hallmark_enrichment_results <- enricher(curr_geneset, TERM2GENE=m_t2g)
  
  # MPs Enrichment
  mp_enrichment_results <- enricher(curr_geneset, TERM2GENE=mp_t2g)
  
  # GO Functional Enrichment
  go_enrichment_results <- enrichGO(gene          = curr_geneset,
                                    OrgDb         = org.Hs.eg.db,
                                    ont           = "BP",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = .01,
                                    qvalueCutoff  = .05,
                                    readable      = TRUE,
                                    keyType = "SYMBOL")
  
  
  all_results[[i]] <- list(hallmark_enrichment_results,
                                              mp_enrichment_results,
                                              go_enrichment_results)
  
}


titles <- list("Cancer Hallmarks","ITH Meta-programs","GO Pathways")
all_plots <- list()

for(i in 1:3){
  
  # For each type of enrichment, create heatmap
  curr_results <- list(all_results[[1]][[i]],all_results[[2]][[i]])
  names(curr_results) <- c("Yeast Orthologs","E. Coli Orthologs")
  
  heatmap <- create_enrichment_heatmap(curr_results, titles[i])
  
  # Add heatmap to curr cell line list of heatmaps
  all_plots <- append(all_plots, heatmap)
}



curr_title <- paste0("Yeast and E. Coli Resistance Orthologs Enrichment")

ht_grob1 = grid.grabExpr(draw(all_plots[[1]], padding = unit(c(0, 80, 0, 0), "mm")))
ht_grob2 = grid.grabExpr(draw(all_plots[[2]], padding = unit(c(0, 80, 0, 0), "mm")))
ht_grob3 = grid.grabExpr(draw(all_plots[[3]], padding = unit(c(0, 80, 0, 0), "mm")))



png(paste0(plotDirectory, "figure_5_supp.png"),
    width=36, height=15, units= "in", res = 300)

grid.newpage()

top.vp <- viewport(layout=grid.layout(2, 3,
                                      widths=unit(c(1,1,1), c("null", "null", "null")),
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

seekViewport("title1")
grid.text(curr_title, gp = gpar(fontsize = 40))

dev.off()


