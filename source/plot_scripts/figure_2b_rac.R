setwd("/data/ruoffcj/projects/drug_treatment/")
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
source("source/cole_functions.R")

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, human_gene_symbol)

mp_t2g <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ith_meta_programs_t2g.rds")

specifc_mps <- c("MP39 Metal-response","MP31 Alveolar","MP29 NPC/OPC","MP28 Oligo normal","MP27 Oligo Progenitor","MP38 Glutathione","MP41 Unassigned","MP35 Hemato-related-I","MP37 Hemato-related-II","MP32 Skin-pigmentation","MP36 IG","MP16 MES (glioma)","MP15 EMT IV")

mp_t2g <- mp_t2g %>% 
  filter(!gs_name %in% specifc_mps)

cell_lines <- c("A549","K562","MCF7")

gene_universe_intersection <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cell_line_gene_universe_intersection.rds")
cell_line_universes <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cell_line_universes.rds")

################################################################################
# Plotting for RAC vs rest cells signature

all_barplots <- list()
for(curr_cell_line in cell_lines){
  y_label <- ""
  
  de_res <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line, "_global_rac_de.rds"))
  
  ranks <- de_res %>% 
    filter(cluster == "rac") %>% 
    arrange(desc(avg_log2FC)) %>% 
    dplyr::select(gene,avg_log2FC) %>% 
    deframe()
  
  # Hallmarks
  if(curr_cell_line == "K562"){
    y_label <- "Cancer Hallmarks"
  }
  
  hallmarks_gsea <- GSEA(ranks, TERM2GENE = m_t2g)
  
  df <- hallmarks_gsea@result %>% 
    filter(p.adjust < 0.05) %>% 
    arrange(desc(NES)) %>% 
    dplyr::select(Description,NES,p.adjust)
  
  
  hallmark_plot <- ggplot(df[1:5,])+
    geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES),width=0.7)+
    scale_fill_continuous(high="firebrick1",low="firebrick1")+
    xlab(y_label)+
    ylab("NES")+
    ggtitle(paste0(curr_cell_line))+
    theme(plot.title = element_text(size=30),
          axis.text = element_text(size=35),
          axis.title = element_text(size=35,face="bold"),
          legend.key.size=unit(1,'cm'),
          legend.title = element_text(size=20),
          legend.text = element_text(size=15))+
    NoLegend()+
    coord_flip()
  
  # Meta-programs
  if(curr_cell_line == "K562"){
    y_label <- "ITH Meta-Programs"
  }
  
  mp_gsea <- GSEA(ranks, TERM2GENE = mp_t2g)
  
  df <- mp_gsea@result %>% 
    filter(p.adjust < 0.05) %>% 
    arrange(desc(NES)) %>% 
    dplyr::select(Description,NES,p.adjust)
  
  
  mp_plot <- ggplot(df[1:5,])+
    geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES),width=0.7)+
    scale_fill_continuous(high="firebrick1",low="firebrick1")+
    xlab(y_label)+
    ylab("NES")+
    ggtitle(paste0(curr_cell_line))+
    theme(plot.title = element_text(size=30),
          axis.text = element_text(size=35),
          axis.title = element_text(size=35,face="bold"),
          legend.key.size=unit(1,'cm'),
          legend.title = element_text(size=20),
          legend.text = element_text(size=15))+
    NoLegend()+
    coord_flip()
  
  
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
  
  
  df <- go_gsea@result %>% 
    filter(p.adjust < 0.05) %>% 
    arrange(desc(NES)) %>% 
    dplyr::select(Description,NES,p.adjust)
  
  go_plot <- ggplot(df[1:5,])+
    geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES),width=0.7)+
    scale_fill_continuous(high="firebrick1",low="firebrick1")+
    xlab(y_label)+
    ylab("NES")+
    ggtitle(paste0(curr_cell_line))+
    theme(plot.title = element_text(size=30),
          axis.text = element_text(size=35),
          axis.title = element_text(size=35,face="bold"),
          legend.key.size=unit(1,'cm'),
          legend.title = element_text(size=20),
          legend.text = element_text(size=15))+
    NoLegend()+
    coord_flip()
  
  
  # all_barplots[["hallmarks"]] <- append(all_barplots[["hallmarks"]], list(hallmark_plot))
  all_barplots[["mps"]] <- append(all_barplots[["mps"]], list(mp_plot))
  all_barplots[["go"]] <- append(all_barplots[["go"]], list(go_plot))
  
}


# hallmark_plt <- ggarrange(plotlist = all_barplots$hallmarks, ncol = 3, common.legend = T,legend=c("right"))


mps_plt <- ggarrange(plotlist = all_barplots$mps, nrow = 3, common.legend = T,legend=c("none"))
# main_title <- paste0("ITH Meta-programs")
# mps_plt <- annotate_figure(mps_plt, top = text_grob(main_title, color = "black", face = "bold", size = 30))


go_plt <- ggarrange(plotlist = all_barplots$go, nrow = 3, common.legend = F,legend=c("none"))
# main_title <- paste0("GO Pathways")
# go_plt <- annotate_figure(go_plt, top = text_grob(main_title, color = "black", face = "bold", size = 30))


plots <- list(mps_plt, go_plt)

png(paste0("/data/ruoffcj/projects/drug_treatment/final_figures/figure_2b.png"),
width=38, height=22, units= "in", res = 300)


figure <- ggarrange(plotlist = plots, ncol=2, common.legend = T, legend=c("none"))

p <- annotate_figure(figure, left = text_grob("", rot = 90, vjust = 1, size=35, face="bold"),
                     bottom = text_grob("", size=35, face="bold"),
                     top=text_grob("Global RACs Enrichment", size=75, face="bold"))

print(p)

dev.off()















