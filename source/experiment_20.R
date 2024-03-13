setwd("/data/ruoffcj/projects/drug_treatment/")
source("source/read_in_all_cell_lines.R")
source("/data/ruoffcj/projects/aucell_scoring/aucell_thresholding.R")
library(AUCell)
library(Seurat)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)

create_barplots <- function(ranks,signature_name){
  m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, human_gene_symbol)
  # hallmark_uni <- unique(c(unlist(all_signatures), m_t2g$human_gene_symbol))
  
  mp_t2g <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ith_meta_programs_t2g.rds")
  
  specifc_mps <- c("MP39 Metal-response","MP31 Alveolar","MP29 NPC/OPC","MP28 Oligo normal","MP27 Oligo Progenitor","MP38 Glutathione","MP41 Unassigned","MP35 Hemato-related-I","MP37 Hemato-related-II","MP32 Skin-pigmentation","MP36 IG","MP16 MES (glioma)","MP15 EMT IV")
  
  mp_t2g <- mp_t2g %>% 
    filter(!gs_name %in% specifc_mps)
  # Meta-programs
  mp_gsea <- GSEA(ranks, TERM2GENE = mp_t2g)
  
  df <- mp_gsea@result %>% 
    filter(p.adjust < 0.05) %>% 
    arrange(desc(NES)) %>% 
    dplyr::select(Description,NES,p.adjust)
  
  
  if(nrow(df)> 10){
    df <- df[1:10,]
  }
  mp_plot <- ggplot(df)+
    geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES))+
    scale_fill_continuous(high="red",low="pink")+
    xlab("Meta-Program")+
    ggtitle(paste0(curr_cell_line, " ", signature_name))+
    theme(plot.title = element_text(size=30),
          axis.text = element_text(size=15),
          axis.title = element_text(size=20,face="bold"))+
    coord_flip()
  
  
  # Hallmarks
  hallmark_gsea <- GSEA(ranks, TERM2GENE = m_t2g)
  
  df <- hallmark_gsea@result %>%
    filter(p.adjust < 0.05) %>%
    arrange(desc(NES)) %>%
    dplyr::select(Description,NES,p.adjust)
  
  
  hallmark_plot <- ggplot(df[1:10,])+
    geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES))+
    scale_fill_continuous(high="red",low="pink")+
    xlab("Cancer Hallmark")+
    ggtitle(paste0(curr_cell_line, " ", signature_name))+
    theme(plot.title = element_text(size=30),
          axis.text = element_text(size=15),
          axis.title = element_text(size=20,face="bold"))+
    coord_flip()
  
  
  # GO Pathways
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
  
  if(nrow(df)> 10){
    df <- df[1:10,]
  }
  
  go_plot <- ggplot(df)+
    geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES))+
    scale_fill_continuous(high="red",low="pink")+
    xlab("GO Pathways")+
    ggtitle(paste0(curr_cell_line, " ", signature_name))+
    theme(plot.title = element_text(size=30),
          axis.text = element_text(size=15),
          axis.title = element_text(size=20,face="bold"))+
    coord_flip()
  
  
  
  
  plot_list <- list(hallmark_plot,mp_plot,go_plot)
  
  
  p <- ggarrange(plotlist = plot_list,ncol=3,nrow=1)
  
  
  return(p)
}

set.seed(42)
curr_cell_line <- cell_lines[1]

cat(curr_cell_line, "\n")

data <- all_data[[curr_cell_line]]  

drug_classes <- as.character(unique(data$pathway_level_1))
drug_classes <- drug_classes[-which(drug_classes == "Vehicle")]
drug_classes <- drug_classes[-which(drug_classes == "Other")]


# drug_classes <- drug_classes[1:3]
all_plots <- list()
for(curr_class in drug_classes){
  cat(curr_class, "\n")
  
  drug_data <- data[,data$treatment_stage =="pre" | data$pathway_level_1 == curr_class]
  
  Idents(drug_data) <- drug_data$treatment_stage
  
  markers <- FindMarkers(drug_data, ident.1 = "post")
  
  temp <- markers %>% 
    arrange(desc(avg_log2FC)) %>% 
    rownames_to_column() %>% 
    dplyr::select(rowname,avg_log2FC) %>% 
    deframe()
  
  signature_name <- paste0(curr_class)
  all_plots <- append(all_plots,list(create_barplots(temp,signature_name)))
}


# saveRDS(all_plots, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/experiment_data/experiment_20_plotlist.rds")
all_plots <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/experiment_data/experiment_20_plotlist.rds")


length(all_plots)
p <- ggarrange(plotlist = all_plots,ncol=5,nrow=3)
p

