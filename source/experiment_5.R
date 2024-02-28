library(forcats)
library(org.Hs.eg.db)
library(clusterProfiler)
library(msigdbr)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

cell_lines <- c("A549","K562","MCF7")

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")

# Create Global RAC Signatures for each cell line
rac_signatures <- list()
for(curr_cell_line in cell_lines){
  cat(curr_cell_line,"\n")
  
  #Read in cell line data
  data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/", curr_cell_line, "_processed_filtered.rds"))
  
  #read in DR signature scores and set active cells
  scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
  threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
  active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
  clusters_of_interest <- RACs[[curr_cell_line]]
  
  #Add metadata for RAC and Cell Group
  data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, "1", ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
  data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, paste0(data$Cluster, "_1"), ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), paste0(data$Cluster, "_2"), paste0(data$Cluster, "_0"))), col.name = "cell_cluster_group")
  
  Idents(data) <- data$rac
  
  #Find global RAC markers for current cell line
  # de_res <- FindAllMarkers(data)
  # saveRDS(de_res, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line, "_global_rac_de.rds"))
  
  de_res <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line, "_global_rac_de.rds"))
  
  curr_rac_up_genes <- de_res %>% 
    filter(cluster == "rac" & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
    arrange(desc(avg_log2FC)) %>% 
    pull(gene)
  
  rac_signatures[[curr_cell_line]] <- curr_rac_up_genes[1:500]
}

saveRDS(rac_signatures, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_rac_signatures.rds")


m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, human_gene_symbol)
# hallmark_uni <- unique(c(unlist(all_signatures), m_t2g$human_gene_symbol))

mp_t2g <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ith_meta_programs_t2g.rds")

specifc_mps <- c("MP39 Metal-response","MP31 Alveolar","MP29 NPC/OPC","MP28 Oligo normal","MP27 Oligo Progenitor","MP38 Glutathione","MP41 Unassigned","MP35 Hemato-related-I","MP37 Hemato-related-II","MP32 Skin-pigmentation","MP36 IG","MP16 MES (glioma)","MP15 EMT IV")

mp_t2g <- mp_t2g %>% 
  filter(!gs_name %in% specifc_mps)


cell_line_universes <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cell_line_universes.rds")


curr_cell_line <- "A549"
curr_geneset <- rac_signatures[[curr_cell_line]]
universe_to_use <- cell_line_universes[[curr_cell_line]]

# Hallmarks Enrichment
hallmark_enrichment_results <- enricher(curr_geneset, TERM2GENE=m_t2g,
                                        universe = universe_to_use)


hallmark_plot <- barplot(hallmark_enrichment_results)


# MPs Enrichment
mp_enrichment_results <- enricher(curr_geneset, TERM2GENE=mp_t2g,
                                  universe = universe_to_use)



mp_plot <- barplot(mp_enrichment_results)

# GO Functional Enrichment
ego <- enrichGO(gene          = curr_geneset,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = .01,
                qvalueCutoff  = .05,
                readable      = TRUE,
                keyType = "SYMBOL",
                universe = universe_to_use)


go_plot <- barplot(ego)


plot_list <- list(hallmark_plot,mp_plot,go_plot)


p <- ggarrange(plotlist=plot_list,ncol=3)

annotate_figure(p, top=curr_cell_line)






curr_cell_line <- "A549"
de_res <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/",curr_cell_line, "_global_rac_de.rds"))

ranks <- de_res %>% 
  filter(cluster == "rac") %>% 
  arrange(desc(avg_log2FC)) %>% 
  dplyr::select(gene,avg_log2FC) %>% 
  deframe()



# Meta-programs
mp_gsea <- GSEA(ranks, TERM2GENE = mp_t2g)

df <- mp_gsea@result %>% 
  filter(p.adjust < 0.05) %>% 
  arrange(desc(NES)) %>% 
  dplyr::select(Description,NES,p.adjust)



mp_plot <- ggplot(df[1:10,])+
  geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES))+
  scale_fill_continuous(high="red",low="pink")+
  xlab("")+
  ggtitle(paste0("ITH Meta-programs"))+
  theme(plot.title = element_text(size=30),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20,face="bold"))+
  coord_flip()


# Hallmarks
# hallmark_gsea <- GSEA(ranks, TERM2GENE = m_t2g)
# 
# df <- hallmark_gsea@result %>% 
#   filter(p.adjust < 0.05) %>% 
#   arrange(desc(NES)) %>% 
#   dplyr::select(Description,NES,p.adjust)
# 
# 
# 
# 
# hallmark_plot <- ggplot(df[1:10,])+
#   geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES))+
#   scale_fill_continuous(high="red",low="pink")+
#   xlab("Cancer Hallmark")+
#   ggtitle(paste0(curr_cell_line, " RAC Enrichment"))+
#   theme(plot.title = element_text(size=30),
#         axis.text = element_text(size=15),
#         axis.title = element_text(size=20,face="bold"))+
#   coord_flip()


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

go_plot <- ggplot(df[1:10,])+
  geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES))+
  scale_fill_continuous(high="red",low="pink")+
  xlab("")+
  ggtitle(paste0("GO Pathways"))+
  theme(plot.title = element_text(size=30),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20,face="bold"))+
  coord_flip()




plot_list <- list(mp_plot,go_plot)

figure <- ggarrange(plotlist = plot_list, nrow=1, common.legend = T,legend=c("right"))

p <- annotate_figure(figure, left = text_grob("", rot = 90, vjust = 1, size=35, face="bold"),
                     bottom = text_grob("", size=35, face="bold"),
                     top=text_grob(paste0(curr_cell_line,"RAC Signatures Enrichment"), size=50, face="bold"))


p

