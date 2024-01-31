library(org.Hs.eg.db)
library(clusterProfiler)
library(msigdbr)
library(ggpubr)
library(tidyverse)

curr_cell_line <- "MCF7"

curr_markers <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_cluster_all_markers_de.rds"))


RACs <- list(c(4,9,12,13,14,16,18),c(4,5,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- RACs[[curr_cell_line]]


dotplots <- list()

for(i in clusters_of_interest){
  
  cat(i,"\n")
  
  geneList <- curr_markers %>% 
    filter(cluster == i & avg_log2FC > 0 & p_val_adj < 0.05) %>% 
    arrange(desc(avg_log2FC)) %>% 
    pull(avg_log2FC)
  
  names(geneList) <- curr_markers %>% 
    filter(cluster == i & avg_log2FC > 0 & p_val_adj < 0.05) %>% 
    arrange(desc(avg_log2FC)) %>% 
    pull(gene)
  
  # if(length(geneList) > 200){
  #   geneList <- geneList[1:200]
  # }
  # 
  # ego3 <- gseGO(geneList     = geneList,
  #               OrgDb        = org.Hs.eg.db,
  #               ont          = "ALL",
  #               minGSSize    = 10,
  #               maxGSSize    = 500,
  #               pvalueCutoff = 0.05,
  #               verbose      = FALSE,
  #               keyType = "SYMBOL")
  # 
  # if(nrow(ego3) > 0){
  #   p <- dotplot(ego3) + ggtitle(paste0("Cluster ", i, " Markers GO GSEA"))
  #   
  #   dotplots <- append(dotplots,list(p))
  # }
  
  gene <- names(geneList)
  
  ego <- enrichGO(gene          = gene,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE,
                  keyType = "SYMBOL")
  
  if(nrow(ego) > 0){
    p <- dotplot(ego) + ggtitle(paste0("Cluster ", i))
    dotplots <- append(dotplots,list(p))
  }
  
}



plt <- egg::ggarrange(plots = dotplots, nrow=3)

main_title <- paste0(curr_cell_line, " Functional Enrichment of Resistance Activated Clusters")
plt <- annotate_figure(plt, top = text_grob(main_title, color = "black", face = "bold", size = 24))

plot(plt)




gene <- names(geneList)

ego <- enrichGO(gene          = gene,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE,
                keyType = "SYMBOL")


dotplot(ego) + ggtitle(paste0("Cluster ", i, " Markers GO ORA Enrichment"))


ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "ALL",
              minGSSize    = 10,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE,
              keyType = "SYMBOL")


dotplot(ego3) + ggtitle(paste0("Cluster ", i, " Markers GO GSEA"))

################################################################################




em <- enricher(gene, TERM2GENE=hallmark_t2g)

dotplot(em)


em2 <- GSEA(geneList, TERM2GENE = hallmark_t2g)


dotplot(em2)









