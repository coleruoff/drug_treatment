
A549_active <- c(4,9,13)

K562_active <- c(5)

MCF7_active <-  c(5,8,17)

signature_length <- 200

curr_cell_line <- "A549"
A549_cluster_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_cluster_signatures_",signature_length,".rds"))


curr_cell_line <- "K562"
K562_cluster_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_cluster_signatures_",signature_length,".rds"))


curr_cell_line <- "MCF7"
MCF7_cluster_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_cluster_signatures_",signature_length,".rds"))




genes_to_use <- A549_cluster_signatures[[3]]

go_enrich <- enrichGO(gene = genes_to_use,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)

dotplot(go_enrich,showCategory=20)+
  ggtitle("A549 Cluster 1 Functional Enrichment")
