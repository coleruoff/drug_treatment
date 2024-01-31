genes_to_use <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_common_resistance_signature.rds")

go_enrich <- enrichGO(gene = genes_to_use,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10,)



dotplot(go_enrich,showCategory=30)+
  ggtitle("Functional Enrichment of Raj Common Resistance Signature")

##########

genes_to_use <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_breast_resistance_signature.rds")

go_enrich <- enrichGO(gene = genes_to_use,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)



dotplot(go_enrich,showCategory=30)+
  ggtitle("Functional Enrichment of Raj Breast Resistance Signature")

##########

genes_to_use <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_common_resistance_signature_union.rds")

go_enrich <- enrichGO(gene = genes_to_use,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)



dotplot(go_enrich,showCategory=30)+
  ggtitle("Functional Enrichment of Raj Common Union Resistance Signature")

##########
raj_resistant_cancer_type <- "common"

genes_to_use <- resistant_vs_control_de_genes <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/raj_", raj_resistant_cancer_type, "_resistant_vs_control_de_signature_500.rds"))

go_enrich <- enrichGO(gene = genes_to_use,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)



dotplot(go_enrich,showCategory=30)+
  ggtitle("Functional Enrichment of Raj Common Resistance Signature")

