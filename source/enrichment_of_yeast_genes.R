source("source/cole_functions.R")
library(Seurat)
library(org.EcK12.eg.db)
library(org.Hs.eg.db)

cell_lines <- c("A549","K562","MCF7")

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
##########################

#read in yeast-human orthologs and save as RDS
upregulated_yeast_human_orthologs_table <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/orthologs/upregulated_yeast_human_orthologs_table.rds")
yeast_human_orthologs_up <- unique(upregulated_yeast_human_orthologs_table$HUMAN_SYMBOL[!is.na(upregulated_yeast_human_orthologs_table$HUMAN_SYMBOL)])
yeast_human_orthologs_up <- list("yeast_human_orthologs_up"=yeast_human_orthologs_up)
saveRDS(yeast_human_orthologs_up,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/yeast_human_orthologs_up.rds")

downregulated_yeast_human_orthologs_table <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/orthologs/downregulated_yeast_human_orthologs_table.rds")
yeast_human_orthologs_down <- unique(downregulated_yeast_human_orthologs_table$HUMAN_SYMBOL[!is.na(downregulated_yeast_human_orthologs_table$HUMAN_SYMBOL)])
yeast_human_orthologs_down <- list("yeast_human_orthologs_down"=yeast_human_orthologs_down)
saveRDS(yeast_human_orthologs_down,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/yeast_human_orthologs_down.rds")

################################################################################
# yeast ORTHOLOGS GSEA USING FIND MARKERS GENES
################################################################################
curr_cell_line <- "A549"
rac_type_de <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_cell_group_de.rds"))


yeast_ortholog_t2g <- data.frame(rbind(cbind("yeast_human_orthologs_up",yeast_human_orthologs_up$yeast_human_orthologs_up),cbind("yeast_human_orthologs_down",yeast_human_orthologs_down$yeast_human_orthologs_down)))
colnames(yeast_ortholog_t2g) <- c("gs_name","human_gene_symbol")

##

rac_type1_genes <- rac_type_de %>% 
  filter(cluster == 1 & p_val_adj < 0.05) %>% 
  pull(avg_log2FC)

names(rac_type1_genes) <- rac_type_de %>% 
  filter(cluster == 1 & p_val_adj < 0.05) %>% 
  pull(gene)

gsea_res <- GSEA(sort(rac_type1_genes, decreasing = T),TERM2GENE = yeast_ortholog_t2g)

gsea_res@result %>% 
  dplyr::select(Description,NES)

##

rac_type2_genes <- rac_type_de %>% 
  filter(cluster == 2 & p_val_adj < 0.05) %>% 
  pull(avg_log2FC)

names(rac_type2_genes) <- rac_type_de %>% 
  filter(cluster == 2 & p_val_adj < 0.05) %>% 
  pull(gene)

gsea_res <- GSEA(sort(rac_type2_genes, decreasing = T),TERM2GENE = yeast_ortholog_t2g)

gsea_res@result %>% 
  dplyr::select(Description,NES)

##

rac_type0_genes <- rac_type_de %>% 
  filter(cluster == 0 & p_val_adj < 0.05) %>% 
  pull(avg_log2FC)

names(rac_type0_genes) <- rac_type_de %>% 
  filter(cluster == 0 & p_val_adj < 0.05) %>% 
  pull(gene)

gsea_res <- GSEA(sort(rac_type0_genes, decreasing = T),TERM2GENE = yeast_ortholog_t2g)

gsea_res@result %>% 
  dplyr::select(Description,NES)


################################################################################
# yeast ORTHOLOGS ORA USING FIND MARKERS GENES
################################################################################
cell_line_universes <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cell_line_universes.rds")

rac_type1_genes_up <- rac_type_de %>% 
  filter(cluster == 1 & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
  pull(gene)


rac_type1_genes_down <- rac_type_de %>% 
  filter(cluster == 1 & p_val_adj < 0.05 & avg_log2FC < 0) %>% 
  pull(gene)

class(yeast_ortholog_t2g)

rac_up_ora <- enricher(rac_type1_genes_up, TERM2GENE = data.frame(yeast_ortholog_t2g),
                       universe = cell_line_universes[[curr_cell_line]])

dotplot(rac_up_ora)


rac_down_ora <- enricher(rac_type1_genes_down, TERM2GENE = data.frame(yeast_ortholog_t2g),
                         universe = cell_line_universes[[curr_cell_line]])


head(rac_down_ora)

rac_type1_genes <- rac_type_de %>% 
  filter(cluster == 1 & p_val_adj < 0.05) %>% 
  arrange(desc(avg_log2FC)) %>% 
  dplyr::select(gene,avg_log2FC) %>% 
  deframe()


yeast_genes <- list()
yeast_genes[["up"]] <- yeast_human_orthologs_up$yeast_human_orthologs_up
yeast_genes[["down"]] <- yeast_human_orthologs_down$yeast_human_orthologs_down


rac_type1_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_type1_signatures.rds")

fgseaRes <- fgsea(yeast_genes, stats = rac_type1_genes, nperm = 1000)





################################################################################
# Enrichment of yeast genes of human-yeast orthologs IN yeast
upregulated_yeast_human_orthologs_table <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/orthologs/upregulated_yeast_human_orthologs_table.rds")
yeast_up <- unique(upregulated_yeast_human_orthologs_table$yeast_SYMBOL[!is.na(upregulated_yeast_human_orthologs_table$yeast_SYMBOL)])


downregulated_yeast_human_orthologs_table <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/orthologs/downregulated_yeast_human_orthologs_table.rds")
yeast_down <- unique(downregulated_yeast_human_orthologs_table$yeast_SYMBOL[!is.na(downregulated_yeast_human_orthologs_table$yeast_SYMBOL)])

# 
# 
# yeast_genename_geneid_conversion_table <- read.table("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/yeast_data/yeast_genename_geneid_conversion_table.txt", fill = T)
# colnames(yeast_genename_geneid_conversion_table) <- yeast_genename_geneid_conversion_table[1,]
# yeast_genename_geneid_conversion_table <- yeast_genename_geneid_conversion_table[-1,]
# genes <- unique(yeast_genename_geneid_conversion_table$gene_name[yeast_genename_geneid_conversion_table$gene_id %in% up_genes])

genes <- yeast_up

res <- enrichGO(genes,
                OrgDb =org.EcK12.eg.db,
                keyType = "SYMBOL")


dotplot(res)+
  ggtitle("Enrichment of Upregulated E. coli genes")


genes <- yeast_down

res <- enrichGO(genes,
                OrgDb =org.EcK12.eg.db,
                keyType = "SYMBOL")


dotplot(res)+
  ggtitle("Enrichment of Downregulated E. coli genes")


################################################################################
# Enrichment of yeast human orthologs in HUMAN
################################################################################

genes <- yeast_human_orthologs_up$yeast_human_orthologs_up

res <- enrichGO(genes,
                OrgDb =org.Hs.eg.db,
                ont="ALL",
                keyType = "SYMBOL")

dotplot(res)+
  ggtitle("Enrichment of Upregulated Yeast-Human Orthologs")

yeast_res <- genesets_characterization(list("yeast_up"=genes))



genes <- yeast_human_orthologs_down$yeast_human_orthologs_down


res <- enrichGO(genes,
                OrgDb =org.Hs.eg.db,
                keyType = "SYMBOL")

dotplot(res)+
  ggtitle("Enrichment of yeast-Human Down orthologs genes")



intersect(yeast_human_orthologs_up$yeast_human_orthologs_up,yeast_human_orthologs_down$yeast_human_orthologs_down)

write.table(sort(yeast_human_orthologs_up$yeast_human_orthologs_up),"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/test.txt",
            quote = F, row.names = F,col.names = F)


write.table(sort(yeast_human_orthologs_down$yeast_human_orthologs_down),"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/test.txt",
            quote = F, row.names = F,col.names = F)


##########################
# Enrichment of shared genes between yeast orthologs and rac type 1 global signatures

curr_cell_line <- "A549"
de_results <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_cell_group_de.rds"))

rac_type1_upregulated_genes <- de_results %>% 
  filter(cluster == 1 | p_val_adj < 0.05, avg_log2FC > 0) %>% 
  pull(gene)


genes <- sort(intersect(rac_type1_upregulated_genes, yeast_human_orthologs_up$yeast_human_orthologs_up))


res <- enrichGO(genes,
                OrgDb =org.Hs.eg.db,
                ont = "BP",
                keyType = "SYMBOL")


dotplot(res)+
  ggtitle("Enrichment of Human orthologs genes shared with RAC type 1")


rac_type1_downregulated_genes <- de_results %>% 
  filter(cluster == 1 | p_val_adj < 0.05, avg_log2FC < 0) %>% 
  pull(gene)


genes <- sort(intersect(rac_type1_downregulated_genes, yeast_human_orthologs_down$yeast_human_orthologs_down))


res <- enrichGO(genes,
                OrgDb =org.Hs.eg.db,
                ont = "BP",
                keyType = "SYMBOL")



dotplot(res)+
  ggtitle("Enrichment of Human orthologs genes shared with RAC type 1")



sort(rac_type1_upregulated_genes)




rac_type1_genes <- de_results %>% 
  filter(cluster == 1 & p_val_adj < 0.05) %>% 
  pull(avg_log2FC)

names(rac_type1_genes) <- de_results %>% 
  filter(cluster == 1 & p_val_adj < 0.05) %>% 
  pull(gene)



ego3 <- gseGO(geneList     = sort(rac_type1_genes, decreasing = T),
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE,
              keyType = "SYMBOL")



ego3@result %>% 
  dplyr::select(Description,NES) %>% 
  arrange(desc(NES))

yeast_ortholog_t2g <- data.frame(rbind(cbind("yeast_human_orthologs_up",yeast_human_orthologs_up$yeast_human_orthologs_up),cbind("yeast_human_orthologs_down",yeast_human_orthologs_down$yeast_human_orthologs_down)))
colnames(yeast_ortholog_t2g) <- c("gs_name","human_gene_symbol")



gsea_res <- GSEA(sort(rac_type1_genes, decreasing = T),TERM2GENE = yeast_ortholog_t2g)


gsea_res@result %>% 
  dplyr::select(Description,NES)

#GO ORA Type 1 UP
rac_type1_genes_up <- de_results %>% 
  filter(cluster == 1 & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
  arrange(desc(avg_log2FC)) %>% 
  pull(gene)


ego <- enrichGO(gene          = rac_type1_genes_up[1:200],
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE,
                keyType = "SYMBOL")
dotplot(ego)



###########

source("source/cole_functions.R")
library(Seurat)
library(org.EcK12.eg.db)
library(org.Hs.eg.db)

cell_lines <- c("A549","K562","MCF7")

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
##########################

#read in yeast-human orthologs and save as RDS
upregulated_yeast_human_orthologs_table <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/orthologs/upregulated_yeast_human_orthologs_table.rds")
yeast_human_orthologs_up <- unique(upregulated_yeast_human_orthologs_table$HUMAN_SYMBOL[!is.na(upregulated_yeast_human_orthologs_table$HUMAN_SYMBOL)])
yeast_human_orthologs_up <- list("yeast_human_orthologs_up"=yeast_human_orthologs_up)


downregulated_yeast_human_orthologs_table <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/orthologs/downregulated_yeast_human_orthologs_table.rds")
yeast_human_orthologs_down <- unique(downregulated_yeast_human_orthologs_table$HUMAN_SYMBOL[!is.na(downregulated_yeast_human_orthologs_table$HUMAN_SYMBOL)])
yeast_human_orthologs_down <- list("yeast_human_orthologs_down"=yeast_human_orthologs_down)


#Find and remove intersecting genes
yeast_human_orthologs_intersect <- intersect(yeast_human_orthologs_up$yeast_human_orthologs_up, yeast_human_orthologs_down$yeast_human_orthologs_down)

yeast_human_orthologs_up$yeast_human_orthologs_up <- yeast_human_orthologs_up$yeast_human_orthologs_up[!yeast_human_orthologs_up$yeast_human_orthologs_up %in% yeast_human_orthologs_intersect]
yeast_human_orthologs_down$yeast_human_orthologs_down <- yeast_human_orthologs_down$yeast_human_orthologs_down[!yeast_human_orthologs_down$yeast_human_orthologs_down %in% yeast_human_orthologs_intersect]

yeast_human_orthologs_intersect <- list("yeast_human_orthologs_intersect"=yeast_human_orthologs_intersect)

saveRDS(yeast_human_orthologs_intersect,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/yeast_human_orthologs_intersect.rds")
saveRDS(yeast_human_orthologs_up,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/yeast_human_orthologs_up.rds")
saveRDS(yeast_human_orthologs_down,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/yeast_human_orthologs_down.rds")

################################################################################
# Enrichment of yeast genes of human-yeast orthologs IN yeast
################################################################################

#Get yeast up genes
upregulated_yeast_human_orthologs_table <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/orthologs/upregulated_yeast_human_orthologs_table.rds")

yeast_up <- upregulated_yeast_human_orthologs_table %>% 
  filter(HUMAN_SYMBOL %in% yeast_human_orthologs_up$yeast_human_orthologs_up) %>% 
  pull()


# Get yeast down genes
downregulated_yeast_human_orthologs_table <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/orthologs/downregulated_yeast_human_orthologs_table.rds")

yeast_down <- downregulated_yeast_human_orthologs_table %>% 
  filter(HUMAN_SYMBOL %in% yeast_human_orthologs_down$yeast_human_orthologs_down) %>% 
  pull(yeast_SYMBOL)


# Get yeast intersection genes
intersect_up <- upregulated_yeast_human_orthologs_table %>% 
  filter(HUMAN_SYMBOL %in% yeast_human_orthologs_intersect$yeast_human_orthologs_intersect) %>% 
  pull(yeast_SYMBOL)
intersect_down <- downregulated_yeast_human_orthologs_table %>% 
  filter(HUMAN_SYMBOL %in% yeast_human_orthologs_intersect$yeast_human_orthologs_intersect) %>% 
  pull(yeast_SYMBOL)


yeast_intersect <- unique(c(intersect_up,intersect_down))

############################################################
# yeast up genes dotplot
genes <- yeast_up
res <- enrichGO(genes,
                OrgDb =org.EcK12.eg.db,
                keyType = "SYMBOL")
dotplot(res)+
  ggtitle("Enrichment of Upregulated yeast Genes")

#yeast down genes dotplot
genes <- yeast_down
res <- enrichGO(genes,
                OrgDb =org.EcK12.eg.db,
                keyType = "SYMBOL")
dotplot(res)+
  ggtitle("Enrichment of Downregulated yeast Genes")

#yeast intersection genes dotplot
genes <- yeast_intersect
res <- enrichGO(genes,
                OrgDb =org.EcK12.eg.db,
                keyType = "SYMBOL")
dotplot(res)+
  ggtitle("Enrichment of Up/Down Intersection yeast Genes")

################################################################################
# Enrichment of yeast human orthologs in HUMAN
################################################################################

# yeast-human up orthologs dotplot
genes <- yeast_human_orthologs_up$yeast_human_orthologs_up
res <- enrichGO(genes,
                OrgDb =org.Hs.eg.db,
                ont = "ALL",
                keyType = "SYMBOL")

dotplot(res)+
  ggtitle("Enrichment of Upregulated Yeast-Human Orthologs")


# yeast-human down orthologs dotplot
genes <- yeast_human_orthologs_down$yeast_human_orthologs_down
res <- enrichGO(genes,
                OrgDb =org.Hs.eg.db,
                ont='ALL',
                keyType = "SYMBOL")

dotplot(res)+
  ggtitle("Enrichment of Downregulated Yeast-Human Orthologs")


# yeast-human up/down intersection orthologs dotplot
genes <- yeast_human_orthologs_intersect$yeast_human_orthologs_intersect
res <- enrichGO(genes,
                OrgDb =org.Hs.eg.db,
                keyType = "SYMBOL")
dotplot(res)+
  ggtitle("Enrichment of Yeast-Human Intersection orthologs genes")




write.table(sort(yeast_human_orthologs_up$yeast_human_orthologs_up),"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/test.txt",
            quote = F, row.names = F,col.names = F)


write.table(sort(yeast_human_orthologs_down$yeast_human_orthologs_down),"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/test.txt",
            quote = F, row.names = F,col.names = F)

################################################################################
# yeast ORTHOLOGS in RAC GSEA USING FIND MARKERS GENES
################################################################################
curr_cell_line <- "A549"
rac_type_de <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_cell_group_de.rds"))


yeast_ortholog_t2g <- data.frame(rbind(cbind("yeast_human_orthologs_up",yeast_human_orthologs_up$yeast_human_orthologs_up),cbind("yeast_human_orthologs_down",yeast_human_orthologs_down$yeast_human_orthologs_down)))
# yeast_ortholog_t2g <- data.frame(rbind(cbind("yeast_human_orthologs_up",yeast_human_orthologs_up$yeast_human_orthologs_up)))
colnames(yeast_ortholog_t2g) <- c("gs_name","human_gene_symbol")

##

rac_type1_genes <- rac_type_de %>% 
  filter(cluster == 1 & p_val_adj < 0.05) %>% 
  pull(avg_log2FC)

names(rac_type1_genes) <- rac_type_de %>% 
  filter(cluster == 1 & p_val_adj < 0.05) %>% 
  pull(gene)

gsea_res <- GSEA(sort(rac_type1_genes, decreasing = T),TERM2GENE = yeast_ortholog_t2g)

gsea_res@result %>% 
  dplyr::select(Description,NES)

##

rac_type2_genes <- rac_type_de %>% 
  filter(cluster == 2 & p_val_adj < 0.05) %>% 
  pull(avg_log2FC)

names(rac_type2_genes) <- rac_type_de %>% 
  filter(cluster == 2 & p_val_adj < 0.05) %>% 
  pull(gene)

gsea_res <- GSEA(sort(rac_type2_genes, decreasing = T),TERM2GENE = yeast_ortholog_t2g)

gsea_res@result %>% 
  dplyr::select(Description,NES)

##

rac_type0_genes <- rac_type_de %>% 
  filter(cluster == 0 & p_val_adj < 0.05) %>% 
  pull(avg_log2FC)

names(rac_type0_genes) <- rac_type_de %>% 
  filter(cluster == 0 & p_val_adj < 0.05) %>% 
  pull(gene)

gsea_res <- GSEA(sort(rac_type0_genes, decreasing = T),TERM2GENE = yeast_ortholog_t2g)

gsea_res@result %>% 
  dplyr::select(Description,NES)


################################################################################
# yeast ORTHOLOGS in RAC ORA USING FIND MARKERS GENES
################################################################################
cell_line_universes <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cell_line_universes.rds")

rac_type1_genes_up <- rac_type_de %>% 
  filter(cluster == 1 & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
  pull(gene)


rac_type1_genes_down <- rac_type_de %>% 
  filter(cluster == 1 & p_val_adj < 0.05 & avg_log2FC < 0) %>% 
  pull(gene)

class(yeast_ortholog_t2g)

rac_up_ora <- enricher(rac_type1_genes_up, TERM2GENE = data.frame(yeast_ortholog_t2g),
                       universe = cell_line_universes[[curr_cell_line]])

dotplot(rac_up_ora)


rac_down_ora <- enricher(rac_type1_genes_down, TERM2GENE = data.frame(yeast_ortholog_t2g),
                         universe = cell_line_universes[[curr_cell_line]])


head(rac_down_ora)

rac_type1_genes <- rac_type_de %>% 
  filter(cluster == 1 & p_val_adj < 0.05) %>% 
  arrange(desc(avg_log2FC)) %>% 
  dplyr::select(gene,avg_log2FC) %>% 
  deframe()


yeast_genes <- list()
yeast_genes[["up"]] <- yeast_human_orthologs_up$yeast_human_orthologs_up
yeast_genes[["down"]] <- yeast_human_orthologs_down$yeast_human_orthologs_down


rac_type1_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_type1_signatures.rds")

fgseaRes <- fgsea(yeast_genes, stats = rac_type1_genes, nperm = 1000)


################################################################################
# Enrichment of shared genes between yeast orthologs and rac type 1 global signatures
################################################################################

shared_with_yeast <- list()

# Get rac de results
curr_cell_line <- "MCF7"
de_results <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_cell_group_de.rds"))

# Shared genes of up rac and up yeast-human orthologs
rac_type1_upregulated_genes <- de_results %>% 
  filter(cluster == 1 & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
  pull(gene)


genes <- sort(intersect(rac_type1_upregulated_genes, yeast_human_orthologs_up$yeast_human_orthologs_up))

shared_with_yeast[[curr_cell_line]] <- sort(genes)

find_consensus_geneset(shared_with_yeast, 3)

res <- enrichGO(genes,
                OrgDb =org.Hs.eg.db,
                keyType = "SYMBOL")

dotplot(res)+
  ggtitle("Enrichment of Human orthologs genes shared with RAC type 1")



write.table(genes,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/test.txt",
            quote = F, row.names = F,col.names = F)


# Shared genes of down rac and down yeast-human orthologs
rac_type1_downregulated_genes <- de_results %>% 
  filter(cluster == 1 & p_val_adj < 0.05 & avg_log2FC < 0) %>% 
  pull(gene)


genes <- sort(intersect(rac_type1_downregulated_genes, yeast_human_orthologs_down$yeast_human_orthologs_down))

res <- enrichGO(genes,
                OrgDb =org.Hs.eg.db,
                ont = "ALL",
                keyType = "SYMBOL")

dotplot(res)+
  ggtitle("Enrichment of downregulated Human orthologs genes shared with RAC type 1 down genes")




###################################################################################
# Supercluster signatures overlap with orthologs
########################################################################################

type1_supercluster_upregulated_genesets <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/type1_superclusters_upregulated_genesets.rds")


shared_supercluster1_genes <- intersect(type1_supercluster_upregulated_genesets[[1]],yeast_human_orthologs_up$yeast_human_orthologs_up)
shared_supercluster2_genes <- intersect(type1_supercluster_upregulated_genesets[[2]],yeast_human_orthologs_up$yeast_human_orthologs_up)

intersect(shared_supercluster1_genes,shared_supercluster2_genes)


length(shared_supercluster2_genes) <- length(shared_supercluster1_genes)
yeast_supercluser_shared_genes <- list(shared_supercluster1_genes,shared_supercluster2_genes)
names(yeast_supercluser_shared_genes) <- c("shared_yeast_supercluster1_genes","shared_yeast_supercluster2_genes")


write.csv(yeast_supercluser_shared_genes,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/fi/test.csv", 
          row.names = F)




genes <- shared_supercluster2_genes

res <- enrichGO(genes,
                OrgDb =org.Hs.eg.db,
                ont = "ALL",
                keyType = "SYMBOL")

dotplot(res)+
  ggtitle("Enrichment of Shared Upregulated Yeast-Human Orthologs and Upregulated Supercluster 1 Genes")



###################################################################################
#s
########################################################################################












