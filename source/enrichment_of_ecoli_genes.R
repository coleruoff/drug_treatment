source("source/cole_functions.R")
library(Seurat)
library(org.EcK12.eg.db)
library(org.Hs.eg.db)

cell_lines <- c("A549","K562","MCF7")

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
##########################

#read in ecoli-human orthologs and save as RDS
upregulated_ecoli_human_orthologs_table <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/orthologs/upregulated_ecoli_human_orthologs_table.rds")
ecoli_human_orthologs_up <- unique(upregulated_ecoli_human_orthologs_table$HUMAN_SYMBOL[!is.na(upregulated_ecoli_human_orthologs_table$HUMAN_SYMBOL)])
ecoli_human_orthologs_up <- list("ecoli_human_orthologs_up"=ecoli_human_orthologs_up)


downregulated_ecoli_human_orthologs_table <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/orthologs/downregulated_ecoli_human_orthologs_table.rds")
ecoli_human_orthologs_down <- unique(downregulated_ecoli_human_orthologs_table$HUMAN_SYMBOL[!is.na(downregulated_ecoli_human_orthologs_table$HUMAN_SYMBOL)])
ecoli_human_orthologs_down <- list("ecoli_human_orthologs_down"=ecoli_human_orthologs_down)


#Find and remove intersecting genes
ecoli_human_orthologs_intersect <- intersect(ecoli_human_orthologs_up$ecoli_human_orthologs_up, ecoli_human_orthologs_down$ecoli_human_orthologs_down)

ecoli_human_orthologs_up$ecoli_human_orthologs_up <- ecoli_human_orthologs_up$ecoli_human_orthologs_up[!ecoli_human_orthologs_up$ecoli_human_orthologs_up %in% ecoli_human_orthologs_intersect]
ecoli_human_orthologs_down$ecoli_human_orthologs_down <- ecoli_human_orthologs_down$ecoli_human_orthologs_down[!ecoli_human_orthologs_down$ecoli_human_orthologs_down %in% ecoli_human_orthologs_intersect]

ecoli_human_orthologs_intersect <- list("ecoli_human_orthologs_intersect"=ecoli_human_orthologs_intersect)

saveRDS(ecoli_human_orthologs_intersect,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ecoli_human_orthologs_intersect.rds")
saveRDS(ecoli_human_orthologs_up,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ecoli_human_orthologs_up.rds")
saveRDS(ecoli_human_orthologs_down,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ecoli_human_orthologs_down.rds")

################################################################################
# Enrichment of ecoli genes of human-ecoli orthologs IN ecoli
################################################################################

#Get ecoli up genes
upregulated_ecoli_human_orthologs_table <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/orthologs/upregulated_ecoli_human_orthologs_table.rds")

ecoli_up <- upregulated_ecoli_human_orthologs_table %>% 
  filter(HUMAN_SYMBOL %in% ecoli_human_orthologs_up$ecoli_human_orthologs_up) %>% 
  pull(ECOLI_SYMBOL)


# Get ecoli down genes
downregulated_ecoli_human_orthologs_table <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/orthologs/downregulated_ecoli_human_orthologs_table.rds")

ecoli_down <- downregulated_ecoli_human_orthologs_table %>% 
  filter(HUMAN_SYMBOL %in% ecoli_human_orthologs_down$ecoli_human_orthologs_down) %>% 
  pull(ECOLI_SYMBOL)


# Get ecoli intersection genes
intersect_up <- upregulated_ecoli_human_orthologs_table %>% 
  filter(HUMAN_SYMBOL %in% ecoli_human_orthologs_intersect$ecoli_human_orthologs_intersect) %>% 
  pull(ECOLI_SYMBOL)
intersect_down <- downregulated_ecoli_human_orthologs_table %>% 
  filter(HUMAN_SYMBOL %in% ecoli_human_orthologs_intersect$ecoli_human_orthologs_intersect) %>% 
  pull(ECOLI_SYMBOL)


ecoli_intersect <- unique(c(intersect_up,intersect_down))

############################################################
# Ecoli up genes dotplot
genes <- ecoli_up
res <- enrichGO(genes,
         OrgDb =org.EcK12.eg.db,
         keyType = "SYMBOL")
dotplot(res)+
ggtitle("Enrichment of Upregulated E. coli Genes")

#Ecoli down genes dotplot
genes <- ecoli_down
res <- enrichGO(genes,
                OrgDb =org.EcK12.eg.db,
                keyType = "SYMBOL")
dotplot(res)+
  ggtitle("Enrichment of Downregulated E. coli Genes")

#Ecoli intersection genes dotplot
genes <- ecoli_intersect
res <- enrichGO(genes,
                OrgDb =org.EcK12.eg.db,
                keyType = "SYMBOL")
dotplot(res)+
  ggtitle("Enrichment of Up/Down Intersection E. coli Genes")

################################################################################
# Enrichment of ecoli human orthologs in HUMAN
################################################################################

# ecoli-human up orthologs dotplot
genes <- ecoli_human_orthologs_up$ecoli_human_orthologs_up
res <- enrichGO(genes,
                OrgDb =org.Hs.eg.db,
                keyType = "SYMBOL")

dotplot(res)+
  ggtitle("Enrichment of Upregulated E coli-Human Orthologs")


# ecoli-human down orthologs dotplot
genes <- ecoli_human_orthologs_down$ecoli_human_orthologs_down
res <- enrichGO(genes,
                OrgDb =org.Hs.eg.db,
                keyType = "SYMBOL")

dotplot(res)+
  ggtitle("Enrichment of Downregulated Ecoli-Human Orthologs")


# ecoli-human up/down intersection orthologs dotplot
genes <- ecoli_human_orthologs_intersect$ecoli_human_orthologs_intersect
res <- enrichGO(genes,
                OrgDb =org.Hs.eg.db,
                keyType = "SYMBOL")
dotplot(res)+
  ggtitle("Enrichment of Ecoli-Human Intersection orthologs genes")




write.table(sort(ecoli_human_orthologs_up$ecoli_human_orthologs_up),"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/test.txt",
            quote = F, row.names = F,col.names = F)


write.table(sort(ecoli_human_orthologs_down$ecoli_human_orthologs_down),"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/test.txt",
            quote = F, row.names = F,col.names = F)

################################################################################
# ECOLI ORTHOLOGS in RAC GSEA USING FIND MARKERS GENES
################################################################################
curr_cell_line <- "A549"
rac_type_de <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_cell_group_de.rds"))


ecoli_ortholog_t2g <- data.frame(rbind(cbind("ecoli_human_orthologs_up",ecoli_human_orthologs_up$ecoli_human_orthologs_up),cbind("ecoli_human_orthologs_down",ecoli_human_orthologs_down$ecoli_human_orthologs_down)))
# ecoli_ortholog_t2g <- data.frame(rbind(cbind("ecoli_human_orthologs_up",ecoli_human_orthologs_up$ecoli_human_orthologs_up)))
colnames(ecoli_ortholog_t2g) <- c("gs_name","human_gene_symbol")

##

rac_type1_genes <- rac_type_de %>% 
  filter(cluster == 1 & p_val_adj < 0.05) %>% 
  pull(avg_log2FC)

names(rac_type1_genes) <- rac_type_de %>% 
  filter(cluster == 1 & p_val_adj < 0.05) %>% 
  pull(gene)

gsea_res <- GSEA(sort(rac_type1_genes, decreasing = T),TERM2GENE = ecoli_ortholog_t2g)

gsea_res@result %>% 
  dplyr::select(Description,NES)

##

rac_type2_genes <- rac_type_de %>% 
  filter(cluster == 2 & p_val_adj < 0.05) %>% 
  pull(avg_log2FC)

names(rac_type2_genes) <- rac_type_de %>% 
  filter(cluster == 2 & p_val_adj < 0.05) %>% 
  pull(gene)

gsea_res <- GSEA(sort(rac_type2_genes, decreasing = T),TERM2GENE = ecoli_ortholog_t2g)

gsea_res@result %>% 
  dplyr::select(Description,NES)

##

rac_type0_genes <- rac_type_de %>% 
  filter(cluster == 0 & p_val_adj < 0.05) %>% 
  pull(avg_log2FC)

names(rac_type0_genes) <- rac_type_de %>% 
  filter(cluster == 0 & p_val_adj < 0.05) %>% 
  pull(gene)

gsea_res <- GSEA(sort(rac_type0_genes, decreasing = T),TERM2GENE = ecoli_ortholog_t2g)

gsea_res@result %>% 
  dplyr::select(Description,NES)


################################################################################
# ECOLI ORTHOLOGS in RAC ORA USING FIND MARKERS GENES
################################################################################
cell_line_universes <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cell_line_universes.rds")

rac_type1_genes_up <- rac_type_de %>% 
  filter(cluster == 1 & p_val_adj < 0.05 & avg_log2FC > 0) %>% 
  pull(gene)


rac_type1_genes_down <- rac_type_de %>% 
  filter(cluster == 1 & p_val_adj < 0.05 & avg_log2FC < 0) %>% 
  pull(gene)

class(ecoli_ortholog_t2g)

rac_up_ora <- enricher(rac_type1_genes_up, TERM2GENE = data.frame(ecoli_ortholog_t2g),
                       universe = cell_line_universes[[curr_cell_line]])

dotplot(rac_up_ora)


rac_down_ora <- enricher(rac_type1_genes_down, TERM2GENE = data.frame(ecoli_ortholog_t2g),
                         universe = cell_line_universes[[curr_cell_line]])


head(rac_down_ora)

rac_type1_genes <- rac_type_de %>% 
  filter(cluster == 1 & p_val_adj < 0.05) %>% 
  arrange(desc(avg_log2FC)) %>% 
  dplyr::select(gene,avg_log2FC) %>% 
  deframe()


ecoli_genes <- list()
ecoli_genes[["up"]] <- ecoli_human_orthologs_up$ecoli_human_orthologs_up
ecoli_genes[["down"]] <- ecoli_human_orthologs_down$ecoli_human_orthologs_down


rac_type1_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_type1_signatures.rds")

fgseaRes <- fgsea(ecoli_genes, stats = rac_type1_genes, nperm = 1000)


################################################################################
# Enrichment of shared genes between ecoli orthologs and rac type 1 global signatures
################################################################################

# Get rac de results
curr_cell_line <- "A549"
de_results <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/de_results/", curr_cell_line, "_cell_group_de.rds"))

# Shared genes of up rac and up ecoli-human orthologs
rac_type1_upregulated_genes <- de_results %>% 
  filter(cluster == 1 | p_val_adj < 0.05, avg_log2FC > 0) %>% 
  pull(gene)


genes <- sort(intersect(rac_type1_upregulated_genes, ecoli_human_orthologs_up$ecoli_human_orthologs_up))

res <- enrichGO(genes,
                OrgDb =org.Hs.eg.db,
                ont="BP",
                keyType = "SYMBOL")


dotplot(res)+
  ggtitle("Enrichment of Human orthologs genes shared with RAC type 1")



write.table(genes,"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/test.txt",
            quote = F, row.names = F,col.names = F)


# Shared genes of down rac and down ecoli-human orthologs
rac_type1_downregulated_genes <- de_results %>% 
  filter(cluster == 1 | p_val_adj < 0.05, avg_log2FC < 0) %>% 
  pull(gene)


genes <- sort(intersect(rac_type1_downregulated_genes, ecoli_human_orthologs_down$ecoli_human_orthologs_down))

res <- enrichGO(genes,
                OrgDb =org.Hs.eg.db,
                keyType = "SYMBOL")

dotplot(res)+
  ggtitle("Enrichment of downregulated Human orthologs genes shared with RAC type 1 down genes")




###################################################################################
# Supercluster signatures overlap with orthologs
########################################################################################


supercluster_signatures_up <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/type1_supercluster_signatures.rds")
supercluster_signatures_down <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/type1_supercluster_down_signatures.rds")

intersect(supercluster_signatures_up[[1]],ecoli_human_orthologs_up$ecoli_human_orthologs_up)
intersect(supercluster_signatures_up[[2]],ecoli_human_orthologs_up$ecoli_human_orthologs_up)

intersect(supercluster_signatures_down[[1]],ecoli_human_orthologs_down$ecoli_human_orthologs_down)
intersect(supercluster_signatures_down[[2]],ecoli_human_orthologs_down$ecoli_human_orthologs_down)


write.table(intersect(supercluster_signatures[[1]],ecoli_human_orthologs_up$ecoli_human_orthologs_up),"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/test.txt",
            quote = F, row.names = F,col.names = F)



write.table(intersect(supercluster_signatures[[1]],ecoli_human_orthologs_down$ecoli_human_orthologs_down),"/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/test.txt",
            quote = F, row.names = F,col.names = F)


###################################################################################
#
########################################################################################




