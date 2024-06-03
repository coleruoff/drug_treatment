source("source/cole_functions.R")
library(ComplexHeatmap)
library(circlize)

#GENESETS 1

curr_cell_line <- "MCF7"

cluster_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_cluster_signatures_200.rds"))


genesets_1 <- cluster_signatures
#GENESETS 2

ecoli_AMR_orthologs <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ecoli_AMR_genesets_orthologs.rds")
ecoli_AMR_consensus_orthologs <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ecoli_AMR_consensus_orthologs.rds")

ecoli_stress_orthologs <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ecoli_stress_orthologs.rds")
ecoli_stress_consensus_orthologs <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ecoli_stress_consensus_orthologs.rds")

yeast_stress_orthologs <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/yeast_upregulated_orthologs.rds")

genesets_2 <- list("yeast_stress_orthologs" = yeast_stress_orthologs)

# genesets_2 <- ecoli_stress_orthologs[grepl("down", names(ecoli_stress_orthologs))]

genesets_title <- "Yeast Stress Upregulated Orthologs"

#BACKGROUND GENES

background_genes <- unique(unlist(cluster_signatures))

##################################################################

or_mat <- compare_genesets_fisher(genesets1 = genesets_1, genesets2 = genesets_2, background_genes = background_genes, signif = F)

heatmap_matrix <- log(or_mat)

colnames(heatmap_matrix)

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- paste0(curr_cell_line,"_cluster",RACs[[curr_cell_line]],"_signature")

rac_ha <- HeatmapAnnotation(RAC = c(ifelse(colnames(heatmap_matrix) %in% clusters_of_interest,"RAC","Non-RAC")),
                            col = list(RAC = c("RAC" = "darkgreen", "Non-RAC" = "lightblue")))


ht <- Heatmap(heatmap_matrix, name="log(OR)", cluster_rows = F, cluster_columns = T,
              bottom_annotation = rac_ha, column_title = "Clusters", column_title_side = "bottom",
              row_title = "", column_names_rot = 45)


draw(ht, column_title = paste0(genesets_title, " ORA (", curr_cell_line, ")\n"), 
     column_title_gp = gpar(fontsize = 26),  padding = unit(c(2, 2, 2, 60), "mm"),
     heatmap_legend_side = "left", annotation_legend_side = "left")



rac_values <- as.vector(heatmap_matrix[,clusters_of_interest])
rest_values <- as.vector(heatmap_matrix[,!colnames(heatmap_matrix) %in% clusters_of_interest])

wilcox_res <- wilcox.test(rac_values,rest_values)
boxplot(rac_values,rest_values, names=c("RACs","Rest"), main=paste0(genesets_title, " ORA (", curr_cell_line, ") pval: ", sprintf("%.4f", wilcox_res$p.value)))



##################################################################

common_genes <- list()

for(trapnell_siganture in genesets_1){
  
  for(ecoli_signature in genesets_2){
    
    curr_overlap <- intersect(trapnell_siganture, ecoli_signature)
    
    common_genes <- append(common_genes, list(curr_overlap))
    
  }
  
}

temp <- unique(unlist(common_genes))

paste(temp,collapse=" ")


genes_to_use <- (unlist(common_genes))

# genes_to_use <- intersect(emergent_signatures[[1]], raj_de_genesets[[1]])

paste(genes_to_use,collapse=" ")

gene_list <- unique(c(unlist(cluster_signatures),unlist(ecoli_AMR_genesets_orthologs)))

go_enrich <- enrichGO(gene = genes_to_use,
                      universe = gene_list,
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)


functional_enrichment_plot_title <- paste0("Functional Enrichment of ", curr_cell_line, " Shared Genes with E Coli AMR Signatures")

dotplot(go_enrich)+
  ggtitle(functional_enrichment_plot_title)


