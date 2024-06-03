library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(fgsea)

create_GSEA_matrix <- function(genesets_with_ranks, genesets2){
  
  #Remove genes from genesets2 that are not in the ranked gene list
  # for(i in 1:length(genesets2)){
  #   genesets2[[i]] <- genesets2[[i]][genesets2[[i]] %in% names(genesets_with_ranks[[1]])]
  # }
  
  gsea_mat <- matrix(NA, nrow=length(genesets_with_ranks), ncol=length(genesets2))
  rownames(gsea_mat) <- names(genesets_with_ranks)
  colnames(gsea_mat) <- names(genesets2)

  gsea_results <- list()
  
  for(i in 1:length(genesets_with_ranks)){

    fgseaRes <- fgsea(pathways = genesets2, 
                      stats    = genesets_with_ranks[[i]],
                      minSize  = 10,
                      maxSize  = max(lengths(genesets2))+1)
    
    gsea_results <- append(gsea_results, list(fgseaRes))
    
    # gsea_mat[i,] <- fgseaRes$NES
    
    gsea_mat[i,(colnames(gsea_mat) %in% fgseaRes$pathway)] <- fgseaRes$NES
    gsea_mat[i,!(colnames(gsea_mat) %in% fgseaRes$pathway)] <- 0
    
    
  }
  
  names(gsea_results) <- sapply(names(genesets_with_ranks), FUN = function(x){gsub("_fc_results", "", x)})

  gsea_mat[is.na(gsea_mat)] <- 0
  
  gsea_mat <- t(gsea_mat)
  return(list(gsea_mat,gsea_results))
  
}

# Compares two genesets by ORA using Fisher test
compare_genesets_fisher <- function(genesets1, genesets2, background_genes, signif = F){
  
  or_mat <- matrix(NA, nrow=length(genesets1), ncol=length(genesets2))
  
  for(i in 1:length(genesets1)){
    
    for(j in 1:length(genesets2)){
      
      geneset1_overlap <- length(intersect(genesets1[[i]], genesets2[[j]]))
      geneset1_nonoverlap <- length(genesets1[[i]]) - geneset1_overlap
      
      background_overlap <- length(intersect(background_genes, genesets2[[j]]))
      background_nonoverlap <- length(background_genes) - background_overlap
      
      fisher_table <- matrix(c(geneset1_overlap,background_overlap,geneset1_nonoverlap,background_nonoverlap), nrow = 2,
                             dimnames =
                               list(c("Geneset1", "Non-emergent Genes"),
                                    c("Overlap with Raj", "Non-overlap with Raj")))
      
      #Run fisher test
      fisher_table <- fisher_table+1
      fisher_results <- fisher.test(fisher_table, alternative = "two.sided")
      
      
      manual_or <- (fisher_table[1,1]/fisher_table[1,2])/(fisher_table[2,1]/fisher_table[2,2])
      or_mat[i,j] <- manual_or
      
      
      # if(signif && fisher_results$p.value > 0.05){
      #   or_mat[i,j] <- 1
      # }
      
      # or_mat[i,j] <- fisher_results$estimate
      
    }
  }
  
  rownames(or_mat) <- names(genesets1)
  colnames(or_mat) <- names(genesets2)
  
  return(t(or_mat))
}

#Compares two genesests using jaccard overlap 
calc_jaccard_matrix <- function(genesets1, genesets2){
  
  jaccard_matrix <- matrix(NA, nrow=length(genesets1), ncol=length(genesets2))
  
  for(i in 1:length(genesets1)){
    for(j in 1:length(genesets2)){
      
      jaccard_matrix[i,j] <- length(intersect(genesets1[[i]],genesets2[[j]]))/length(union(genesets1[[i]],genesets2[[j]]))
    }
  }
  
  colnames(jaccard_matrix) <- names(genesets2)
  rownames(jaccard_matrix) <- names(genesets1)
  
  return(jaccard_matrix)
}


plot_functional_enrichment_dotplot <- function(genes_to_use, background_genes, plot_title){
  
  go_enrich <- enrichGO(gene = genes_to_use,
                        universe = background_genes,
                        OrgDb = "org.Hs.eg.db",
                        keyType = 'SYMBOL',
                        readable = T,
                        ont = "ALL",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)
  
  
  dotplot(go_enrich)+
    ggtitle(plot_title)
}



find_consensus_geneset <- function(genesets, n){
  
  all_genes <- unique(unlist(genesets))
  
  consensus_geneset <- c()
  
  for(curr_gene in all_genes){
    
    if(sum(sapply(genesets, function(x) curr_gene %in% x)) >=n){
      consensus_geneset <- append(consensus_geneset, curr_gene)
    }
  }
  
  return(consensus_geneset)
}


plot_pretty_heatmap <- function(heatmap_data, title, xlab, ylab, legend_title, col_fun){
  
  
  ht <- Heatmap((heatmap_data), col = col_fun, column_title = plot_title,
                name=legend_title, column_names_rot = 45, column_title_gp = gpar(fontsize = 30))
  
  

  draw(ht, padding = unit(c(2, 20, 2, 2), "mm"))
  
}




#Calcualates intra and inter distances of cells in PC space for given group Ident
calc_group_distances <- function(data, ident.1, ident.2 = NULL, sample_size){
  
  if(is.null(ident.2)){
    
    cat("Calculating intra-distances\n")
    curr_cell_names <- colnames(data)[Idents(data) == ident.1]
    
    if(sample_size > length(curr_cell_names)){
      ident1_cell_names <- sample(curr_cell_names, length(curr_cell_names))
    } else {
      ident1_cell_names <- sample(curr_cell_names, sample_size)
    }
    
    res <- as.matrix(dist(data@reductions$PCA@cell.embeddings[ident1_cell_names, 1:10]))
    
    intra_distances <- as.vector(res)
    
    cat("Calculating inter-distances\n")
    
    other_cell_names <-  colnames(data)[Idents(data) != ident.1]
    
    if(sample_size > length(curr_cell_names)){
      other_cell_names <- sample(other_cell_names,  length(curr_cell_names)/2)
      ident1_cell_names <- sample(ident1_cell_names,  length(curr_cell_names)/2)
    } else {
      other_cell_names <- sample(other_cell_names, sample_size/2)
      ident1_cell_names <- sample(ident1_cell_names, sample_size/2)
    }
    
    
    cells_to_use <- c(ident1_cell_names,other_cell_names)
    
    
    res <- as.matrix(dist(data@reductions$PCA@cell.embeddings[cells_to_use,1:10]))
    
    inter_distances <- res[ident1_cell_names,other_cell_names]
    
    inter_distances <- as.vector(inter_distances)
    
    wilcox_result <- wilcox.test(intra_distances, inter_distances)
    
    distances_list <- list("intra_distances" = intra_distances, "inter_distances" = inter_distances)
    
  } else{
    
  }
  
  return(list("distances_list" = distances_list, "wilcox_result" = wilcox_result))
  
}


plot_km_plot <- function(fit, data_to_use, plot_title,xlab,ylab){
  
  p <- ggsurvplot(fit, data = data_to_use,
                  pval = T,
                  pval.size = 10,
                  pval.coord = c(0, 0.03),
                  font.title=c(28),
                  font.x=c(28),
                  font.y=c(28),
                  font.tickslab=c(20),
                  font.legned=c(50),
                  font.caption=c(30),
                  tables.theme = clean_theme(),
                  legend.title="",
                  legend.labs= c("High","Low"),)+
    xlab(xlab)+
    ylab(ylab)
  
  
  p$plot <- p$plot + 
    theme(legend.text = element_text(size = 24, color = "black", face = "bold"),
          legend.title = element_text(size = 28, color = "black", face = "bold"),
          legend.key.height = unit(1.5,"cm"),
          legend.key.width = unit(1.5,"cm"),
          legend.position = c(0.1, 0.25))
  
  
  p <- p + ggtitle(plot_title)
  
  
  return(p)
}

