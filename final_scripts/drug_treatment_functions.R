library(Seurat)
library(fgsea)
library(ComplexHeatmap)
library(circlize)

create_GSEA_matrix <- function(genesets_with_ranks, genesets2){
  
  #Remove genes from genesets2 that are not in the ranked gene list
  for(i in 1:length(genesets2)){
    genesets2[[i]] <- genesets2[[i]][genesets2[[i]] %in% names(genesets_with_ranks[[1]])]
  }
  
  # Create matrix of NES. rows = ranked geneset, cols = genesets of interest
  gsea_mat <- matrix(NA, nrow=length(genesets_with_ranks), ncol=length(genesets2))
  rownames(gsea_mat) <- names(genesets_with_ranks)
  colnames(gsea_mat) <- names(genesets2)

  # Run GSEA for each ranked geneset  
  gsea_results <- list()
  for(i in 1:length(genesets_with_ranks)){
    
    fgseaRes <- fgsea(pathways = genesets2, 
                      stats    = genesets_with_ranks[[i]],
                      minSize  = 10,
                      maxSize  = max(lengths(genesets2))+1)
    
    gsea_results <- append(gsea_results, list(fgseaRes))
    
    
    gsea_mat[i,(colnames(gsea_mat) %in% fgseaRes$pathway)] <- fgseaRes$NES
    gsea_mat[i,!(colnames(gsea_mat) %in% fgseaRes$pathway)] <- 0
  }
  
  names(gsea_results) <- sapply(names(genesets_with_ranks), FUN = function(x){gsub("_fc_results", "", x)})
  
  gsea_mat[is.na(gsea_mat)] <- 0
  gsea_mat <- t(gsea_mat)
  
  #Return GSEA NES matrix and all GSEA results objects
  return(list(gsea_mat,gsea_results))
}

watermelon_validation <- function(curr_signature, plot_title=NULL){
  
  if(!exists("watermelon_data")){
    watermelon_data <- readRDS(paste0(dataDirectory, "processed_data/watermelon_data/watermelon_pc9_processed.rds"))
  }
  
  ################################################################################
  # Calculate FC for each time point
  ################################################################################
  
  Idents(watermelon_data) <- watermelon_data$time_point
  time_points <- levels(Idents(watermelon_data))
  time_points_fc <- list()
  for(curr_time_point in time_points){
    
    cat(curr_time_point,"\n")
    curr_fc_result <- FoldChange(watermelon_data, ident.1 = curr_time_point)
    
    time_points_fc <- append(time_points_fc, list(curr_fc_result))
  }
  
  names(time_points_fc) <- paste0("pc9_day",time_points,"_fc")
  
  ################################################################################
  # Run GSEA of current signature for each time point
  ################################################################################
  
  #Create ranked genelists for each time point
  genesets_with_ranks <- list()
  for(curr_fc_res in time_points_fc){
    
    ranks <- curr_fc_res$avg_log2FC
    
    names(ranks) <- rownames(curr_fc_res)
    
    ranks <- sort(ranks, decreasing = T)
    
    genesets_with_ranks <- append(genesets_with_ranks, list(ranks))
    
  }
  
  names(genesets_with_ranks) <- paste0("pc9_day",time_points,"_fc")
  
  
  # Run GSEA
  result <- create_GSEA_matrix(genesets_with_ranks, curr_signature)
  
  gsea_matrix <- result[[1]]
  
  colnames(gsea_matrix) <- paste0("Day ", time_points)
  
  if(nrow(gsea_matrix) == 1){
    rownames(gsea_matrix) <- ""  
  }
  
  
  return(gsea_matrix)
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

geneset_group_matrix <- function(data, ident_to_use, genesets_to_use){
  
  Idents(data) <- ident_to_use
  
  num_groups <- nlevels(data)
  
  groups <- levels(data)
  
  cat(genesets_to_use, "\n")
  
  scores <- readRDS(paste0(dataDirectory, "data/aucell_score_objects/watermelon_pc9_processed_", genesets_to_use, "_aucell_scores.rds"))
  threshold <- readRDS(paste0(dataDirectory, "data/aucell_score_objects/watermelon_pc9_processed_", genesets_to_use, "_aucell_thresholds.rds"))
  
  # scores <- scale(scores)
  
  geneset_names <- colnames(scores)
  
  score_matrix <- matrix(NA, ncol=length(groups), nrow=length(geneset_names))
  pvalue_matrix <- matrix(NA, ncol=length(groups), nrow=length(geneset_names))
  
  for(curr_group in groups){
    
    cat(curr_group, "\n")
    
    curr_cell_names <- colnames(data)[data[[ident_to_use]] == curr_group]
    rest_names <- colnames(data)[data[[ident_to_use]] != curr_group]
    
    for(curr_geneset in geneset_names){
      
      curr_scores <- scores[curr_cell_names, curr_geneset]
      
      rest_scores <- scores[rest_names, curr_geneset]
      
      wilcox_res <- wilcox.test(curr_scores,rest_scores)
      
      i <- which(curr_geneset == geneset_names)
      j <- which(curr_group == groups)
      
      score_matrix[i,j] <- mean(curr_scores)
      pvalue_matrix[i,j] <- wilcox_res$p.value
      
    }
  }
  
  scaled_matrix <- t(scale(t(score_matrix)))
  
  colnames(scaled_matrix) <- groups
  rownames(scaled_matrix) <- geneset_names
  
  colnames(score_matrix) <- groups
  rownames(score_matrix) <- geneset_names
  
  return(list(scaled_matrix, pvalue_matrix,score_matrix))
}

genelist_to_table <- function(curr_list){
  
  max_len <- max(lengths(curr_list))
  
  curr_list <- lapply(curr_list, sort)
  
  curr_list <- lapply(curr_list, FUN = function(x) append(x,rep("",max_len-length(x))))
  
  curr_list <- data.frame(curr_list)
  
  return(curr_list)
}

create_barplot <- function(result, y_label, plot_title){
  df <- result@result %>% 
    filter(p.adjust < 0.05) %>% 
    arrange(desc(NES)) %>% 
    dplyr::select(Description,NES,p.adjust)
  
  
  plot <- ggplot(df[1:5,])+
    geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES),width=0.7)+
    scale_fill_continuous(high="firebrick1",low="firebrick1")+
    xlab(y_label)+
    ylab("NES")+
    ggtitle(paste0(plot_title))+
    theme(plot.title = element_text(size=30),
          axis.text = element_text(size=35),
          axis.title = element_text(size=35,face="bold"),
          legend.key.size=unit(1,'cm'),
          legend.title = element_text(size=20),
          legend.text = element_text(size=15))+
    NoLegend()+
    coord_flip()
  
  return(plot)
}

genesets_characterization <- function(genesets_to_use, universe_to_use, num_pathways=10){
  
  m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, human_gene_symbol)
  
  mp_t2g <- readRDS(paste0(dataDirectory, "genesets/ith_meta_programs_t2g.rds"))
  
  specifc_mps <- c("MP39 Metal-response","MP31 Alveolar","MP29 NPC/OPC","MP28 Oligo normal","MP27 Oligo Progenitor","MP38 Glutathione","MP41 Unassigned","MP35 Hemato-related-I","MP37 Hemato-related-II","MP32 Skin-pigmentation","MP36 IG","MP16 MES (glioma)","MP15 EMT IV")
  
  mp_t2g <- mp_t2g %>% 
    filter(!gs_name %in% specifc_mps)
  
  
  hallmarks_dotplots <- list()
  mps_dotplots <- list()
  go_dotplots <- list()
  go_results <- list()
  
  for(i in 1:length(genesets_to_use)){
    
    curr_cluster <- names(genesets_to_use)[i]
    cat(curr_cluster, "\n")
    
    curr_geneset <- genesets_to_use[[i]]
    
    
    # Hallmarks Enrichment
    hallmark_enrichment_results <- enricher(curr_geneset, TERM2GENE=m_t2g,
                                            universe = universe_to_use)
    
    # MPs Enrichment
    mp_enrichment_results <- enricher(curr_geneset, TERM2GENE=mp_t2g,
                                      universe = universe_to_use)
    
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
    
    go_results <- append(go_results,list(ego))
    
    
    if(nrow(hallmark_enrichment_results) > 0){
      p <- barplot(hallmark_enrichment_results,
                   showCategory = num_pathways, font.size=20) + 
        ggtitle(paste0("", curr_cluster))+
        theme(plot.title = element_text(size=20))
      
      hallmarks_dotplots <- append(hallmarks_dotplots, list(p))
    }
    if(!is.null(mp_enrichment_results) && nrow(mp_enrichment_results) > 0){
      p <- barplot(mp_enrichment_results,
                   showCategory = num_pathways, font.size=20)+
        ggtitle(paste0("", curr_cluster))+
        theme(plot.title = element_text(size=20))
      
      mps_dotplots <- append(mps_dotplots, list(p))
    }
    if(nrow(ego) > 0){
      p <- barplot(ego,
                   showCategory = num_pathways, font.size=20) + 
        ggtitle(paste0("", curr_cluster))+
        theme(plot.title = element_text(size=20))
      
      go_dotplots <- append(go_dotplots, list(p))
    }
  }
  
  
  ret_list <- list(hallmarks_dotplots, mps_dotplots, go_dotplots,go_results)
  names(ret_list) <- c("hallmarks_dotplots", "mps_dotplots", "go_dotplots","go_results")
  
  return(ret_list) 
}

get_tcga_project <- function(cell_line){
  
  if(cell_line == "A549"){
    tcga_project <- "TCGA-LUAD"
  } else if(cell_line == "K562"){
    tcga_project <- "TCGA-LAML"
  } else{
    tcga_project <- "TCGA-BRCA"
  }
  
  return(tcga_project)
}

get_read_count_data <- function(i){
  ################################################################################
  # Get read count data
  ################################################################################
  
  query_tcga <- GDCquery(project = i,
                         data.category = 'Transcriptome Profiling',
                         data.type = 'Gene Expression Quantification',
                         experimental.strategy = "RNA-Seq",
                         workflow.type = 'STAR - Counts',
                         access = 'open')
  
  output_query_tcga <- getResults(query_tcga)
  
  
  all_samples <- output_query_tcga[grepl("Primary",output_query_tcga$sample_type), "cases"]
  
  query_tcga <- GDCquery(project = i,
                         data.category = 'Transcriptome Profiling',
                         data.type = 'Gene Expression Quantification',
                         experimental.strategy = "RNA-Seq",
                         workflow.type = 'STAR - Counts',
                         access = 'open',
                         barcode = all_samples)
  
  
  output_query_tcga <- getResults(query_tcga)
  
  # Download data - GDCdownload
  GDCdownload(query_tcga)
  
  # Prepare data
  tcga_data <- GDCprepare(query_tcga, summarizedExperiment = T)
  
  tcga_matrix <- assay(tcga_data, 'unstranded')
  
  
  # extract gene and sample metadata from summarizedExperiment object
  gene_metadata <- as.data.frame(rowData(tcga_data))
  coldata <- as.data.frame(colData(tcga_data))
  
  # Normalize data with DESeq
  dds <- DESeqDataSetFromMatrix(countData = tcga_matrix,
                                colData = coldata,
                                design = ~ 1)
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  # variance stabalizing transformation
  vsd <- vst(dds, blind = F)
  tcga_matrix_vst <- assay(vsd)
  
  
  temp <- gene_metadata %>% 
    filter(gene_id %in% rownames(tcga_matrix_vst)) %>% 
    select(gene_id,gene_name)
  
  
  if(all(temp$gene_id == rownames(tcga_matrix_vst))){
    rownames(tcga_matrix_vst) <- temp$gene_name
  }
  
  return(tcga_matrix_vst)
}

get_clinical_data <- function(project){
  
  clinical_data <- read_xlsx("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/TCGA_clinical_data.xlsx")
  colnames(clinical_data)[2] <- "submitter_id"
  
  clinical_data$age <- clinical_data$age_at_initial_pathologic_diagnosis
  clinical_data$sex <- clinical_data$gender
  
  tcga_type <- strsplit(project, "-")[[1]][2]
  if(tcga_type %in% Tumor.purity$Cancer.type){
    
    purity_data <- Tumor.purity
    
    colnames(purity_data)[1] <- "submitter_id"
    
    purity_data$submitter_id <- sapply(purity_data$submitter_id, function(x) substring(x, 1,12))
    
    clinical_data <- merge(clinical_data,purity_data,by = "submitter_id")
    
    clinical_data$CPE <- as.numeric(sapply(clinical_data$CPE, function(x) gsub(",",".",x)))
    
    clinical_data$purity <- clinical_data$CPE
  }
  
  return(clinical_data)
}


#results = list of enrichment results that will be the columns in the heatmap
create_enrichment_heatmap <- function(results, title){
  
  heatmap_df <- list()
  for(curr_cluster in names(results)){
    
    curr_result <- results[[curr_cluster]]
    
    if(!is.null(curr_result)){
      
      if("NES" %in% colnames(curr_result@result)){ # is GSEA results
        legend_title <-  "NES"
        
        curr_result <- curr_result@result %>% 
          filter(p.adjust < 0.05) %>% 
          arrange(desc(NES)) %>% 
          dplyr::select(Description,NES)
        
      } else { # is ORA results
        legend_title <-  "Count"
        
        curr_result <- curr_result@result %>% 
          filter(p.adjust < 0.05) %>% 
          arrange(desc(Count)) %>% 
          dplyr::select(Description,Count)
      }
      
      colnames(curr_result) <- c("Description","value")
      
      if(nrow(curr_result) > 5){
        curr_result <- curr_result[1:5,]
      }    
      
      
      heatmap_df[["cluster"]] <- append(heatmap_df[["cluster"]], rep(curr_cluster,nrow(curr_result)))
      heatmap_df[["pathway"]] <- append(heatmap_df[["pathway"]], curr_result$Description)
      heatmap_df[["value"]] <- append(heatmap_df[["value"]], curr_result$value)
    }
    
    
    
  }
  
  
  heatmap_df <- data.frame(heatmap_df)  
  
  
  heatmap_df <- heatmap_df %>% 
    pivot_wider(names_from = "cluster",values_from = "value") %>% 
    column_to_rownames("pathway")
  
  
  min_value <- min(heatmap_df[!is.na(heatmap_df)])
  max_value <- max(heatmap_df[!is.na(heatmap_df)])
  
  if(min_value > 0){
    col_fun = colorRamp2(c(0,max_value), c("white", "red1"))  
  } else if(max_value < 0){
    col_fun = colorRamp2(c(min_value,0), c("royalblue", "white"))
  } else{
    col_fun = colorRamp2(c(min_value,0, max_value), c("royalblue", "white", "red1"))
  }
  
  
  
  curr_title <- paste0(title)
  
  curr_ht <- Heatmap(as.matrix(heatmap_df),cluster_rows = F,cluster_columns = F, col = col_fun,
                     name=legend_title, column_names_rot = 45, column_title = curr_title,
                     column_names_gp = gpar(fontsize=15),row_names_gp = gpar(fontsize=20),
                     column_title_gp = gpar(fontsize=22),
                     row_names_max_width = max_text_width(rownames(heatmap_df)),
                     row_names_side = "left",
                     heatmap_legend_param = list(title_gp = gpar(fontsize = 22),legend_height = unit(3, "cm"), grid_width=unit(1,"cm"),
                                                 labels_gp = gpar(fontsize = 14)))
  
  # curr_ht <- draw(curr_ht, heatmap_legend_side  = "left")
  
  return(curr_ht)
}


create_consensus_ranks <- function(ranks_list){
  
  all_genes <- sort(unique(names(unlist(ranks_list))))
  
  avg_ranks <- c()
  for(curr_gene in all_genes){
    
    has_gene <- unlist(lapply(ranks_list, FUN = function(x) curr_gene %in% names(x)))
    
    curr_gene_ranks <- c()
    for(i in ranks_list[has_gene]){
      
      curr_gene_ranks <- append(curr_gene_ranks, i[which(curr_gene == names(i))])
      
    }
    
    avg_ranks <- append(avg_ranks, mean(curr_gene_ranks))
    
  }
  
  names(avg_ranks) <- all_genes
  
  # Create consensus ranks
  
  gene_list <- list()
  for(i in ranks_list){
    
    gene_list <- append(gene_list, list(names(i)))
    
  }
  
  shared_genes <- find_consensus_geneset(gene_list,3)
  
  
  consensus_ranks <- sort(avg_ranks[names(avg_ranks) %in% shared_genes])
  
  return(consensus_ranks)
}
