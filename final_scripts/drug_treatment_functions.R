library(Seurat)
library(tidyverse)
library(fgsea)
library(monocle3)
library(AUCell)
library(GSEABase)
library(data.table)
library(foreach)
library(doMC)
library(ggpubr)
library(ComplexHeatmap)
library(matrixStats)
library(circlize)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(cowplot)
library(presto)
library(decoupleR)
library(OmnipathR)
library(GSVA)
library(readxl)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(survminer)
library(survival)
library(patchwork)
library(latex2exp)
library(grid)
library(scales)
library(org.Hs.eg.db)
library(org.Sc.sgd.db)
library(org.EcK12.eg.db)
library(rstatix)
library(AnnotationDbi)
library(lmtest)
library(xlsx)
library(poolr) 
library(readxl)
library(decoupleR)
library(rstatix) 
library(RobustRankAggreg)
library(effsize)
library(factoextra)
library(NbClust)

set.seed(42)

plot_enrichment_ora <- function(genesets){
  # Hallmarks term2gene list
  m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, gene_symbol)
  
  new_geneset_names <- sapply(m_t2g$gs_name, FUN = function(x) gsub("HALLMARK_", "", x))
  new_geneset_names <- sapply(new_geneset_names, FUN = function(x) gsub("_", " ", x))
  m_t2g$gs_name <- new_geneset_names
  
  # ITH Meta-programs term2gene list
  mp_t2g <- readRDS(paste0(dataDirectory, "genesets/ith_meta_programs_t2g.rds"))
  
  #Remove 'specific' MPs
  specifc_mps <- c("MP39 Metal-response","MP31 Alveolar","MP29 NPC/OPC","MP28 Oligo normal","MP27 Oligo Progenitor","MP38 Glutathione","MP41 Unassigned","MP35 Hemato-related-I","MP37 Hemato-related-II","MP32 Skin-pigmentation","MP36 IG","MP16 MES (glioma)","MP15 EMT IV")
  mp_t2g <- mp_t2g %>% 
    filter(!gs_name %in% specifc_mps)
  
  ######
  
  titles <- list("Cancer Hallmarks","ITH Meta-programs","GO Pathways")
  
  
  all_plots <- list()
  for(j in 1:length(genesets)){
    
    curr_plots <- list()
    for(i in 1:3){
      
      curr_geneset <- genesets[[j]]
      
      if(i == 1){
        # Hallmarks Enrichment
        curr_enrichment_results <- enricher(curr_geneset, TERM2GENE=m_t2g)
        
        
      } else if(i == 2){
        # MPs Enrichment
        curr_enrichment_results <- enricher(curr_geneset, TERM2GENE=mp_t2g)
        
        
      } else{
        # GO Functional Enrichment
        curr_enrichment_results <- enrichGO(gene          = curr_geneset,
                                            OrgDb         = org.Hs.eg.db,
                                            ont           = "BP",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff  = .01,
                                            qvalueCutoff  = .05,
                                            readable      = TRUE,
                                            keyType = "SYMBOL")
        
      }
      
      if(!(is.null(curr_enrichment_results))){
        if(!nrow(curr_enrichment_results) == 0){
          if(nrow(curr_enrichment_results) == 1){
            df <- as.data.frame(curr_enrichment_results)
            # Convert enrichment result to data frame
            
            # Keep only necessary columns for plotting
            df_plot <- df[, c("Description", "Count","p.adjust")]
            
            # Add dummy rows with matching column names
            dummy_top <- data.frame(Description = "dummy1", Count = 0,p.adjust=df_plot$p.adjust[1]+.00000001)
            dummy_bottom <- data.frame(Description = "dummy2", Count = 0,p.adjust=df_plot$p.adjust[1]+.0000001)
            
            # Combine them with the real data
            df_full <- rbind(dummy_top, df_plot, dummy_bottom)
            
            # Set Description as a factor to preserve order
            df_full$Description <- factor(df_full$Description, levels = df_full$Description)
            
            df_full$p.adjust <- as.numeric(df_full$p.adjust)
            
            
            # Plot
            p <- ggplot(df_full, aes(x = Description, y = Count,fill=p.adjust)) +
              geom_bar(stat = "identity", width = 0.5, na.rm = TRUE) +
              coord_flip() +
              ggtitle(paste0(titles[i]))+
              theme_classic()+
              scale_fill_gradient(low="firebrick4", high="royalblue4")+
              theme(axis.text.y = element_text(size = 8 , color = ifelse(df_full$Description %in% c("dummy1", "dummy2"), "transparent", "black")),
                    plot.title = element_text(size=10),
                    axis.text.x = element_text(size = 8),
                    legend.text = element_text(size=4),
                    legend.title = element_text(size=6),
                    legend.key.size = unit(2, "mm"),
                    axis.text = element_text(size=8),
                    axis.title = element_text(size=8),
                    axis.line = element_line(linewidth=.2),
                    axis.ticks = element_line(linewidth = .2))+ 
              scale_size_area(max_size = 3)+
              labs(x="")
            
            
            
          } else{
            
            df <- as.data.frame(curr_enrichment_results)
            # Convert enrichment result to data frame
            
            # Keep only necessary columns for plotting
            df_plot <- df[, c("Description", "Count","p.adjust")]
            
            df_plot <- df_plot %>% 
              arrange((p.adjust)) %>% 
              head(5)
            
            p <- ggplot(df_plot, aes(x = reorder(Description,desc(p.adjust)), y = Count,fill=p.adjust)) +
              geom_bar(stat = "identity", width = 0.5, na.rm = TRUE) +
              coord_flip() +
              ggtitle(paste0(titles[i]))+
              theme_classic()+
              scale_fill_gradient(low="firebrick4", high="royalblue4")+
              theme(axis.text.y = element_text(size=8,color = "black"),
                    plot.title = element_text(size=10),
                    axis.text.x = element_text(size = 8),
                    legend.text = element_text(size=4),
                    legend.title = element_text(size=6),
                    legend.key.size = unit(2, "mm"),
                    axis.text = element_text(size=8),
                    axis.title = element_text(size=8),
                    axis.line = element_line(linewidth=.2),
                    axis.ticks = element_line(linewidth = .2))+ 
              scale_size_area(max_size = 3)+
              labs(x="")
            
          }
          
          curr_plots <- append(curr_plots, list(p))
          
        } 
      }
    }
    
    if(length(curr_plots) > 0){
      # plot <- ggarrange(plotlist = curr_plots, nrow=1, common.legend = T, legend = "right",
      #                   widths = c(1,1,1),heights = c(1,1,1))
      
      
      if(length(curr_plots) == 3){
        plot <- (curr_plots[[1]]|curr_plots[[2]]|curr_plots[[3]]) + plot_annotation(names(genesets)[j],theme=theme(plot.title=element_text(hjust=0.5,size=12)))
      } else {
        plot <- (curr_plots[[1]]|curr_plots[[2]]) + plot_annotation(names(genesets)[j],theme=theme(plot.title=element_text(hjust=0.5,size=12)))
      }
      
      
      
      
      # plot <- annotate_figure(plot, top = text_grob(names(genesets)[j], size = 12))
      
      all_plots <- append(all_plots, list(plot))
    }
    
  }
  
  
  
  p <- ggarrange(plotlist = all_plots, nrow=length(all_plots))
  
  
  return(p)
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
          arrange(desc(FoldEnrichment)) %>% 
          dplyr::select(Description,Count)
      }
      
      colnames(curr_result) <- c("Description","value")
      
      if(nrow(curr_result) > 5){
        curr_result <- curr_result[1:5,]
      } 
      
      
      heatmap_df[["cluster"]] <- append(heatmap_df[["cluster"]], rep(curr_cluster,nrow(curr_result)))
      
      curr_result$Description <- sapply(curr_result$Description, FUN = function(x) paste(strwrap(x, 20),collapse="\n"))
      heatmap_df[["pathway"]] <- append(heatmap_df[["pathway"]], unlist(curr_result$Description))
      
      heatmap_df[["value"]] <- append(heatmap_df[["value"]], curr_result$value)
    }
    
    
    
  }
  
  heatmap_df <- data.frame(heatmap_df)  
  
  
  heatmap_df <- heatmap_df %>% 
    pivot_wider(names_from = "cluster",values_from = "value") %>% 
    column_to_rownames("pathway")
  
  
  min_value <- min(heatmap_df[!is.na(heatmap_df)])
  max_value <- max(heatmap_df[!is.na(heatmap_df)])
  min_max <- c(min_value,max_value)
  
  
  curr_title <- paste0(title)
  
  if(nrow(as.matrix(heatmap_df)) == 1){
    curr_ht <- Heatmap(as.matrix(heatmap_df),cluster_rows = F,cluster_columns = F, col = "red",rect_gp = gpar(col = "white", lwd = 2),
                       name=legend_title, column_names_rot = 45, column_title = curr_title,
                       column_names_gp = gpar(fontsize=10), row_names_gp = gpar(fontsize=10),
                       column_title_gp = gpar(fontsize=16),
                       row_names_max_width = max_text_width(rownames(heatmap_df)),
                       row_names_side = "left", show_heatmap_legend = FALSE,
                       heatmap_legend_param = list(title_gp = gpar(fontsize = 10),legend_height = unit(5, "mm"), grid_width=unit(1.5,"mm"),
                                                   labels_gp = gpar(fontsize = 8)))
    
  } else {
    curr_ht <- Heatmap(as.matrix(heatmap_df),cluster_rows = F,cluster_columns = F,rect_gp = gpar(col = "white", lwd = 2),
                       name=legend_title, column_names_rot = 45, column_title = curr_title,
                       column_names_gp = gpar(fontsize=10), row_names_gp = gpar(fontsize=10),
                       column_title_gp = gpar(fontsize=16),
                       row_names_max_width = max_text_width(rownames(heatmap_df)),
                       row_names_side = "left", show_heatmap_legend = FALSE,
                       heatmap_legend_param = list(title_gp = gpar(fontsize = 10),legend_height = unit(5, "mm"), grid_width=unit(1.5,"mm"),
                                                   labels_gp = gpar(fontsize = 8)))
  }
  
  return(list(curr_ht, min_max))
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

genelist_to_table <- function(curr_list){
  
  max_len <- max(lengths(curr_list))
  
  curr_list <- lapply(curr_list, sort)
  
  curr_list <- lapply(curr_list, FUN = function(x) append(x,rep("",max_len-length(x))))
  
  curr_list <- data.frame(curr_list)
  
  return(curr_list)
}
