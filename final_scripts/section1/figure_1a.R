args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/revision_figures/"
################################################################################

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
                      maxSize  = max(lengths(genesets2))+1,
                      nPermSimple = 10000)
    
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
  
  # If watermelon timepoints genesets with ranks already exists, just read it in
  # Else create it and save it
  if (!file.exists(paste0(dataDirectory, "processed_data/oren_data/oren_time_points_genesets_with_ranks.rds"))) {
    if(!exists("watermelon_data")){
      watermelon_data <- readRDS(paste0(dataDirectory, "processed_data/oren_data/oren_pc9_processed.rds"))
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
    saveRDS(genesets_with_ranks, paste0(dataDirectory, "processed_data/oren_data/oren_time_points_genesets_with_ranks.rds"))
  } else {
    genesets_with_ranks <- readRDS(paste0(dataDirectory, "processed_data/oren_data/oren_time_points_genesets_with_ranks.rds"))
    time_points <- c(0,3,7,14)
  }
  
  
  # Run GSEA
  result <- create_GSEA_matrix(genesets_with_ranks, curr_signature)
  
  gsea_results <- do.call(rbind,result[[2]])
  
  gsea_results[["day"]] <- time_points
  
  return(gsea_results)
}

################################################################################

resistance_signature <- readRDS(paste0(dataDirectory, "genesets/resistance_signature.rds"))

# Calculate GSEA enrichment at each time point for resistance signature
gsea_results <- watermelon_validation(resistance_signature)

# Plot NES of resistance signature at each time point as line plot
df <- as.data.frame(gsea_results)
# colnames(df) <- c("day","value")



# df$value <- as.numeric(df$value)
df$day <- factor(df$day, levels = df$day)


# df <- df %>% 
#   filter(day %in% c(0,14))

df$dot_color <- ifelse(df$NES > 0, "pos","neg")

p <- ggplot(df, aes(x=day, y=NES, fill=dot_color)) + 
  geom_bar(stat = "identity",color="black",linewidth=.2)+
  geom_text(aes(x=day,label=paste0("p-val: ", sprintf(fmt = "%.5f", padj))), y= 2.2, size=2)+
  scale_fill_manual(values=c("pos" = "red","neg" ="blue"))+
  xlab("Time Point")+
  ylab("NES")+
  geom_hline(yintercept = 0,color="black", linewidth = .2)+
  theme_classic()+
  theme(axis.title = element_text(size=8, color = "black"),
        axis.text = element_text(size=6, color = "black"),
        title = element_text(size=8),
        axis.line = element_line(linewidth=.2),
        axis.ticks = element_line(linewidth = .2))+
  NoLegend()+
  geom_errorbar(aes(ymin=NES-log2err, ymax=NES+log2err), size=.2,width=.1,
                position=position_dodge(.9)) +
  scale_y_continuous(name="NES", limits=c(-3.1, 3))

p

jpeg(paste0(plotDirectory,"figure_1a.jpg"),width=80, height = 60, units = "mm", res = 1000)
print(p)
dev.off()



