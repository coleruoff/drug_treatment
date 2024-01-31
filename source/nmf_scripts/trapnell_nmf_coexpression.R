library(Seurat)
library(NNLM)
library(cluster)
library(ComplexHeatmap)
library(tidyverse)
library(ggpubr)
library(UCell)

A549.data <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/A549_processed_filtered.rds")

K562.data <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/K562_processed_filtered.rds")

MCF7.data <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/MCF7_processed_filtered.rds")

subset_cell_names <- c()

subset_percent <- .4

for(i in 1:nlevels(A549.data)){
    subset_cell_names <- append(subset_cell_names,sample(colnames(A549.data)[A549.data$Cluster == i], size = subset_percent*sum(A549.data$Cluster == i)))
}

A549.data <- A549.data[,subset_cell_names]

subset_cell_names <- c()

for(i in 1:nlevels(K562.data)){
    subset_cell_names <- append(subset_cell_names,sample(colnames(K562.data)[K562.data$Cluster == i], size = subset_percent*sum(K562.data$Cluster == i)))
}

K562.data <- K562.data[,subset_cell_names]

subset_cell_names <- c()

for(i in 1:nlevels(MCF7.data)){
    subset_cell_names <- append(subset_cell_names,sample(colnames(MCF7.data)[MCF7.data$Cluster == i], size = subset_percent*sum(MCF7.data$Cluster == i)))
}

MCF7.data <- MCF7.data[,subset_cell_names]

cat("test1")

# tumor_obj <- readRDS("/data/CDSL_hannenhalli/Cole/BrM2_10x_cluster_object.pc20_res0.2_BR13.rds")

# dim(tumor_obj)

trapnell.combined  <- merge(A549.data, y = c(K562.data, MCF7.data), add.cell.ids = c("A549", "K562", "MCF7"), project = "trapnell.combined")

trapnell.combined <- RenameAssays(trapnell.combined, RNA = "RNA")

dim(trapnell.combined)
rm(A549.data)
rm(K562.data)
rm(MCF7.data)

# saveRDS(trapnell.combined, "/data/CDSL_hannenhalli/Cole/trapnell_combined.rds")

# tumor_obj <- readRDS("/data/CDSL_hannenhalli/Cole/BrM2_10x_cluster_object.pc20_res0.2_BR13.rds")

# tumor_obj <- readRDS("/data/CDSL_hannenhalli/Cole/trapnell_combined.rds")

tumor_obj <- trapnell.combined
rm(trapnell.combined) 

tumor_obj$sample_type <- tumor_obj$cell_type
DefaultAssay(tumor_obj) <- "RNA"
tumor_obj <- NormalizeData(tumor_obj) 


cat("test2")
gene_exp_mat <- as.matrix(tumor_obj[["RNA"]]@data)
gene_exp_frac <- rowMeans(gene_exp_mat > 0)

# gene_exp_mat <- as.matrix(tumor_obj[["RNA"]]@data)
# gene_exp_frac <- rowMeans(gene_exp_mat > 0)

frac_exp_threshold <- 0.05 #Genes expressed in less than this fraction of cells will not be considered for NMF
genes_to_keep <- names(which(gene_exp_frac >= frac_exp_threshold))
gene_exp_mat <- gene_exp_mat[genes_to_keep,]

dim(gene_exp_mat)

num_factors <- 30   
nmf_list <- list()

#In the loop below, we carry out separate NMFs for each tumor model (cell line). The resulting factor matrices are stored in a list of 
# lists named after each cell line
for (sample_type_ in c("A549","K562","MCF7")) {
  cells <- WhichCells(tumor_obj, expression = sample_type == sample_type_ )
  sample_type_gene_exp_mat <- gene_exp_mat[,cells]
    
  print(num_factors)
  flush.console()
    
  nmf_out <- NNLM::nnmf(sample_type_gene_exp_mat,k=num_factors,n.threads=6)
  factor_loading_mat <- nmf_out$W
  factor_score_mat <- nmf_out$H
    
  colnames(factor_loading_mat) <- paste(sample_type_,"Factor",1:num_factors,sep="_")
  rownames(factor_score_mat) <- paste(sample_type_,"Factor",1:num_factors,sep="_")
  
  nmf_list[[sample_type_]][["Loadings"]] <- factor_loading_mat
  nmf_list[[sample_type_]][["Score"]] <- factor_score_mat
  
}

saveRDS(nmf_list, "/data/CDSL_hannenhalli/Cole/nmf_files/nmf_list.rds")

# nmf_list <- readRDS("/data/CDSL_hannenhalli/Cole/nmf_files/nmf_list.rds")

# rm(trapnell.combined)
ls()


#################Factor co-expression test
z_score_thresh <- 1.5              #Z-score threshold above which a gene will be assigned to that factor
min_factor_size <- 10            #Minimum number of genes in each factor. Factors with fewer genes than this threshold will be discarded
factor_gene_list <- list()          #List of factors, with each factor being a vector of genes.
cor_threshold <- 0                #Threshold above which two genes are considered to be co-expressed.

coexp_factor_gene_list <- list()    #List of factors, with each factor being a vector of genes where each gene is co-expressed with at least one other gene in the factor
A549_cells <- WhichCells(tumor_obj,expression = sample_type == "A549")
K562_cells <- WhichCells(tumor_obj,expression = sample_type == "K562")
MCF7_cells <- WhichCells(tumor_obj,expression = sample_type == "MCF7")

A549_exp_mat <- gene_exp_mat[,A549_cells]
K562_exp_mat <- gene_exp_mat[,K562_cells]
MCF7_exp_mat <- gene_exp_mat[,MCF7_cells]

for (sample_type_ in c("A549","K562","MCF7")) {
  factor_loading_mat <- nmf_list[[sample_type_]]$Loadings
  z_scored_factor_loading_mat <- t(scale(t(factor_loading_mat)))
  factor_names <- colnames(factor_loading_mat)
  
  factor_gene_list[[sample_type_]] <- list()
  coexp_factor_gene_list[[sample_type_]] <- list()
  
  for (factor in factor_names) {
    print(factor)
    flush.console()
      
    #Assign genes with the highest z-scores to the factor
    factor_genes <- names(which(z_scored_factor_loading_mat[,factor] > z_score_thresh))#[1:num_top_genes] %>% na.omit
    factor_gene_list[[sample_type_]][[factor]] <- factor_genes
      
    factor_genes <- factor_genes[factor_genes %in% rownames(A549_exp_mat)]  
    factor_genes <- factor_genes[factor_genes %in% rownames(K562_exp_mat)]
    factor_genes <- factor_genes[factor_genes %in% rownames(MCF7_exp_mat)]
    
    A549_co_exp_mat <- cor(t(A549_exp_mat[factor_genes,]),method="spearman")
    K562_co_exp_mat <- cor(t(K562_exp_mat[factor_genes,]),method="spearman")
    MCF7_co_exp_mat <- cor(t(MCF7_exp_mat[factor_genes,]),method="spearman")
    
    #Additionally, create another gene list involving factors where gene pairs are highly co-expressed
    A549_co_exp_df <- reshape2::melt(A549_co_exp_mat,value.name="cor_val") %>% 
      na.omit %>% dplyr::filter(Var1 != Var2) %>% 
      dplyr::filter(cor_val >= cor_threshold)
      
    K562_co_exp_df <- reshape2::melt(K562_co_exp_mat,value.name="cor_val") %>% 
      na.omit %>% dplyr::filter(Var1 != Var2) %>% 
      dplyr::filter(cor_val >= cor_threshold)
      
    MCF7_co_exp_df <- reshape2::melt(MCF7_co_exp_mat,value.name="cor_val") %>% 
      na.omit %>% dplyr::filter(Var1 != Var2) %>% 
      dplyr::filter(cor_val >= cor_threshold)
      
    co_exp_df <- merge(A549_co_exp_df,K562_co_exp_df,by=c("Var1","Var2"))
    co_exp_df <- merge(co_exp_df,MCF7_co_exp_df,by=c("Var1","Var2"))
      
    genes_to_retain <- unique(as.character(co_exp_df$Var1))
    
    
    if (length(genes_to_retain) < min_factor_size){
        cat('TOO FEW GENES\n')
        next
    }
      
    
    coexp_factor_gene_list[[sample_type_]][[factor]] <- genes_to_retain
  }
}

combined_factor_list <- c(coexp_factor_gene_list[["A549"]],coexp_factor_gene_list[["K562"]],coexp_factor_gene_list[["MCF7"]])


saveRDS(combined_factor_list, "/data/CDSL_hannenhalli/Cole/nmf_files/combined_factor_list.rds")

#####Score each tumor cell for expression of each factor. Change "ncores" argument based on # of cores available to R.

tumor_ucell_obj <- AddModuleScore_UCell(tumor_obj, features = combined_factor_list, ncores = 4,maxRank = 4500)
for (factor in names(combined_factor_list)) {
  col_name <- grep(paste0(factor,"_"),colnames(tumor_ucell_obj@meta.data),value=T)
  tumor_ucell_obj[[factor]] <- tumor_ucell_obj[[col_name]]
  tumor_ucell_obj[[col_name]] <- NULL
}

saveRDS(tumor_ucell_obj, "/data/CDSL_hannenhalli/Cole/nmf_files/tumor_ucell_obj_combined.rds")

##### Plotting co-expression heatmap of factors.
factor_exp_mat <- dplyr::select(tumor_ucell_obj@meta.data,all_of(names(combined_factor_list))) %>% as.matrix
cor_mat <- cor(factor_exp_mat,method="pearson")     #Note that we use Pearson instead of Spearman for factors!
cor_df <- reshape2::melt(cor_mat) %>% dplyr::rename(cor=value)

# jaccard_df <- reshape2::melt(jaccard_mat) %>% dplyr::rename(jaccard_index=value)
# merged_df <- merge(cor_df,jaccard_df,by=c("Var1","Var2")) %>%
#   dplyr::filter(Var1 != Var2)

colors <- circlize::colorRamp2(c(min(cor_mat),0,max(cor_mat)),c("blue","white","red"))
num_split <- 7          #Here, I've chosen this value arbitrarily.
clust_obj <- hclust(as.dist(1-cor_mat),method="ward.D")     #Compute clustering of co-expression matrix.
cluster_assignment_vec <- cutree(clust_obj,k=num_split)

Heatmap(cor_mat, col=colors, row_split=num_split, column_split=num_split, 
         cluster_columns=clust_obj, cluster_rows=clust_obj,
         border_gp=gpar(lwd=0.5,col="black"),name="Pearson")


##### We explore how to merge co-expressed factors together. We use silhouette coefficients to explore 
#the number of clusters to pick for merging factors together.

cat("NUMBER OF COMBINED FACTORS: ",length(combined_factor_list),"\n")

num_split_vec <- 2:(length(combined_factor_list) - 1)     #Alter range of parameter sweep based on number of factors. Minimum is 2 and maximum is (# of factors - 1)
sil_df <- data.frame()
dist_obj <- as.dist(1-cor_mat)

for (num_split in num_split_vec) {
  cluster_assignment_vec <- cutree(clust_obj,k=num_split)
  sil_out <- summary(silhouette(cluster_assignment_vec,dist=dist_obj))
    
  cat(num_split,"\n")
    
  mean_sil_coef <- as.numeric(sil_out$si.summary["Mean"])
  sil_df <- rbind(sil_df,data.frame(num_split=num_split,sil=mean_sil_coef))
}

ggplot(sil_df) + geom_line(aes(x=num_split,y=sil,group=1)) +
  geom_point(aes(x=num_split,y=sil)) + theme_pubr(base_size=15)

# Set best split based on max sil coef
best_split <- sil_df$num_split[which(sil_df$sil == max(sil_df$sil))]

num_split <- best_split   #This number was chosen based on the parameter sweep above where the # of clusters with the highest silhouette coefficient was chosen.

cluster_assignment_vec <- cutree(clust_obj,k=num_split)
merged_factor_gene_list <- list()

# Merge together factors from the clustering. Uncomment internal code in the loop to add a further round of co-expression
#filtering. 
for (num in 1:num_split) {
  print(num)
  flush.console()
    
  factors_to_merge <- names(which(cluster_assignment_vec == num))
  merged_genes <- unique(unlist(combined_factor_list[factors_to_merge]))
  #     co_exp_mat <- cor(t(gene_exp_mat[merged_genes,]),method="spearman")
  #     co_exp_df <- reshape2::melt(co_exp_mat,value.name="cor_val") %>% na.omit %>% dplyr::filter(Var1 != Var2) %>% 
  #     dplyr::filter(cor_val >= cor_threshold)
  genes_to_retain <- merged_genes#unique(as.character(co_exp_df$Var1))
  merged_factor_name  <- paste("Merged_Factor",num,sep="_")
  
  if (length(genes_to_retain) >= min_factor_size) {
    merged_factor_gene_list[[merged_factor_name]] <- genes_to_retain
  }
}
#Print out the number of genes in each factor.
sapply(merged_factor_gene_list,length)

max_rank_to_use <- max(lengths(merged_factor_gene_list))+1

saveRDS(merged_factor_gene_list, "/data/CDSL_hannenhalli/Cole/nmf_files/merged_factor_list.rds")

#Score each tumor cell for expression of merged factors.

tumor_ucell_obj <- AddModuleScore_UCell(tumor_obj,features = merged_factor_gene_list,ncores = 6,maxRank = max_rank_to_use)
for (factor in names(merged_factor_gene_list)) {
  col_name <- grep(paste0(factor,"_"),colnames(tumor_ucell_obj@meta.data),value=T)
  tumor_ucell_obj[[factor]] <- tumor_ucell_obj[[col_name]]
  tumor_ucell_obj[[col_name]] <- NULL
}

saveRDS(tumor_ucell_obj, "/data/CDSL_hannenhalli/Cole/nmf_files/tumor_ucell_obj_merged.rds")

#Sample box-plot of expression of each merged factor.

df <- FetchData(tumor_ucell_obj,vars=c(names(merged_factor_gene_list),"sample_type","cell_type")) %>%
  tidyr::pivot_longer(cols=names(merged_factor_gene_list),names_to="factor",values_to="activity")
options(repr.plot.width=12,repr.plot.height=12)

ggplot( df ) + 
geom_boxplot(aes(x=sample_type,y=activity,fill=cell_type)) +
  facet_wrap(~factor,nrow=3) + theme_pubr(base_size=15) +
  theme(panel.grid.major=element_line(color="gray",linetype="dashed"))

options(repr.plot.width=7,repr.plot.height=7)


saveRDS(df, "/data/CDSL_hannenhalli/Cole/nmf_files/last_plot_df.rds")

