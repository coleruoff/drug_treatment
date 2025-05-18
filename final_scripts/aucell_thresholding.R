library(data.table)
library(Seurat)
library(dplyr)
library(AUCell)

#Construct a background gene set collection for a given foreground gene set collection.
#Inputs :
#gene_exp_mat : A matrix of genes x cells. 
#num_bins : Number of bins to distribution genes based on their mean expression across cells. Default : 10. Tuned to
#typical 10X datsets, and can be increased if the number of expressed genes/cell is large, such as in SMART-seq data.
#gene_universe : The set of genes from which the background must be chosen. Default : NULL. In some applications, it may
#be desirable to set this to be identical to the set of genes amongst the foreground gene set collection.
#
#Returns : 
#A list of gene sets where the following two conditions are satisified : 
#1. Each background gene set is the same size as its corresponding foreground gene set.
#2. The mean expression of genes within each background gene set are matched to the corresponding foreground gene set.
find_control_gene_sets <- function(gene_exp_mat,gene_sets,num_bins=10,gene_universe=NULL) {
    #Compute mean expression of each gene across cells
    mean_gene_exp_vec <- rowMeans(gene_exp_mat)
    genes_by_bin <- list()
    if (is.null(gene_universe)) {
        gene_universe <- rownames(gene_exp_mat)#rownames(gene_exp_mat)
    }

    #Assign bins to each gene based on its mean expression==
    
    gene_exp_df <- tibble::enframe(mean_gene_exp_vec,name="gene",value="gene_exp") %>% 
    mutate(bin=ntile(gene_exp,num_bins)) %>%
    dplyr::filter(gene %in% gene_universe)

    for (exp_bin in 1:num_bins) {
        genes_by_bin[[exp_bin]] <- gene_exp_df %>% dplyr::filter(bin == exp_bin) %>% pull(gene)
    }

    control_gene_sets <- list()
    for (gene_set_name in names(gene_sets)) {
        gene_set <- gene_sets[[gene_set_name]]
        
        bin_dist_df <- gene_exp_df %>% dplyr::filter(gene %in% gene_set) %>% group_by(bin) %>% dplyr::count()
        control_genes <- c()
        for (idx in 1:nrow(bin_dist_df)) {
            bin_num <- bin_dist_df[idx,] %>% pull(bin)
            num_genes <- bin_dist_df[idx,] %>% pull(n)
            control_genes <- c(control_genes,sample(genes_by_bin[[bin_num]],size=num_genes))
        }
        control_gene_sets[[paste(gene_set_name,"Control",sep="_")]] <- control_genes
    }
    
    return(control_gene_sets) 
}


#Compute AUcell scores of a given gene set collection.
#Arguments : 
#seurat_obj_ - A Seurat object. NormalizeData/ScaleData/FindVariableFeatures/PCA does not have to be run. 
#gene_sets - A named list of gene sets, with each list member containing a vector of gene names.
#compute_thresholds - Whether to ask AUCell to compute thresholds. Default : False.
#threshold_type - If compute_thresholds is set to True, this specifies the type of threshold returned by AUCell. Default
#: Global_k1. This is computed by fitting a bimodal distribution to AUCell scores. See AUCell documentation for other
#options.
#rankings_obj - This is the ranking object created by AUCell during calculation. This is useful when repeatedly running
#compute_AUCell_scores in a loop, as ranking object construction can be an expensive step. Default : NULL
#assay_to_use - Seurat object assay to use. Default : RNA. Using the 'integrated' assay is not recommended as those
#read counts are not meant for any downstream analyses.
#nCores - Number of cores to use. Default : 3. Increasing this number drastically increases the memory usage of the code.
#Returns : 
#A list with two members.
#-   auc_mat : Contains a matrix of AUCell scores 
#-  thresholds : If compute_threholds is set to True, then a vector of thresholds, one per gene set, is returned. Else,
#   NULL is returned.
compute_AUCell_scores <-
function(seurat_obj=NULL,gene_sets=NULL,compute_thresholds=F,threshold_type="Global_k1",rankings_obj=NULL,assay_to_use="RNA",
nCores=3) {
    
    gene_set_names <- names(gene_sets)
    to_return <- list()
    if (!is.null(rankings_obj)) {
        cells_rankings <- rankings_obj
    } else {
        
        cells_rankings <- AUCell_buildRankings(seurat_obj[[assay_to_use]]$counts, nCores=nCores, plotStats=F,verbose=F,splitByBlocks = T)
        
        to_return$rankings <- cells_rankings
    }
    cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings,verbose=F,nCores=nCores)
    aucell_scores_mat <- t(getAUC(cells_AUC))
    to_return[["auc_mat"]] =aucell_scores_mat

    to_return[["thresholds"]] <- NULL
    if (compute_thresholds == T) {
        cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F )
        auc_thr <- sapply(cells_assignment, function(x){return(x$aucThr$thresholds[threshold_type,"threshold"])})
        to_return[["thresholds"]] <- auc_thr
    }

    return(to_return)
}

#Compute AUcell scores of randomly generated gene sets in a single-cell data set for a given gene set collection.
#Arguments : 
#seurat_obj_ - A Seurat object. NormalizeData needs to be run as normalized counts are required for background gene set
#construction. 
#gene_sets - A named list of gene sets, with each list member containing a vector of gene names.
#sample_info_column - The column in the Seurat meta.data that represents the 'batch' variable/co-variate in the analysis. This could
#be donor ID, sample number, etc. Default : orig.ident. 
#do_sample_wise - If set to TRUE, compute AUCell thresholds separately for each batch of cells. Default : TRUE
#num_controls - Number of background gene sets to generate.
#num_bins - Number of gene expression bins within which to classify genes based on their mean expression across all
#cells. Default : 10. This is tuned to 10X data, and can be increased if the # of expressed genes/per cell tends to be
#over 5000 genes/cell.
#q_thresh - Percentile of background gene set activity distribution to be used as a threshold for each gene set. Default : 1 i.e., the
#maximum value observed.
#nCores - Number of cores to use. Default : 3. Increasing this number drastically increases the memory usage of the code.
#
#Returns : 
#A data frame of AUCell thresholds each data set or (data set, batch) combination depending on whether do_sample_wise
#was set to F or T respectively.
compute_shuffled_gene_set_AUCell_scores <-
function(seurat_obj_,gene_sets,sample_info_column="orig.ident",do_sample_wise=T, num_controls=100, num_bins=10,
q_thresh=1.0,nCores=3,gene_universe=NULL, assay_to_use="RNA") {
    control_sd_df <- data.frame()
    rankings_obj_list <- list()
    control_gene_set_list <- list()

    samples <- unique(seurat_obj_@meta.data[[sample_info_column]])
    for (control in 1:num_controls) {
        print(control)
        flush.console()

        if (do_sample_wise) {
            sd_df <- data.frame()
            for (sample_name in samples) {
                cells <- rownames(seurat_obj_@meta.data %>% dplyr::filter(!!sym(sample_info_column) == sample_name))
                if (!sample_name %in% names(control_gene_set_list)) {
                    control_gene_set_list[[sample_name]] <- find_control_gene_sets(seurat_obj_[[assay_to_use]]$data[,cells],gene_sets)
                }

                if (control == 1) {
                    auc_output <- compute_AUCell_scores(seurat_obj_[,cells],
                                                        compute_thresholds = F,
                                                        control_gene_set_list[[sample_name]],rankings_obj=NULL,
                                                        nCores=nCores, assay_to_use = assay_to_use)
                    rankings_obj_list[[sample_name]] <- auc_output$rankings
                } else {
                    auc_output <- compute_AUCell_scores(rankings_obj=rankings_obj_list[[sample_name]],
                                                        compute_thresholds = F,gene_sets=control_gene_set_list[[sample_name]],
                                                        nCores=nCores, assay_to_use = assay_to_use)
                }

                temp_df <- apply(as.matrix(auc_output$auc_mat[cells,]),2,sd) %>%
                tibble::enframe(.,name="gene_set",value="stdev") %>%
                mutate(shuffle=control,sample=sample_name,mean=apply(as.matrix(auc_output$auc_mat[cells,]),2,
                function(x){return(quantile(x,q_thresh))}))
                sd_df <- rbind(sd_df,temp_df)
            }
        } else {
                control_gene_set_list <-
                find_control_gene_sets(seurat_obj_[[assay_to_use]]$data,gene_sets,gene_universe=gene_universe)
                if (control == 1) {
                    auc_output <- compute_AUCell_scores(seurat_obj_,
                                                        compute_thresholds = F,
                                                        control_gene_set_list,rankings_obj=NULL,
                                                        nCores=nCores, assay_to_use = assay_to_use)
                    rankings_obj <- auc_output$rankings
                } else {
                    auc_output <- compute_AUCell_scores(rankings_obj=rankings_obj,
                                                        compute_thresholds = F,gene_sets=control_gene_set_list,
                                                        nCores=nCores, assay_to_use = assay_to_use)
                }
                sd_df <- apply(as.matrix(auc_output$auc_mat),2,sd) %>%
                tibble::enframe(.,name="gene_set",value="stdev") %>% 
                mutate(mean=apply(as.matrix(auc_output$auc_mat),2,function(x){return(quantile(x,q_thresh))}),shuffle=control)
        }
        control_sd_df <- rbind(control_sd_df,sd_df)
    }

    if (do_sample_wise) {
       control_sd_df <- group_by(control_sd_df,gene_set,sample) %>% summarize(threshold=max(mean)) %>%
       mutate(gene_set=gsub("_Control","",gene_set)) %>% dplyr::rename(!!sym(sample_info_column):=sample)
    } else {
       control_sd_df <- group_by(control_sd_df,gene_set) %>% summarize(threshold=max(mean)) %>%
       mutate(gene_set=gsub("_Control","",gene_set))
    }

    return(control_sd_df)
}

