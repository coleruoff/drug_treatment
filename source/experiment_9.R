library(tidyverse)
library(ggpubr)
library(Seurat)
# install.packages("devtools")
# devtools::install_local("/data/CDSL_hannenhalli/Cole/CytoTRACE_0.3.3.tar.gz")
# library(CytoTRACE)

CytoTRACE <- function (mat, batch = NULL, enableFast = TRUE, ncores = 1, 
                       subsamplesize = 1000) {
  range01 <- function(x) {
    (x - min(x))/(max(x) - min(x))
  }
  a1 <- mat
  a2 <- batch
  if (ncol(mat) < 3000) {
    enableFast = FALSE
    message("The number of cells in your dataset is less than 3,000. Fast mode has been disabled.")
  }
  else {
    message("The number of cells in your dataset exceeds 3,000. CytoTRACE will now be run in fast mode (see documentation). You can multi-thread this run using the 'ncores' flag. To disable fast mode, please indicate 'enableFast = FALSE'.")
  }
  pqgenes <- is.na(Matrix::rowSums(mat > 0)) | apply(mat, 1, var) == 
    0
  num_pqgenes <- length(which(pqgenes == TRUE))
  mat <- mat[!pqgenes, ]
  if (num_pqgenes > 0) {
    warning(paste(num_pqgenes, "genes have zero expression in the matrix and were filtered"))
  }
  if (enableFast == FALSE) {
    size <- ncol(mat)
  }
  else if (enableFast == TRUE & subsamplesize < ncol(mat)) {
    size <- subsamplesize
  }
  else if (enableFast == TRUE & subsamplesize >= ncol(mat)) {
    stop("Please choose a subsample size less than the number of cells in dataset.")
  }
  chunk <- round(ncol(mat)/size)
  subsamples <- split(1:ncol(mat), sample(factor(1:ncol(mat)%%chunk)))
  message(paste("CytoTRACE will be run on", chunk, "sub-sample(s) of approximately", 
                round(mean(unlist(lapply(subsamples, length)))), "cells each using", 
                min(chunk, ncores), "/", ncores, "core(s)"))
  message(paste("Pre-processing data and generating similarity matrix..."))
  batches <- parallel::mclapply(subsamples, mc.cores = min(chunk, 
                                                           ncores), function(subsample) {
                                                             mat <- mat[, subsample]
                                                             batch <- batch[subsample]
                                                             if (max(mat) < 50) {
                                                               mat <- 2^mat - 1
                                                             }
                                                             if (length(grep("ERCC-", rownames(mat))) > 0) {
                                                               mat <- mat[-grep("ERCC-", rownames(mat)), ]
                                                             }
                                                             mat <- t(t(mat)/apply(mat, 2, sum)) * 1e+06
                                                             pqcells <- is.na(apply(mat > 0, 2, sum)) | apply(mat > 
                                                                                                                0, 2, sum) <= 10
                                                             num_pqcells <- length(which(pqcells == TRUE))
                                                             mat <- mat[, !pqcells]
                                                             mat <- log(mat + 1, 2)
                                                             mat <- data.matrix(mat)
                                                             counts <- apply(mat > 0, 2, sum)
                                                             if (ncol(a1) == length(a2)) {
                                                               batch <- batch[!pqcells]
                                                               suppressMessages(mat <- sva::ComBat(mat, batch, 
                                                                                                   c()))
                                                               mat <- data.matrix(mat)
                                                               mat[which(mat < 0)] <- 0
                                                             }
                                                             census_normalize <- function(mat, counts) {
                                                               xnl <- 2^data.matrix(mat) - 1
                                                               rs <- apply(xnl, 2, sum)
                                                               rnorm <- t(t(xnl) * counts/rs)
                                                               A <- log(rnorm + 1, 2)
                                                               return(A)
                                                             }
                                                             mat2 <- census_normalize(mat, counts)
                                                             mvg <- function(matn) {
                                                               A <- matn
                                                               n_expr <- Matrix::rowSums(A > 0)
                                                               A_filt <- A[n_expr >= 0.05 * ncol(A), ]
                                                               vars <- apply(A_filt, 1, var)
                                                               means <- apply(A_filt, 1, mean)
                                                               disp <- vars/means
                                                               last_disp <- tail(sort(disp), 1000)[1]
                                                               A_filt <- A_filt[disp >= last_disp, ]
                                                               return(A_filt)
                                                             }
                                                             mat2.mvg <- mvg(mat2)
                                                             rm1 <- colSums(mat2.mvg) == 0
                                                             mat2 <- mat2[, !rm1]
                                                             counts <- counts[!rm1]
                                                             similarity_matrix_cleaned <- function(similarity_matrix) {
                                                               D <- similarity_matrix
                                                               cutoff <- mean(as.vector(D))
                                                               diag(D) <- 0
                                                               D[which(D < 0)] <- 0
                                                               D[which(D <= cutoff)] <- 0
                                                               Ds <- D
                                                               D <- D/Matrix::rowSums(D)
                                                               D[which(Matrix::rowSums(Ds) == 0), ] <- 0
                                                               return(D)
                                                             }
                                                             D <- similarity_matrix_cleaned(HiClimR::fastCor(mvg(mat2)))
                                                             return(list(mat2 = mat2, counts = counts, D = D))
                                                           })
  cat(class(batches[1]),"\n")
  cat(str(batches[1]))
  mat2 <- do.call(cbind, lapply(batches, function(x) x$mat2))
  counts <- do.call(c, lapply(batches, function(x) x$counts))
  filter <- colnames(a1)[-which(colnames(a1) %in% colnames(mat2))]
  if (length(filter) > 0) {
    warning(paste(length(filter), "poor quality cells were filtered based on low or no expression. See 'filteredCells' in returned object for names of filtered cells."))
  }
  message("Calculating gene counts signature...")
  ds2 <- sapply(1:nrow(mat2), function(x) ccaPP::corPearson(mat2[x, 
  ], counts))
  names(ds2) <- rownames(mat2)
  gcs <- apply(mat2[which(rownames(mat2) %in% names(rev(sort(ds2))[1:200])), 
  ], 2, mean)
  samplesize <- unlist(lapply(lapply(batches, function(x) x$counts), 
                              length))
  gcs2 <- split(gcs, as.numeric(rep(names(samplesize), samplesize)))
  D2 <- lapply(batches, function(x) x$D)
  regressed <- function(similarity_matrix_cleaned, score) {
    out <- nnls::nnls(similarity_matrix_cleaned, score)
    score_regressed <- similarity_matrix_cleaned %*% out$x
    return(score_regressed)
  }
  diffused <- function(similarity_matrix_cleaned, score, ALPHA = 0.9) {
    vals <- score
    v_prev <- rep(vals)
    v_curr <- rep(vals)
    for (i in 1:10000) {
      v_prev <- rep(v_curr)
      v_curr <- ALPHA * (similarity_matrix_cleaned %*% 
                           v_curr) + (1 - ALPHA) * vals
      diff <- mean(abs(v_curr - v_prev))
      if (diff <= 1e-06) {
        break
      }
    }
    return(v_curr)
  }
  message("Smoothing values with NNLS regression and diffusion...")
  cytotrace <- parallel::mclapply(1:length(D2), mc.cores = ncores, 
                                  function(i) {
                                    gcs_regressed <- regressed(D2[[i]], gcs2[[i]])
                                    gcs_diffused <- diffused(D2[[i]], gcs_regressed)
                                    cytotrace <- rank(gcs_diffused)
                                  })
  cytotrace <- cytotrace_ranked <- unlist(cytotrace)
  cytotrace <- range01(cytotrace)
  cytogenes <- sapply(1:nrow(mat2), function(x) ccaPP::corPearson(mat2[x, 
  ], cytotrace))
  names(cytogenes) <- rownames(mat2)
  message("Calculating genes associated with CytoTRACE...")
  names(cytotrace) <- names(cytotrace_ranked) <- names(gcs) <- names(counts) <- colnames(mat2)
  cytotrace <- cytotrace[colnames(a1)]
  cytotrace_ranked <- cytotrace_ranked[colnames(a1)]
  gcs <- gcs[colnames(a1)]
  counts <- counts[colnames(a1)]
  mat2 <- t(data.frame(t(mat2))[colnames(a1), ])
  names(cytotrace) <- names(cytotrace_ranked) <- names(gcs) <- names(counts) <- colnames(mat2) <- colnames(a1)
  message("Done")
  return(list(CytoTRACE = cytotrace, CytoTRACErank = cytotrace_ranked, 
              cytoGenes = sort(cytogenes, decreasing = T), GCS = gcs, 
              gcsGenes = sort(ds2, decreasing = T), Counts = counts, 
              filteredCells = filter, exprMatrix = mat2))
}

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"
cell_lines <- c("A549","K562","MCF7")

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- cell_lines

emergent <- list(c(14:19),c(9,11),c(13,15,18))
names(emergent) <- cell_lines

curr_cell_line <- "MCF7"
data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/",curr_cell_line, "_processed_filtered.rds"))

#read in DR signature scores and set active cells
scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
threshold <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
clusters_of_interest <- RACs[[curr_cell_line]]

#Add metadata for RAC and Cell Group
data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, "1", ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), 2, 0)), col.name = "cell_group")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & colnames(data) %in% active_cell_names, paste0(data$Cluster, "_1"), ifelse(data$rac == "rac" & (!colnames(data) %in% active_cell_names), paste0(data$Cluster, "_2"), paste0(data$Cluster, "_0"))), col.name = "cell_cluster_group")
data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & data$Cluster %in% emergent[[curr_cell_line]], "emergent_rac", ifelse(data$rac == "rac" & (!data$Cluster %in% emergent[[curr_cell_line]]), "non_emergent_rac", "non_rac")), col.name = "emergent_rac")


results <- CytoTRACE(as.matrix(data@assays$RNA$data), ncores=10)
saveRDS(results, paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cytotrace_results/",curr_cell_line, "_cytotrace_results.rds"))

# results <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cytotrace_results/",curr_cell_line, "_cytotrace_results.rds"))
# 
# data <- AddMetaData(data,results$CytoTRACE,col.name="cytotrace")
# 
# #############
# 
# df <- data@meta.data
# 
# p <- ggboxplot(df, x = "emergent_rac", y = "cytotrace",fill = "rac")
# 
# my_comparisons <- list( c("rac", "nonrac"))
# 
# plot_title <- "CytoTRACE Score Distributions"
# 
# p <- p + stat_compare_means(comparisons = my_comparisons,label = "p.format", method = "wilcox", label.x = 2.2, size=8)+
#   ggtitle(plot_title)+
#   xlab("")+
#   ylab("G1S Score")+
#   scale_fill_manual(values=c("lightblue", "pink"),name = "Cell Groups")+
#   theme(legend.position="right",
#         title = element_text(size=20, face = "bold"),
#         axis.text = element_text(size=20),
#         legend.text = element_text(size=24),
#         legend.title = element_text(size=26),
#         legend.key.height = unit(1.5,"cm"),
#         legend.key.width = unit(1.5,"cm"))
# 
# p
