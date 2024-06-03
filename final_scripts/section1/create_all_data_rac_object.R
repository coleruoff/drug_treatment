args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
setwd(args[1])

library(Seurat)
library(tidyverse)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"


CellCycleScoring <- function (object, s.features, g2m.features, g0.features, ctrl = NULL, set.ident = FALSE, 
                              ...) 
{
  name <- "Cell.Cycle"
  features <- list(S.Score = s.features, G2M.Score = g2m.features, G0.Score = g0.features)
  if (is.null(x = ctrl)) {
    ctrl <- min(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1)))
  }
  object.cc <- AddModuleScore(object = object, features = features, 
                              name = name, ctrl = ctrl, ...)
  cc.columns <- grep(pattern = name, x = colnames(x = object.cc[[]]), 
                     value = TRUE)
  cc.scores <- object.cc[[cc.columns]]
  rm(object.cc)
  CheckGC()
  assignments <- apply(X = cc.scores, MARGIN = 1, FUN = function(scores, 
                                                                 first = "S", second = "G2M", third = "G0", null = "G1") {
    if (all(scores < 0)) {
      return(null)
    }
    else {
      if (length(which(x = scores == max(scores))) > 1) {
        return("Undecided")
      }
      else {
        return(c(first, second, third)[which(x = scores == 
                                               max(scores))])
      }
    }
  })
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments), 
                     by = 0)
  colnames(x = cc.scores) <- c("rownames", "S.Score", "G2M.Score", "G0.Score", 
                               "Phase")
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c("S.Score", "G2M.Score", "G0.Score", "Phase")]
  object[[colnames(x = cc.scores)]] <- cc.scores
  if (set.ident) {
    object[["old.ident"]] <- Idents(object = object)
    Idents(object = object) <- "Phase"
  }
  return(object)
}

################################################################################

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
g0.genes <- readRDS(paste0(dataDirectory, "genesets/g0_signature.rds"))
g0.genes <- g0.genes$g0_signature

cell_lines <- c("A549","K562","MCF7")

RACs <- readRDS(paste0(dataDirectory, "processed_data/all_RACs.rds"))

# emergent <- list(c(14:19),c(9),c(13,15,18))
# names(emergent) <- cell_lines

all_data <- list()
for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line,"\n")
  
  data <- readRDS(paste0(dataDirectory, "processed_data/sciPlex_data/",curr_cell_line, "_processed_filtered.rds"))
  
  #read in DR signature scores and set active cells
  scores <- readRDS(paste0(dataDirectory, "aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_scores.rds"))
  threshold <- readRDS(paste0(dataDirectory, "aucell_score_objects/", curr_cell_line, "_processed_filtered_raj_watermelon_resistance_signature_aucell_thresholds.rds"))
  active_cell_names <- rownames(scores)[scores[,1] > threshold$threshold]
  clusters_of_interest <- RACs[[curr_cell_line]]
  
  #Add metadata for RAC and emergent clusters
  data <- AddMetaData(data, metadata = scores, col.name = "resistance_score")
  data <- AddMetaData(data, metadata = ifelse(colnames(data) %in% active_cell_names, "active","inactive"),col.name = "resistance_active")
  data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% clusters_of_interest, "rac","nonrac"), col.name = "rac")
  # data <- AddMetaData(data, metadata = ifelse(data$rac == "rac" & data$Cluster %in% emergent[[curr_cell_line]], "emergent_rac", ifelse(data$rac == "rac" & (!data$Cluster %in% emergent[[curr_cell_line]]), "non_emergent_rac", "non_rac")), col.name = "emergent_rac")
  # data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% emergent[[curr_cell_line]], "emergent", "non_emergent"), col.name = "emergent")
  
  data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, g0.features = g0.genes, set.ident = F)
  
  all_data[[curr_cell_line]] <- data
}


saveRDS(all_data, paste0(dataDirectory, "processed_data/sciPlex_data/all_cell_lines_data.rds"))

