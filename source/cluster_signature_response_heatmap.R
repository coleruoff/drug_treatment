library(ComplexHeatmap)

cell_lines <- c("A549","K562","MCF7")

RACs <- list(c(4,9,12,13,14,16,18),c(4,5,11),c(5,8,12,13,17))
names(RACs) <- c("A549","K562","MCF7")
clusters_of_interest <- RACs[[curr_cell_line]]

#################################################################################
signature_length <- 200

curr_cell_line <- "A549"
A549_cluster_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_cluster_signatures_",signature_length,".rds"))


curr_cell_line <- "K562"
K562_cluster_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_cluster_signatures_",signature_length,".rds"))


curr_cell_line <- "MCF7"
MCF7_cluster_signatures <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/", curr_cell_line, "_cluster_signatures_",signature_length,".rds"))

cell_line_cluster_signatures <- c(list(A549_cluster_signatures), list(K562_cluster_signatures), list(MCF7_cluster_signatures))
names(cell_line_cluster_signatures) <- cell_lines

#################################################################################

enlight_response_data <- read.csv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/drug_response_classifications_with_type.csv")

all_drugs <- unique(enlight_response_data$Dataset)

data_list <- list()
for(curr_drug in all_drugs){
  data_list[[curr_drug]] <- read.csv(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/enlight_data/", curr_drug,".csv"), row.names = 1)
}

curr_data <- data_list[[1]]

#################################################################################



for(curr_cell_line in cell_lines){
  
  curr_RACs <- RACs[[curr_cell_line]]
  
  
  for(curr_cluster in curr_RACs){
    
    curr_signature <- cell_line_cluster_signatures[[curr_cell_line]][[curr_cluster]]
    
    for(curr_data in data_list){
      
      filtered_data <- curr_data[rownames(curr_data) %in% curr_signature,]
      
      response_vector <- enlight_response_data %>% 
        filter(Sample.ID %in% gsub("\\.","-",colnames(filtered_data))) %>% 
        pull(Response)
      
      names(response_vector) <- enlight_response_data %>% 
        filter(Sample.ID %in% gsub("\\.","-",colnames(filtered_data))) %>% 
        pull(Sample.ID)
      
      
      colnames(filtered_data) <- gsub("\\.","-",colnames(filtered_data))
      rownames(filtered_data) <- NULL
      
      filtered_data <- filtered_data[,names(response_vector)]
      
      response_ha <- HeatmapAnnotation(response = response_vector,
                                       col = list(response = c("Responder" = "darkgreen", "Non-responder" = "lightblue")))
      
      
      
      Heatmap(filtered_data,bottom_annotation = response_ha)
      
      
    }
  }
}














