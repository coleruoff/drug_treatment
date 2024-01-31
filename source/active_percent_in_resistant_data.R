
# Watermelon data
data <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/watermelon_data/watermelon_pc9_processed.rds")

DimPlot(data[,data$time_point == 0])

scores <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/watermelon_pc9_processed_raj_watermelon_resistance_signature_aucell_scores.rds")
threshold <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/watermelon_pc9_processed_raj_watermelon_resistance_signature_aucell_thresholds.rds")

#Percent active of day 0 cells
sum(scores[rownames(scores) %in% colnames(data[,data$time_point == 0]),] > threshold$threshold)/ncol(data[,data$time_point == 0])

#Percent active of day 3 cells
sum(scores[rownames(scores) %in% colnames(data[,data$time_point == 3]),] > threshold$threshold)/ncol(data[,data$time_point == 3])

#Percent active of day 7 cells
sum(scores[rownames(scores) %in% colnames(data[,data$time_point == 7]),] > threshold$threshold)/ncol(data[,data$time_point == 7])

#Percent active of day 14 cells
sum(scores[rownames(scores) %in% colnames(data[,data$time_point == 14]),] > threshold$threshold)/ncol(data[,data$time_point == 14])

################################################################################
# Raj breast data
data <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/diverse_clonal_fates_data/raj_resistant_breast_processed.rds")

scores <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/raj_resistant_breast_processed_raj_watermelon_resistance_signature_aucell_scores.rds")
threshold <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/raj_resistant_breast_processed_raj_watermelon_resistance_signature_aucell_thresholds.rds")

sum(scores[,1] > threshold$threshold)/ncol(data)


################################################################################
# Raj melanoma data
data <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/diverse_clonal_fates_data/raj_resistant_melanoma_processed.rds")

scores <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/raj_resistant_melanoma_processed_raj_watermelon_resistance_signature_aucell_scores.rds")
threshold <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/raj_resistant_melanoma_processed_raj_watermelon_resistance_signature_aucell_thresholds.rds")

sum(scores[,1] > threshold$threshold)/ncol(data)


