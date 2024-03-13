library(Seurat)

raj_resistant_breast <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/diverse_clonal_fates_data/raj_resistant_breast_processed.rds")
raj_resistant_melanoma <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/drug_treatment/processed_data/diverse_clonal_fates_data/raj_resistant_melanoma_processed.rds"))

raj_control_breast <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/drug_treatment/processed_data/memorySeq_data/MDA_control_tpm.rds"))
raj_control_melanoma <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/drug_treatment/processed_data/memorySeq_data/WM989_control_tpm.rds"))

files <- list.files("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/watermelon_data/")

files <- files[1:8]
total_df <- list()
common_cols <- c()
for(curr_file in files){
  cat(curr_file, "\n")
  watermelon_data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/watermelon_data/",curr_file))
  
  watermelon_data <- watermelon_data[VariableFeatures(watermelon_data),]
  watermelon_data$time_point <- ifelse(grepl("_10_", watermelon_data$sample_pool), 1, 0)
  
  curr_df <- cbind(t(as.matrix(watermelon_data@assays$RNA$data)),watermelon_data$time_point)
  colnames(curr_df)[ncol(curr_df)] <- "labels"
  
  total_df <- append(total_df, list(curr_df))
  
  if(which(curr_file==files) == 1){
    common_cols <- colnames(curr_df)
  } else {
    common_cols <- intersect(common_cols,colnames(curr_df))
  }
  
}

#Add raj data
temp <- cbind(t(raj_control_breast),0)
colnames(temp)[ncol(temp)] <- "labels"
total_df <- append(total_df, list(temp))
common_cols <- intersect(common_cols,colnames(temp))

temp <- cbind(t(raj_control_melanoma),0)
colnames(temp)[ncol(temp)] <- "labels"
total_df <- append(total_df, list(temp))
common_cols <- intersect(common_cols,colnames(temp))

temp <- cbind(t(as.matrix(raj_resistant_breast@assays$RNA@data)),1)
colnames(temp)[ncol(temp)] <- "labels"
total_df <- append(total_df, list(temp))
common_cols <- intersect(common_cols,colnames(temp))

temp <- cbind(t(as.matrix(raj_resistant_melanoma@assays$RNA@data)),1)
colnames(temp)[ncol(temp)] <- "labels"
total_df <- append(total_df, list(temp))
common_cols <- intersect(common_cols,colnames(temp))

# add day 14 watermelon data
watermelon_data <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/watermelon_data/watermelon_pc9_processed.rds")
watermelon_data <- watermelon_data[,watermelon_data$time_point %in% c(0,14)]
temp <- cbind(t(as.matrix(watermelon_data@assays$RNA@data)),ifelse(watermelon_data$time_point == 0, 0, 1))
colnames(temp)[ncol(temp)] <- "labels"
total_df <- append(total_df, list(temp))

common_cols <- intersect(common_cols,colnames(temp))

total_df <- lapply(total_df, FUN = function(x) subset(x, select=common_cols))

total_df <- do.call("rbind",total_df)

cancer_data <- total_df
colnames(cancer_data)[ncol(cancer_data)] <- "labels"

set.seed(42)  # for reproducibility
train_idx <- sample(nrow(cancer_data), 0.7 * nrow(cancer_data))
train_data <- cancer_data[train_idx, ]
test_data <- cancer_data[-train_idx, ]

train_data <- data.frame(train_data)

train_data <- lapply(train_data, unlist)

# Create the logistic regression model
model <- glm(labels ~ ., data = train_data, family = binomial(link='logit'))

saveRDS(model, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/experiment_data/experiment_22_model.rds")

# model <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/experiment_data/experiment_22_model.rds")
# 
# 
# summary(model)
# 
# cancer_data["labels",]
# # Make predictions on the testing set
# predictions <- predict(model, newdata = data.frame(test_data), type = "response")
# 
# # Convert probabilities to binary predictions
# threshold <- 0.5
# binary_predictions <- ifelse(predictions > threshold, 1, 0)
# 
# # Evaluate the model
# accuracy <- mean(binary_predictions == test_data[,"labels"])
# 
# 
# 
# temp <- summary(model)
# names(sort(temp$coefficients[,"Estimate"],decreasing = T)[1:200])


