source("../survival_analysis/cox_regression.R")
library(GSVA)

# clinical_sample <- read.table("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/metabric_data/brca_metabric/data_clinical_sample.txt",
#                               fill=T)
# colnames(clinical_sample) <- clinical_sample[1,]
# clinical_sample <- clinical_sample[-1,]


clinical_data <- read_tsv("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/metabric_data/brca_metabric_clinical_data.tsv")

patient_sample <- read.table("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/metabric_data/brca_metabric/data_clinical_patient.txt",
           fill=T)

colnames(patient_sample) <- patient_sample[1,]
patient_sample <- patient_sample[-1,]

temp <- matrix(unlist(patient_sample),nrow=nrow(patient_sample))
colnames(temp) <- colnames(patient_sample)
rownames(temp) <- rownames(patient_sample)

patient_sample <- temp

patient_sample <- as.data.frame(patient_sample) %>% 
  filter(CLAUDIN_SUBTYPE == "LumA")


microarray_data <- read.table("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/metabric_data/brca_metabric/data_mrna_illumina_microarray.txt",
                              fill=T)

colnames(microarray_data) <- microarray_data[1,]
microarray_data <- microarray_data[-1,]

rownames(microarray_data) <- make.names(microarray_data[,1],unique=T)
microarray_data <- microarray_data[,-c(1,2)]

class(temp[,1])

class(microarray_data)

temp <- matrix(as.numeric(unlist(microarray_data)),nrow=nrow(microarray_data))
colnames(temp) <- colnames(microarray_data)
rownames(temp) <- rownames(microarray_data)

microarray_data <- temp

##########################################################################
type_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_rac_type_signatures.rds")

type_signatures <- type_signatures[grepl("type1",names(type_signatures))]

curr_signature <- type_signatures[3]
##########################################################################

ssgsea_result <- gsva(microarray_data, curr_signature, method='ssgsea')

gene_set_scores_df <- as.data.frame(t(ssgsea_result)) %>%
  tibble::rownames_to_column("PATIENT_ID")


patient_sample$OS_MONTHS <- as.numeric(patient_sample$OS_MONTHS)
patient_sample$OS_STATUS <- as.numeric(gsub(":[A-Z]*", "", patient_sample$OS_STATUS))
patient_sample$AGE_AT_DIAGNOSIS <- as.numeric(patient_sample$AGE_AT_DIAGNOSIS)


cox_regression_info <- cox_regression(sample_survival_df = patient_sample,
                                      survival_data_id_col = "PATIENT_ID",
                                      duration_str = "OS_MONTHS",
                                      feature_data_id_col = "PATIENT_ID",
                                      model_covariates = c("AGE_AT_DIAGNOSIS"),
                                      status_str="OS_STATUS",
                                      sample_features_df = gene_set_scores_df, 
                                      km_features_to_plot=names(curr_signature))


modform <- as.formula(paste0("Surv(OS_MONTHS, OS_STATUS) ~ ", names(curr_signature)))

fit <- eval(substitute(survfit(modform, data = cox_regression_info$km_df), list(modform = modform)))


cox_regression_info$km_df[[names(curr_signature)]] <- factor(cox_regression_info$km_df[[names(curr_signature)]], levels = c("High","Medium","Low"))


# cox_regression_info$km_df$MCF7_global_rac_type1_signature <- "Low"

p <- ggsurvplot(fit, data = cox_regression_info$km_df,
                pval = T,
                pval.size = 6,
                legend.title="Expression Level:",
                legend.labs= c("High","Medium","Low"),
                font.title=c(28),
                font.x=c(28),
                font.y=c(28),
                font.tickslab=c(20),
                font.legned=c(50),
                font.caption=c(30),
                tables.theme = clean_theme())+
  xlab("")+
  ylab("")


p$plot <- p$plot + 
  theme(legend.text = element_text(size = 24, color = "black", face = "bold"),
        legend.title = element_text(size = 28, color = "black", face = "bold"),
        legend.key.height = unit(1.5,"cm"),
        legend.key.width = unit(1.5,"cm"))


p

cox_regression_info$regression_df$hazard_ratio















