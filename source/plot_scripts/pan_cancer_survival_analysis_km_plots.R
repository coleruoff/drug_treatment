library(tidyverse)
library(gt)
library(GSVA)
library(readxl)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(survminer)
library(survival)
library(patchwork)
source("/data/ruoffcj/projects/survival_analysis/cox_regression.R")


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


i <- "TCGA-LUAD"
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
  GDCdownload(query_tcga,
              directory = "/data/CDSL_hannenhalli/Cole/GDCdata/")
  
  # Prepare data
  tcga_data <- GDCprepare(query_tcga, summarizedExperiment = T,
                          directory = "/data/CDSL_hannenhalli/Cole/GDCdata/")
  
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

################################################################################
# Get list of projects
##################################################################################

gdcprojects <- getGDCprojects()

tcga_projects <- gdcprojects[grepl("TCGA", gdcprojects$id),]$id

tcga_projects <- tcga_projects

##################################################################################


type1_supercluster_down_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/type1_supercluster_down_signatures.rds")

type1_supercluster_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/type1_supercluster_signatures.rds")

all_signatures <- type1_supercluster_signatures

##################################################################################

metric_to_use <- "OS"

hazard_ratio_df <- list()

curr_project <- tcga_projects[1]

curr_project <- "TCGA-UCEC"
for(curr_project in tcga_projects){
  
  cat(curr_project, "\n")
  
  #Get read count data
  # tcga_matrix_vst <- get_read_count_data(curr_project)
  tcga_matrix_vst <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/TCGA_processed_data/",curr_project,"_data.rds"))
  
  # Score samples
  ssgsea_result <- gsva(tcga_matrix_vst,all_signatures,method='ssgsea')
  
  gene_set_scores_df <- as.data.frame(t(ssgsea_result)) %>%
    tibble::rownames_to_column("submitter_id")
  
  gene_set_scores_df$submitter_id <- sapply(gene_set_scores_df$submitter_id, function(x) gsub("\\.", "-",x))
  gene_set_scores_df$submitter_id <- sapply(gene_set_scores_df$submitter_id, function(x) substring(x, 1,12))
  
  # Get clinical data
  clinical_df <- get_clinical_data(curr_project)
  
  #Run cox regression
  if("purity" %in% colnames(clinical_df)){
    covariates <- c("purity", "sex", "age")
  } else {
    covariates <- c("sex", "age")
  }
  
  
  clinical_df$sex <- factor(clinical_df$sex, levels = c("FEMALE","MALE"))
 
  
  cox_regression_info <- cox_regression(sample_survival_df = clinical_df,
                                        survival_data_id_col = "submitter_id",
                                        duration_str = paste0(metric_to_use, ".time"),
                                        feature_data_id_col = "submitter_id",
                                        model_covariates = covariates,
                                        status_str=metric_to_use,
                                        sample_features_df = gene_set_scores_df, 
                                        km_features_to_plot=names(all_signatures))
  
  # supercluster1_hazard_ratios <- append(supercluster1_hazard_ratios, cox_regression_info$regression_df$hazard_ratio[1])
  # supercluster2_hazard_ratios <- append(supercluster2_hazard_ratios, cox_regression_info$regression_df$hazard_ratio[2])
  
  
  plots <- list()
  for(curr_signature in 1:nrow(cox_regression_info$regression_df)){
    
    signature_to_use <- names(all_signatures)[curr_signature]
    signatures_title <- gsub("_"," ",signature_to_use)
    signatures_title <- gsub("type1 ","",signatures_title)
    signatures_title <- gsub(" signature","",signatures_title)
    signatures_title <- gsub("supercluster","Supercluster ",signatures_title)
    
    
    # Cox regression
    regression_str <- paste0("Surv(",metric_to_use,".time, ",metric_to_use,") ~ ", signature_to_use, " + ", paste(covariates,collapse = " + "))
    data_df <- as.data.frame(merge(clinical_df,gene_set_scores_df, by='submitter_id'))

    data_df <- data_df %>%
      select(covariates,names(all_signatures), metric_to_use, paste0(metric_to_use,".time"))

    data_df$sex <- factor(data_df$sex, levels = c("FEMALE","MALE"))
    
    model <- coxph(as.formula(regression_str), data = data_df)
    model_sum <- summary(model)
    
    hazard_ratio_df[["project"]] <- append(hazard_ratio_df[["project"]],curr_project)
    hazard_ratio_df[["hazard_ratio"]] <- append(hazard_ratio_df[["hazard_ratio"]], model_sum$conf.int[signature_to_use,"exp(coef)"])
    hazard_ratio_df[["hazard_ratio.high"]] <- append(hazard_ratio_df[["hazard_ratio.high"]],model_sum$conf.int[signature_to_use,"upper .95"])
    hazard_ratio_df[["hazard_ratio.low"]] <- append(hazard_ratio_df[["hazard_ratio.low"]],model_sum$conf.int[signature_to_use,"lower .95"])
    hazard_ratio_df[["p_value"]] <- append(hazard_ratio_df[["p_value"]], model_sum$coefficients[signature_to_use, "Pr(>|z|)"])
    hazard_ratio_df[["signature"]] <- append(hazard_ratio_df[["signature"]], paste0("Supercluster ",curr_signature))
    
    
    
    #KM Plots
    modform <- as.formula(paste0("Surv(",metric_to_use,".time, ",metric_to_use,") ~ ", signature_to_use))
    
    fit <- eval(substitute(survfit(modform, data = cox_regression_info$km_df), list(modform = modform)))
    
    
    cox_regression_info$km_df[[signature_to_use]] <- factor(cox_regression_info$km_df[[signature_to_use]], levels = c("High","Medium","Low"))
    
    p <- ggsurvplot(fit, data = cox_regression_info$km_df,
                    pval = T, 
                    legend.title="Expression Level:",
                    legend.labs= c("High","Medium","Low"),
                    pval.size = 6,
                    font.title=c(28),
                    font.x=c(28),
                    font.y=c(28),
                    font.tickslab=c(20),
                    font.legned=c(50),
                    font.caption=c(30),
                    tables.theme = clean_theme())
    
    p$plot <- p$plot + 
      theme(legend.text = element_text(size = 14, color = "black", face = "bold"),
            legend.title = element_text(size = 14, color = "black", face = "bold"))
    
    p <- p + ggtitle(signatures_title)
    
    plots <- append(plots, list(p$plot))
    
  }
  
  if(metric_to_use == "OS"){
    main_title <- paste0("RAC Type 1 Supercluster Signatures - Overall Survival (",curr_project,")")
    file_name <- paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/survival_plots/pan_cancer/",curr_project,"_OS.png")
  } else {
    main_title <- paste0("RAC Type 1 Supercluster Signatures - Progression Free Survival (",curr_project,")")
    file_name <- paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/survival_plots/pan_cancer/",curr_project,"_PFS.png")
  }
  
  # png(file_name, height = 1000, width = 2000)
  
  final_plot <- ggarrange(plotlist = plots,ncol = 2,nrow=1, common.legend = T)
  
  print(annotate_figure(final_plot, top = text_grob(main_title, face = "bold", size = 36)))
  
  dev.off()
  
}


hazard_ratio_df <- bind_rows(hazard_ratio_df)


# saveRDS(hazard_ratio_df, "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/survival_analysis_hazard_ratios/pan_cancer_hazard_ratios.rds")


hazard_ratio_df <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/survival_analysis_hazard_ratios/pan_cancer_hazard_ratios.rds")

df <- hazard_ratio_df


df$hazard_ratio <- log(df$hazard_ratio)
df$hazard_ratio.high <- log(df$hazard_ratio.high)
df$hazard_ratio.low <- log(df$hazard_ratio.low)


df$project_signature <- paste0(df$project,df$signature)



p <- ggplot(df, aes(x=hazard_ratio,color=signature,y=project))+
  geom_point(position = position_dodge(width = 1))+
  geom_errorbar(aes(xmin=hazard_ratio.low, xmax=hazard_ratio.high), width=1, position = "dodge") +
  geom_vline(xintercept = 0, linetype="dashed")+
  theme_classic()+
  coord_cartesian(ylim=c(1,34),xlim=c(-20,25))+
  labs(x="Log(Hazard Ratio)",y="")
  


p

p_mid <- p +
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank())+
   scale_color_discrete(name = "")
  


# wrangle results into pre-plotting table form
res_plot <- df  |>
  # round estimates and 95% CIs to 2 decimal places for journal specifications
  mutate(across(
    c(hazard_ratio, hazard_ratio.low, hazard_ratio.high),
    ~ str_pad(
      round(.x, 2),
      width = 4,
      pad = "0",
      side = "right"
    )
  ),
  # add an "-" between HR estimate confidence intervals
  estimate_lab = paste0(hazard_ratio, " (", hazard_ratio.low, "-", hazard_ratio.high, ")")) |>
  # round p-values to two decimal places, except in cases where p < .001
  mutate(p_value = case_when(
    p_value < .001 ~ "<0.001",
    round(p_value, 2) == .05 ~ as.character(round(p_value,3)),
    p_value < .01 ~ str_pad( # if less than .01, go one more decimal place
      as.character(round(p_value, 3)),
      width = 4,
      pad = "0",
      side = "right"
    ),
    TRUE ~ str_pad( # otherwise just round to 2 decimal places and pad string so that .2 reads as 0.20
      as.character(round(p_value, 2)),
      width = 4,
      pad = "0",
      side = "right")
  )) |>
  # add a row of data that are actually column names which will be shown on the plot in the next step
  bind_rows(
    data.frame(
      project = "Project",
      estimate_lab = "Hazard Ratio (95% CI)",
      conf.low = "",
      conf.high = "",
      p_value = "p-value"
    )
  ) |>
  mutate(project_signature = fct_rev(fct_relevel(project_signature, "Model")))

glimpse(res_plot)

new_project <- c()
for(i in 1:length(res_plot$project)){
  if(i%%2==1){
    new_project <- append(new_project, res_plot$project[i])
  } else{
    new_project <- append(new_project,"")
  }
}

res_plot$project <- new_project

p_left <-
  res_plot  |>
  ggplot(aes(y = project_signature))

p_left <-
  p_left +
  geom_text(aes(x = 0, label = project), hjust = 0, fontface = "bold")

p_left <-
  p_left +
  geom_text(
    aes(x = 1, label = estimate_lab),
    hjust = 0,
    fontface = ifelse(res_plot$estimate_lab == "Hazard Ratio (95% CI)", "bold", "plain")
  )

p_left <-
  p_left +
  theme_void() +
  coord_cartesian(xlim = c(0, 4))

p_left

# right side of plot - pvalues
p_right <-
  res_plot  |>
  ggplot() +
  geom_text(
    aes(x = 0, y = project_signature, label = p_value),
    hjust = 0,
    fontface = ifelse(res_plot$p_value == "p-value", "bold", "plain")) +
  theme_void()

p_right



layout <- c(
  area(t = 0, l = 0, b = 30, r = 3), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
  area(t = 1, l = 4, b = 30, r = 9), # middle plot starts a little lower (t=1) because there's no title. starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30 units), and 6 units further over from the left plot (r=9 whereas left plot is r=3)
  area(t = 0, l = 9, b = 30, r = 11) # right most plot starts at top of page, begins where middle plot ends (l=9, and middle plot is r=9), goes to bottom of page (b=30), and extends two units wide (r=11)
)
# final plot arrangement


jpeg("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/survival_plots/pan_cancer/figure_4c.jpg",
     width=1400,
     height = 800)

p_left + p_mid + p_right + plot_layout(design = layout)

dev.off()

