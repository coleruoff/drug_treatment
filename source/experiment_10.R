library(Seurat)
library(presto)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(GSVA)
library(tidyverse)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

cell_lines <- c("A549","K562","MCF7")

RACs <- list(c(4,9,12,13,14,16,18,19),c(4,5,9,11),c(5,8,12,13,17))
names(RACs) <- cell_lines

emergent <- list(c(14:19),c(9),c(13,15,18))
names(emergent) <- cell_lines

all_data <- list()
for(curr_cell_line in cell_lines){
  
  cat(curr_cell_line,"\n")
  
  data <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/sciPlex_data/",curr_cell_line, "_processed_filtered.rds"))
  
  data <- data[,data$treatment_stage == 'post']
  
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
  data <- AddMetaData(data, metadata = ifelse(data$Cluster %in% emergent[[curr_cell_line]], "emergent", "non_emergent"), col.name = "emergent")
  
  data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  
  all_data[curr_cell_line] <- data
}
################################################################################
################################################################################

curr_cell_line <- "A549"

data <- all_data[[curr_cell_line]]

Idents(data) <- data$emergent_rac

de_results <- FindMarkers(data, ident.1 = "emergent_rac",ident.2 = "non_emergent_rac")
################################################################################
# Cell Cycle Phases Proportions
################################################################################

df <- data@meta.data
temp <- df %>% 
  count(rac,emergent_rac, Phase)

percents <- df %>% 
  count(emergent_rac)

temp <- merge(temp,percents,by="emergent_rac")

temp <- temp %>% 
  mutate("percent" = n.x/n.y)

ggplot(temp)+
  geom_col(aes(x=emergent_rac,y=percent,fill=Phase), position = "dodge")+
  ggtitle(paste0(curr_cell_line, " Cell Cycle Phase Percentages"))

################################################################################
# GSEA of DE genes
################################################################################

emergent_rac_ranks <- de_results %>% 
  arrange(desc(avg_log2FC)) %>% 
  rownames_to_column() %>% 
  dplyr::select(rowname,avg_log2FC) %>% 
  deframe()

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, human_gene_symbol)
# hallmark_uni <- unique(c(unlist(all_signatures), m_t2g$human_gene_symbol))

mp_t2g <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/ith_meta_programs_t2g.rds")

specifc_mps <- c("MP39 Metal-response","MP31 Alveolar","MP29 NPC/OPC","MP28 Oligo normal","MP27 Oligo Progenitor","MP38 Glutathione","MP41 Unassigned","MP35 Hemato-related-I","MP37 Hemato-related-II","MP32 Skin-pigmentation","MP36 IG","MP16 MES (glioma)","MP15 EMT IV")

mp_t2g <- mp_t2g %>% 
  filter(!gs_name %in% specifc_mps)

ranks <- emergent_rac_ranks

# Meta-programs
mp_gsea <- GSEA(ranks, TERM2GENE = mp_t2g)

df <- mp_gsea@result %>% 
  filter(p.adjust < 0.05) %>% 
  arrange(desc(NES)) %>% 
  dplyr::select(Description,NES,p.adjust)



mp_plot <- ggplot(df[1:10,])+
  geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES))+
  scale_fill_continuous(high="red",low="pink")+
  xlab("Meta-Program")+
  ggtitle(paste0(curr_cell_line, " non-emergent RAC Enrichment"))+
  theme(plot.title = element_text(size=30),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20,face="bold"))+
  coord_flip()


# Hallmarks
hallmark_gsea <- GSEA(ranks, TERM2GENE = m_t2g)

df <- hallmark_gsea@result %>% 
  filter(p.adjust < 0.05) %>% 
  arrange(desc(NES)) %>% 
  dplyr::select(Description,NES,p.adjust)


hallmark_plot <- ggplot(df[1:10,])+
  geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES))+
  scale_fill_continuous(high="red",low="pink")+
  xlab("Cancer Hallmark")+
  ggtitle(paste0(curr_cell_line, " emergent RAC Enrichment"))+
  theme(plot.title = element_text(size=30),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20,face="bold"))+
  coord_flip()


# GO Pathways
go_gsea <- gseGO(geneList     = ranks,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 minGSSize    = 100,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.01,
                 verbose      = FALSE,
                 keyType = "SYMBOL")


df <- go_gsea@result %>% 
  filter(p.adjust < 0.05) %>% 
  arrange(desc(NES)) %>% 
  dplyr::select(Description,NES,p.adjust)

go_plot <- ggplot(df[1:10,])+
  geom_col(aes(x=(fct_reorder(Description,NES)),y=NES, fill=NES))+
  scale_fill_continuous(high="red",low="pink")+
  xlab("GO Pathways")+
  ggtitle(paste0(curr_cell_line, " non-emergent RAC Enrichment"))+
  theme(plot.title = element_text(size=30),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20,face="bold"))+
  coord_flip()




plot_list <- list(hallmark_plot,mp_plot,go_plot)


p <- ggarrange(plotlist = plot_list,ncol=3,nrow=1)


p

################################################################################
# ORA of genes
################################################################################
cell_line_universes <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/cell_line_universes.rds")


pre_rac_genes <- de_results %>% 
  filter(p_val_adj < 0.05 & avg_log2FC > 0 & cluster == "non_emergent_rac") %>% 
  arrange(desc(avg_log2FC)) %>% 
  pull(gene)

emergent_rac_genes <- de_results %>% 
  filter(p_val_adj < 0.05 & avg_log2FC > 0 & cluster == "emergent_rac") %>% 
  arrange(desc(avg_log2FC)) %>% 
  pull(gene)



sort(emergent_rac_genes)

sort(pre_rac_genes)

curr_geneset <- emergent_rac_genes
universe_to_use <- cell_line_universes[[curr_cell_line]]

# Hallmarks Enrichment
hallmark_enrichment_results <- enricher(curr_geneset, TERM2GENE=m_t2g,
                                        universe = universe_to_use)


hallmark_plot <- barplot(hallmark_enrichment_results)


# MPs Enrichment
mp_enrichment_results <- enricher(curr_geneset, TERM2GENE=mp_t2g,
                                  universe = universe_to_use)



mp_plot <- barplot(mp_enrichment_results)

# GO Functional Enrichment
ego <- enrichGO(gene          = curr_geneset,
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = .01,
                qvalueCutoff  = .05,
                readable      = TRUE,
                keyType = "SYMBOL",
                universe = universe_to_use)


go_plot <- barplot(ego)


################################################################################

################################################################################

ggboxplot(data@meta.data, x="emergent_rac",y="proliferation_index",fill="emergent_rac")


data <- AddMetaData(data, scores, col.name ="resistance_score")


ggboxplot(data@meta.data, x="emergent_rac",y="resistance_score",fill="emergent_rac")


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

ggboxplot(data@meta.data, x="emergent_rac",y="percent.mt",fill="emergent_rac")


temp <- data@meta.data %>% 
  dplyr::count(emergent_rac,dose)


temp$dose <- factor(temp$dose)

ggplot(temp)+
  geom_col(aes(x=emergent_rac,y=n,fill=dose),position = "dodge")


scores <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/aucell_score_objects/",curr_cell_line,"_processed_filtered_kegg_endocytosis_aucell_scores.rds"))
data <- AddMetaData(data, scores,col.name = "endocytosis_score")

ggboxplot(data@meta.data, x="emergent_rac",y="endocytosis_score",fill="emergent_rac")

################################################################################
# Survival Analysis
################################################################################
curr_signature <- list("emergent"=emergent_rac_genes)
curr_signature <- list("pre"=pre_rac_genes)
metric_to_use <- "OS"
# metric_to_use <- "PFI"

hazard_ratios <- c()


plots <- list()
for(curr_cell_line in cell_lines){
  cat(curr_cell_line, "\n")
  
  # Get TCGA sample count data  
  tcga_project <- get_tcga_project(curr_cell_line)
  
  # Get read count data
  # tcga_matrix_vst <- get_read_count_data(tcga_project)
  tcga_matrix_vst <- readRDS(paste0("/data/CDSL_hannenhalli/Cole/TCGA_processed_data/",tcga_project,"_data.rds"))
  
  # Score samples
  # curr_signature <- all_signatures[grepl(curr_cell_line,names(all_signatures))]
  
  ssgsea_result <- gsva(tcga_matrix_vst,curr_signature,method='ssgsea')
  
  gene_set_scores_df <- as.data.frame(t(ssgsea_result)) %>%
    tibble::rownames_to_column("submitter_id")
  
  gene_set_scores_df$submitter_id <- sapply(gene_set_scores_df$submitter_id, function(x) gsub("\\.", "-",x))
  gene_set_scores_df$submitter_id <- sapply(gene_set_scores_df$submitter_id, function(x) substring(x, 1,12))
  
  # Get clinical data
  clinical_df <- get_clinical_data(tcga_project)
  
  #Run cox regression
  if("purity" %in% colnames(clinical_df)){
    covariates <- c("purity", "sex", "age")
  } else {
    covariates <- c("sex", "age")
  }
  
  
  cox_regression_info <- cox_regression(sample_survival_df = clinical_df,
                                        survival_data_id_col = "submitter_id",
                                        duration_str = paste0(metric_to_use, ".time"),
                                        feature_data_id_col = "submitter_id",
                                        model_covariates = covariates,
                                        status_str=metric_to_use,
                                        sample_features_df = gene_set_scores_df, 
                                        km_features_to_plot=names(curr_signature),
                                        low_high_percentiles = c(.5,.49))
  
  
  if(!all(is.na(cox_regression_info$km_df[[metric_to_use]]))){
    
    hazard_ratios <- append(hazard_ratios, cox_regression_info$regression_df$hazard_ratio)
    
    modform <- as.formula(paste0("Surv(",metric_to_use,".time, ",metric_to_use,") ~ ", names(curr_signature)))
    
    fit <- eval(substitute(survfit(modform, data = cox_regression_info$km_df), list(modform = modform)))
    
    
    if(metric_to_use == "OS"){
      plot_title <- paste0("\n",curr_cell_line, " RAC Signature Overall Survival (",tcga_project,")")
    } else {
      plot_title <- paste0("\n",curr_cell_line, " RAC Signature Progression Free Survival (",tcga_project,")")
    }
    
    
    cox_regression_info$km_df[[names(curr_signature)]] <- factor(cox_regression_info$km_df[[names(curr_signature)]], levels = c("High","Medium","Low"))
    
    p <- ggsurvplot(fit, data = cox_regression_info$km_df,
                    pval = T,
                    pval.size = 6,
                    legend.title="Expression Level:",
                    legend.labs= c("High","Low"),
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
    
    
    p <- p + ggtitle(plot_title)
    
    plots <- append(plots, list(p$plot))
    
  }
}

################################################################################
# Plot KM plots
################################################################################
figure <- ggarrange(plotlist = plots, ncol=3, common.legend = T, legend=c("right"))

p <- annotate_figure(figure, left = text_grob("Survival Probability", rot = 90, vjust = 1, size=35, face="bold"),
                     bottom = text_grob("Time", size=35, face="bold"),
                     top=text_grob(paste0("RAC Signatures Overall Survival in TCGA Samples"), size=40, face="bold"))


p

hazard_ratios












