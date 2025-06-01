args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])

source("final_scripts/drug_treatment_functions.R")
set.seed(42)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
# plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/"

################################################################################
all_data <- readRDS(paste0(dataDirectory, "processed_data/sciplex_data/all_cell_lines_data.rds"))

cell_lines <- c("A549","K562","MCF7")

supercluster_components <- readRDS(paste0(dataDirectory, "processed_data/supercluster_components.rds"))

organism_to_use <- "ecoli"

df <- list()

for(curr_cell_line in cell_lines){
  data <- all_data[[curr_cell_line]]
  
  data <- AddMetaData(data, ifelse(data$Cluster == supercluster_components[[1]][[curr_cell_line]], 1,0), col.name = "supercluster1")
  data <- AddMetaData(data, ifelse(data$Cluster == supercluster_components[[2]][[curr_cell_line]], 1,0), col.name = "supercluster2")
  data <- AddMetaData(data, ifelse(data$Cluster == supercluster_components[[3]][[curr_cell_line]], 1,0), col.name = "supercluster3")
  
  scores <- readRDS(paste0(dataDirectory, "aucell_score_objects/",curr_cell_line,"_processed_filtered_",organism_to_use,"_human_orthologs_up_aucell_scores.rds"))
  data <- AddMetaData(data, scores[,1], col.name = colnames(scores))
  
  # Supercluster 1
  curr_scores <- data@meta.data %>% 
    filter(supercluster1 == 1 & treatment_stage == "post") %>% 
    pull(paste0(organism_to_use,"_human_orthologs_up"))
  
  df[["scores"]] <- append(df[["scores"]],curr_scores)
  df[["sc"]] <- append(df[["sc"]],rep("Supercluster 1", length(curr_scores)))
  df[["cell_line"]] <- append(df[["cell_line"]],rep(curr_cell_line, length(curr_scores)))
  
  # Supercluster 2
  curr_scores <- data@meta.data %>% 
    filter(supercluster2 == 1 & treatment_stage == "post") %>% 
    pull(paste0(organism_to_use,"_human_orthologs_up"))
  
  df[["scores"]] <- append(df[["scores"]],curr_scores)
  df[["sc"]] <- append(df[["sc"]],rep("Supercluster 2", length(curr_scores)))
  df[["cell_line"]] <- append(df[["cell_line"]],rep(curr_cell_line, length(curr_scores)))
  
  # Supercluster 3
  curr_scores <- data@meta.data %>%
    filter(supercluster3 == 1 & treatment_stage == "post") %>%
    pull(paste0(organism_to_use,"_human_orthologs_up"))

  df[["scores"]] <- append(df[["scores"]],curr_scores)
  df[["sc"]] <- append(df[["sc"]],rep("Supercluster 3", length(curr_scores)))
  df[["cell_line"]] <- append(df[["cell_line"]],rep(curr_cell_line, length(curr_scores)))
  
  # Non-RACs
  curr_scores <- data@meta.data %>%
    filter(supercluster1 == 0 & supercluster2 == 0 & supercluster3 == 0 & treatment_stage == "post") %>%
    pull(paste0(organism_to_use,"_human_orthologs_up"))
  
  df[["scores"]] <- append(df[["scores"]],curr_scores)
  df[["sc"]] <- append(df[["sc"]],rep("Non-RAC", length(curr_scores)))
  df[["cell_line"]] <- append(df[["cell_line"]],rep(curr_cell_line, length(curr_scores)))
  
}

df <- data.frame(df)
df$rac <- ifelse(df$sc == "Non-RAC", "Non-RAC", "RAC")

saveRDS(df, paste0(dataDirectory, "figure_data/ecoli_boxplot_data.rds"))

my_comparisons <- list(c("Supercluster 1","Non-RAC"),c("Supercluster 2","Non-RAC"),c("Supercluster 3","Non-RAC"))


p <- ggboxplot(df, x="sc",y="scores",fill="sc", outlier.shape = NA, size = .2)+
  stat_compare_means(comparisons = my_comparisons, method.args = list(alternative="greater"), 
                     label = "p.format", label.y = c(.0675,.075,.085), size=2)+
  xlab("")+
  ylab("AUCell Score")+
  ylim(0,.095)+
  theme(legend.position="right",
        axis.text = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.title = element_text(size=6),
        legend.text = element_text(size=6),
        legend.title = element_text(size=6),
        legend.key.height = unit(1.5,"mm"),
        legend.key.width = unit(1.5,"mm"),
        axis.line = element_line(linewidth=.2),
        axis.ticks = element_line(linewidth = .2))+
  NoLegend()

p

jpeg(paste0(plotDirectory,"figure_4c.jpg"), width=80, height = 60, units = "mm", res = 600)
print(p)
dev.off()

################################################################################
# Calculate Effect Size
################################################################################
sc2_scores <- df %>% 
  filter(sc == "Supercluster 2") %>% 
  pull(scores)

sc3_scores <- df %>% 
  filter(sc == "Supercluster 3") %>% 
  pull(scores)

nr_scores <- df %>% 
  filter(sc == "Non-RAC") %>% 
  pull(scores)


sc2_d <- cohen.d(sc2_scores, nr_scores)
cat("Supercluster 2: ", sc2_d$estimate)

sc3_d <- cohen.d(sc3_scores, nr_scores)
cat("Supercluster 3: ", sc3_d$estimate)
