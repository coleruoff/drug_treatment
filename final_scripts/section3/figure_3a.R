args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(readxl)
library(tidyverse)
library(ggpubr)
library(decoupleR)

# dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

#################################################################################

treatment_title <- "OVCAR8 IC50 Prexasertib 20 days"

# net <- get_collectri(organism='human', split_complexes=FALSE)
# full_tf_list <- unique(net$source)

supercluster_top_tfs <- readRDS(paste0(dataDirectory, "genesets/rac_supercluster_top_tf_list.rds"))
supercluster_bottom_tfs <- readRDS(paste0(dataDirectory, "genesets/rac_supercluster_bottom_tf_list.rds"))

# supercluster_top_tfs <- readRDS(paste0(dataDirectory, "genesets/supercluster_ranks_top_tfs.rds"))
# supercluster_bottom_tfs <- readRDS(paste0(dataDirectory, "genesets/supercluster_ranks_bottom_tfs.rds"))

final_df <- list()

crispr_ko_data <- read_xlsx(paste0(dataDirectory, "gottesman_crispr_data/CRISPR_KO_summary_negative.xlsx"), sheet = 5)

for(j in 1:2){
  crispr_ko_ranks <- crispr_ko_data %>% 
    arrange(`neg|fdr`) %>% 
    dplyr::select(id,`neg|rank`) 
  
  crispr_ko_ranks$ranks <- 1:nrow(crispr_ko_ranks)
  
  up_ranks <- crispr_ko_ranks %>% 
    filter(id %in% (supercluster_top_tfs[[j]])) %>% 
    pull(ranks)
  
  down_ranks <- crispr_ko_ranks %>%
    filter(id %in% (supercluster_bottom_tfs[[j]])) %>%
    pull(ranks)
  
  
  final_df[["rank"]] <- append(final_df[["rank"]], c(up_ranks,down_ranks))
  final_df[["top"]] <- append(final_df[["top"]], c(rep("Top",length(up_ranks)),rep("Bottom",length(down_ranks))))
  final_df[["geneset"]] <- append(final_df[["geneset"]], rep(paste0("Supercluster ", j),length(c(up_ranks,down_ranks))))
  
}


final_df <- data.frame(final_df)

final_df$top <- factor(final_df$top,levels = c("Top","Bottom"))

saveRDS(final_df, paste0(dataDirectory, "supercluster_tf_top_bottom_crispr_ranks_data.rds"))

p <- ggboxplot(final_df, x = "geneset", y = "rank",
               fill = "top", short.panel.labs = T, ncol=1)

p <-p + stat_compare_means(aes(group = top), label = "p.format")+
  scale_y_reverse()+
  scale_fill_discrete(name = "TFs")+
  ylab("Gene Rank")+
  xlab("")+
  ggtitle("Supercluster Top/Bottom TFs Ranks in CRISPR Screen", subtitle = "(OVCAR8 IC50 Prexasertib 20 days)")+
  theme(axis.text.x = element_text(angle=0, vjust = 1, hjust=.5),
        strip.text = element_text(size=20),
        legend.position = "right",
        legend.title = element_text(size=26),
        legend.text = element_text(size=20),
        plot.title = element_text(size=30))


png(paste0(plotDirectory,"figure_3a.png"),
    width=30,height=12, units = "in", res = 300)

print(p)

dev.off()

