library(readxl)
library(tidyverse)
library(ggpubr)
library(decoupleR)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/"

all_treatments <- c("RKO 1.9 uM Oxaliplatin 10 days","SW620 0.38 uM Oxaliplatin 21 days","SW620 0.5 uM Oxaliplatin 21 days","OVCAR8 IC25 Prexasertib 20 days","OVCAR8 IC50 Prexasertib 20 days","OVCAR8 IC75 Prexasertib 20 days")

# net <- get_collectri(organism='human', split_complexes=FALSE)
# full_tf_list <- unique(net$source)

supercluster_top_tfs <- readRDS(paste0(dataDirectory, "genesets/rac_supercluster_top_tf_list.rds"))
supercluster_bottom_tfs <- readRDS(paste0(dataDirectory, "genesets/rac_supercluster_bottom_tf_list.rds"))

# supercluster_top_tfs <- readRDS(paste0(dataDirectory, "genesets/supercluster_ranks_top_tfs.rds"))
# supercluster_bottom_tfs <- readRDS(paste0(dataDirectory, "genesets/supercluster_ranks_bottom_tfs.rds"))

final_df <- list()

i <- 5

crispr_ko_data <- read_xlsx(paste0(dataDirectory, "gottesman_crispr_data/CRISPR_KO_summary_negative.xlsx"), sheet = i)

for(j in 1:2){
  crispr_ko_ranks <- crispr_ko_data %>% 
    arrange(`neg|fdr`) %>% 
    select(id,`neg|rank`) 
  
  
  crispr_ko_ranks$ranks <- 1:nrow(crispr_ko_ranks)
  
  up_ranks <- crispr_ko_ranks %>% 
    filter(id %in% (supercluster_top_tfs[[j]])) %>% 
    pull(ranks)
  
  down_ranks <- crispr_ko_ranks %>%
    filter(id %in% (supercluster_bottom_tfs[[j]])) %>%
    pull(ranks)
  
  
  # rest_tfs <- full_tf_list[(!full_tf_list %in% supercluster_top_tfs[[j]])]
  # down_ranks <- crispr_ko_ranks %>% 
  #   filter(id %in% rest_tfs) %>% 
  #   pull(ranks)
  
  
  
  final_df[["rank"]] <- append(final_df[["rank"]], c(up_ranks,down_ranks))
  final_df[["top"]] <- append(final_df[["top"]], c(rep("Top",length(up_ranks)),rep("Bottom",length(down_ranks))))
  # final_df[["treatment"]] <- append(final_df[["treatment"]], rep(all_treatments[i],length(c(up_ranks,down_ranks))))
  final_df[["geneset"]] <- append(final_df[["geneset"]], rep(paste0("Supercluster ", j),length(c(up_ranks,down_ranks))))
  
}


temp <- data.frame(final_df)

temp$top <- factor(temp$top,levels = c("Top","Bottom"))



p <- ggboxplot(temp, x = "geneset", y = "rank",
               fill = "top", short.panel.labs = T, ncol=1)

p + stat_compare_means(aes(group = top), label = "p.format")+
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





