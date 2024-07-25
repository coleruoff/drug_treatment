args = commandArgs(trailingOnly=TRUE)
dataDirectory <- paste0(args[1],"final_data/")
plotDirectory <- paste0(args[1],"final_figures/")
setwd(args[1])
library(readxl)
library(tidyverse)
library(ggpubr)
library(decoupleR)
library(poolr)
library(rstatix) 
set.seed(42)

dataDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_data/"
plotDirectory <- "/data/CDSL_hannenhalli/Cole/projects/drug_treatment/final_figures/"

#################################################################################


# net <- get_collectri(organism='human', split_complexes=FALSE)
# full_tf_list <- unique(net$source)

supercluster_top_tfs <- readRDS(paste0(dataDirectory, "genesets/rac_supercluster_top_tf_list.rds"))
supercluster_bottom_tfs <- readRDS(paste0(dataDirectory, "genesets/rac_supercluster_bottom_tf_list.rds"))

# supercluster_top_tfs <- readRDS(paste0(dataDirectory, "genesets/supercluster_ranks_top_tfs.rds"))
# supercluster_bottom_tfs <- readRDS(paste0(dataDirectory, "genesets/supercluster_ranks_bottom_tfs.rds"))

final_df <- list()

crispr_ko_data <- read_xlsx("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/gottesman_crispr_data/CRISPR_KO_summary_negative.xlsx", sheet = 5)

crispr_ko_data[,5]

for(j in 1:2){
  crispr_ko_ranks <- crispr_ko_data %>% 
    arrange(`neg|score`) %>% 
    dplyr::select(id,`neg|score`) %>% 
    mutate(ranks = 1:nrow(.))  
  
  
  top_ranks <- crispr_ko_ranks %>%
    filter(id %in% supercluster_top_tfs[[j]]) %>%
    pull(`neg|score`)
  
  rest_ranks <- crispr_ko_ranks %>%
    filter(id %in% supercluster_bottom_tfs[[j]]) %>%
    pull(`neg|score`)

  # rest_ranks <- crispr_ko_ranks %>%
  #   filter(id %in% full_tf_list & !id %in% supercluster_top_tfs[[j]]) %>%
  #   pull(`neg|score`)
  
  # rest_ranks <- crispr_ko_ranks %>%
  #   filter(!id %in% supercluster_top_tfs[[j]]) %>%
  #   pull(`neg|score`)
  
  
  final_df[["rank"]] <- append(final_df[["rank"]], c(top_ranks,rest_ranks))
  final_df[["top"]] <- append(final_df[["top"]], c(rep("Top",length(top_ranks)),rep("Bottom",length(rest_ranks))))
  final_df[["geneset"]] <- append(final_df[["geneset"]], rep(paste0("Supercluster ", j),length(c(top_ranks,rest_ranks))))
  
}


final_df <- data.frame(final_df)

final_df$top <- factor(final_df$top,levels = c("Top","Bottom"))

# saveRDS(final_df, paste0(dataDirectory, "supercluster_tf_top_bottom_crispr_ranks_data.rds"))
# 


p <- ggviolin(final_df, x = "top", y = "rank",
               fill = "top", short.panel.labs = T, ncol=1,
              add="boxplot")+
  facet_wrap(~geneset, strip.position = "bottom")


stat.test <- final_df %>%
  group_by(geneset) %>%
  wilcox_test(rank ~ top, alternative = "less") %>%
  adjust_pvalue() %>%
  add_significance()

# stat.test$p.adj <- paste0("p = ", sprintf("%.2f", stat.test$p.adj))
stat.test$p.adj[2] <- "ns"


p <- p + stat_pvalue_manual(stat.test, label = "p.adj", y.position = 0, size=2, vjust = -2, label.y = .5)+
  scale_y_continuous(breaks=c(-.5,0,.5,1))+
  scale_y_reverse()+
  scale_fill_discrete(name = "TFs")+
  ylab("CRISPR Score")+
  xlab("")+
  theme(strip.text = element_text(size=10),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid"),
        legend.position = "right",
        legend.title = element_text(size=10),
        legend.text = element_text(size=8),
        axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


p

tiff(paste0(plotDirectory,"figure_3.tiff"), width=100, height = 80, units = "mm", res = 1000)

print(p)

dev.off()

# png(paste0(plotDirectory,"figure_3a.png"),
#     width=10,height=6, units = "in", res = 300)
# 
# print(p)
# 
# dev.off()





