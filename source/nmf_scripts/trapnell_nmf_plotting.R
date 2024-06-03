library(ggpubr)
library(devtools)
library(clusterProfiler)

nmf_list <- readRDS("/data/CDSL_hannenhalli/Cole/nmf_files/nmf_list.rds")
nmf_list <- readRDS("//hpcdrive.nih.gov/CDSL_hannenhalli/Cole/nmf_files/nmf_list.rds")

combined_factor_list <- readRDS("/data/CDSL_hannenhalli/Cole/nmf_files/combined_factor_list.rds")
combined_factor_list <- readRDS("//hpcdrive.nih.gov/CDSL_hannenhalli/Cole/nmf_files/combined_factor_list.rds")

tumor_ucell_obj <- readRDS("/data/CDSL_hannenhalli/Cole/nmf_files/tumor_ucell_obj_combined.rds")
tumor_ucell_obj <- readRDS("//hpcdrive.nih.gov/CDSL_hannenhalli/Cole/nmf_files/tumor_ucell_obj_combined.rds")

merged_factor_gene_list <- readRDS("/data/CDSL_hannenhalli/Cole/nmf_files/merged_factor_list.rds")
merged_factor_gene_list <- readRDS("//hpcdrive.nih.gov/CDSL_hannenhalli/Cole/nmf_files/merged_factor_list.rds")

tumor_ucell_obj <- readRDS("/data/CDSL_hannenhalli/Cole/nmf_files/tumor_ucell_obj_merged.rds")
tumor_ucell_obj <- readRDS("//hpcdrive.nih.gov/CDSL_hannenhalli/Cole/nmf_files/tumor_ucell_obj_merged.rds")

df <- readRDS("/data/CDSL_hannenhalli/Cole/nmf_files/last_plot_df.rds")
df <- readRDS("//hpcdrive.nih.gov/CDSL_hannenhalli/Cole/nmf_files/last_plot_df.rds")



df <- FetchData(tumor_ucell_obj,vars=c(names(merged_factor_gene_list),"sample_type","cell_type")) %>%
  tidyr::pivot_longer(cols=names(merged_factor_gene_list),names_to="factor",values_to="activity")

options(repr.plot.width=12,repr.plot.height=12)

ggplot( df ) + 
  geom_boxplot(aes(x=sample_type,y=activity,fill=cell_type)) +
  facet_wrap(~factor,nrow=3) + theme_pubr(base_size=15) +
  theme(panel.grid.major=element_line(color="gray",linetype="dashed"))

options(repr.plot.width=7,repr.plot.height=7)


merged_factor_gene_list$Merged_Factor_2


"MALAT1.1" %in% rownames(data)


grep("\\.", rownames(data), value = T)


nrow(data)


genes_to_use <- merged_factor_gene_list[["Merged_Factor_1"]]
paste(genes_to_use,collapse=" ")

gene_list <- genes_to_keep
  
go_enrich <- enrichGO(gene = genes_to_use,
      
                      OrgDb = "org.Hs.eg.db",
                      keyType = 'SYMBOL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)



dotplot(go_enrich)



