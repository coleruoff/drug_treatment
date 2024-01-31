library(readxl)
library(GSVA)
# Read in data
data <- read.table("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/mapki_melanoma/salmon.fastq.tpm.GSE75299-genes.txt")

# remove rows of all 0
data <- data[apply(data,MARGIN = 1, FUN = function(x) !all(x==0)),]

# Read in patient relapse data
relapse_data <- read_xlsx("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/mapki_melanoma/mapki-response.masterData.xlsx",
                           sheet = 3)

patient_sample_names <- relapse_data %>% 
  filter(grepl("Pt",relapse_data$ID)) %>% 
  pull(Run)

#Filter data to only the patient samples
data <- data[,colnames(data) %in% patient_sample_names]

# Read in metadata (sheet 2)
data_metadata <- read_xlsx("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/processed_data/mapki_melanoma/mapki-response.masterData.xlsx",
                           sheet = 2)


data_metadata <- merge(data_metadata,relapse_data,by="Run")

data_metadata <- data_metadata %>% 
  filter(data_metadata$Run %in% colnames(data)) %>% 
  mutate("treatment_stage" = ifelse(.$Timeline == "baseline", "pre","post"))


length(c(258,383,145,1095,300,62,50))
length(unique(data_metadata$ID))
PFS_data <- cbind(c(258,383,145,1095,300,62,50),unique(data_metadata$ID))
colnames(PFS_data) <- c("PFS","ID")

data_metadata <- merge(data_metadata,PFS_data,by='ID')

data_metadata <- data_metadata %>% 
  mutate("response" = ifelse(as.numeric(data_metadata$PFS) < 100, "short","long"))


#################################################################################
# Read in gene signatures

rac_type_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/global_rac_type_signatures.rds")
rac_type_signatures <- rac_type_signatures[grepl("type1",names(rac_type_signatures))]

rac_supercluster_consensus_signature <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/rac_supercluster_consensus_signature.rds")
type1_supercluster_signatures <- readRDS("/data/CDSL_hannenhalli/Cole/projects/drug_treatment/data/genesets/type1_supercluster_signatures.rds")
type1_supercluster_signatures <- type1_supercluster_signatures[1:2]

geneset_to_use <- rac_type_signatures
geneset_title <- "RAC Supercluster Signature"

#################################################################################
# Score samples

ssgsea_res <- gsva(as.matrix(data), geneset_to_use, method="ssgsea")

df <- data.frame(ssgsea_res) %>%
  rownames_to_column("geneset") %>% 
  pivot_longer(!geneset, names_to = "Run",values_to = "score")


df <- merge(df,data_metadata,by = "Run")


# Plot boxplots

p <- ggboxplot(df , x = "response", y = "score",
               color = "response", palette = "jco",
               facet.by = "geneset", short.panel.labs = FALSE,
               add="")

p + stat_compare_means(label = "p.format", method = "wilcox")+
  ggtitle(paste0(geneset_title," ssGSEA Scores"))+
  xlab("")+
  ylab("Score")+
  theme(legend.position="right",
        strip.text = element_text(size=20),
        title = element_text(size=30))


# dev.off()
